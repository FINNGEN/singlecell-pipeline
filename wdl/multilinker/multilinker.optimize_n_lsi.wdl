version 1.0

workflow multilinker {
    input {
        String chromosome = "chr1"
        File cell_type_list
        String gex_h5ad
        String atac_h5ad_pattern
        String all_cell_annot
        String obs_sample_id
        String covar_pattern
        Array[String] covariates
        String offset = "NULL"
        String atac_lsi_pattern
        String gex_nonzero_prop_pattern
        Array[Int] k_values = [0, 5, 10, 15, 20, 25, 30, 50]
        Float nonzero_prop_threshold = 0.05
        Int n_genes = 1000
        Int chunk_size = 100
        Int seed = 42
        Boolean debug = false
        String docker_suite
        String zones
    }

    Array[String] cell_types = if debug then [read_lines(cell_type_list)[0]] else read_lines(cell_type_list)
    String atac_h5ad = sub(atac_h5ad_pattern, "\\{CHR\\}", chromosome)

    scatter (cell_type in cell_types) {
        # Step 1: Stratified sampling by sparsity + chunking
        call sample_genes {
            input:
                cell_type = cell_type,
                all_cell_annot = all_cell_annot,
                gex_nonzero_prop = sub(gex_nonzero_prop_pattern, "\\{CELL_TYPE\\}", cell_type),
                nonzero_prop_threshold = nonzero_prop_threshold,
                n_genes = n_genes,
                chunk_size = chunk_size,
                seed = seed,
                docker = docker_suite,
                zones = zones
        }

        # Step 2: Scatter hurdle fits over gene chunks
        scatter (gene_chunk in sample_genes.out_gene_chunks) {
            call run_optimize_lsi {
                input:
                    gex_h5ad = gex_h5ad,
                    atac_h5ad = atac_h5ad,
                    cell_type = cell_type,
                    all_cell_annot = all_cell_annot,
                    obs_sample_id = obs_sample_id,
                    covar = sub(covar_pattern, "\\{CELL_TYPE\\}", cell_type),
                    covariates = covariates,
                    offset = offset,
                    atac_lsi = sub(atac_lsi_pattern, "\\{CELL_TYPE\\}", cell_type),
                    gene_list = gene_chunk,
                    k_values = k_values,
                    docker = docker_suite,
                    zones = zones
            }
        }

        # Step 3: Combine chunks and make scree plot
        call summarize_optimize_lsi {
            input:
                results = run_optimize_lsi.out_results,
                cell_type = cell_type,
                docker = docker_suite,
                zones = zones
        }
    }

    output {
        Array[File] out_results = summarize_optimize_lsi.out_results
        Array[File] out_summary = summarize_optimize_lsi.out_summary
        Array[File] out_plot = summarize_optimize_lsi.out_plot
    }
}

# ---------------------------------------------------------------------------
# Task 1: Stratified gene sampling by nonzero proportion + chunking
# ---------------------------------------------------------------------------
task sample_genes {
    input {
        String cell_type
        String all_cell_annot
        File gex_nonzero_prop
        Float nonzero_prop_threshold = 0.05
        Int n_genes = 1000
        Int chunk_size = 100
        Int seed = 42
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        cat << "__EOF__" > script.R
        library(dplyr)
        set.seed(~{seed})

        cell_type_full <- "~{cell_type}"
        cell_type <- if ("~{cell_type}" == "~{all_cell_annot}") {
            NULL
        } else {
            stringr::str_replace(cell_type_full, "^.*\\.([^\\.]+)$", "\\1")
        }

        # Read nonzero proportions
        nzp <- data.table::fread("~{gex_nonzero_prop}", data.table = FALSE) %>%
            dplyr::filter(nonzero_prop >= ~{nonzero_prop_threshold}) %>%
            dplyr::mutate(decile = cut(nonzero_prop,
                breaks = quantile(nonzero_prop, probs = seq(0, 1, 0.1)),
                include.lowest = TRUE, labels = FALSE))

        message(sprintf("Expressed genes above threshold: %d", nrow(nzp)))

        # Stratified sampling by nonzero proportion deciles
        n_genes <- min(~{n_genes}, nrow(nzp))
        n_per_decile <- ceiling(n_genes / 10)
        sampled <-
            dplyr::group_by(nzp, decile) %>%
            dplyr::mutate(.rand = sample(dplyr::n())) %>%
            dplyr::filter(.rand <= n_per_decile) %>%
            dplyr::select(-.rand) %>%
            dplyr::ungroup()
        message(sprintf("Sampled %d genes across %d deciles", nrow(sampled), length(unique(sampled$decile))))

        # Write sampled gene list
        gene_ids <- sampled$gene_id

        # Chunk into files
        n_chunks <- ceiling(length(gene_ids) / ~{chunk_size})
        for (i in seq_len(n_chunks)) {
            start <- (i - 1) * ~{chunk_size} + 1
            end <- min(i * ~{chunk_size}, length(gene_ids))
            writeLines(gene_ids[start:end], sprintf("genes.chunk.%d.txt", i - 1))
        }
        message(sprintf("Written %d chunks of up to %d genes", n_chunks, ~{chunk_size}))
        __EOF__

        Rscript script.R && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        Array[File] out_gene_chunks = glob("genes.chunk.*.txt")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "8 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}

# ---------------------------------------------------------------------------
# Task 2: Fit covariate-only hurdle models for a chunk of genes at all K
# ---------------------------------------------------------------------------
task run_optimize_lsi {
    input {
        String gex_h5ad
        String atac_h5ad
        String cell_type
        String all_cell_annot
        String obs_sample_id
        File covar
        Array[String] covariates
        String offset = "NULL"
        File atac_lsi
        File gene_list
        Array[Int] k_values
        String prefix = basename(gene_list, ".txt")
        String docker
        String zones
    }

    String offset_var = if offset != "NULL" then "\"~{offset}\"" else "NULL"

    command <<<
        set -e -o pipefail

        cat << "__EOF__" > script.R
        library(multilinkeR)
        library(fasthurdle)
        library(dplyr)
        library(data.table)

        # --- Parse inputs ---
        cell_type_full <- "~{cell_type}"
        cell_type_col <- stringr::str_replace(cell_type_full, "\\.[^\\.]+$", "")
        cell_type <- if ("~{cell_type}" == "~{all_cell_annot}") {
            NULL
        } else {
            stringr::str_replace(cell_type_full, "^.*\\.([^\\.]+)$", "\\1")
        }

        base_covariates <- unlist(stringr::str_split("~{sep=',' covariates}", ","))
        offset <- ~{offset_var}
        k_values <- as.integer(unlist(stringr::str_split("~{sep=',' k_values}", ",")))

        # Ensure offset is not duplicated in covariates
        if (!is.null(offset) && offset %in% base_covariates) {
            base_covariates <- setdiff(base_covariates, offset)
        }

        # --- Load gene list for this chunk ---
        gene_ids <- readLines("~{gene_list}")
        message(sprintf("Processing %d genes", length(gene_ids)))

        # --- Load data via multilinkeR (same pattern as run_multilinker) ---
        obj <- multilinkeR::create_multilinker(
            rna = "~{gex_h5ad}",
            atac = "~{atac_h5ad}",
            cell_type_col = cell_type_col
        )

        covar_df <- data.table::fread("~{covar}", data.table = FALSE)
        if ("IID" %in% colnames(covar_df)) {
            covar_df <- dplyr::rename(covar_df, ~{obs_sample_id} = "IID")
        }
        if (any(stringr::str_starts(colnames(covar_df), "PEER"))) {
            base_covariates <- c(base_covariates, colnames(covar_df)[stringr::str_starts(colnames(covar_df), "PEER")])
        }

        atac_lsi_df <- data.table::fread("~{atac_lsi}", data.table = FALSE)

        if (nrow(atac_lsi_df) < 30000) {
            # rare cell workaround
            k_values <- seq(0, 10)
        }

        max_k <- min(max(k_values), ncol(atac_lsi_df) - 1)
        lsi_cols_all <- paste0("ATAC_LSI_", seq_len(max_k))
        k_values <- k_values[k_values <= max_k]

        # Build cell-level metadata: ATAC QC + sample covariates + LSI
        atac_obs <- dplyr::transmute(
                obj$atac$obs,
                ~{obs_sample_id},
                log_total_counts_atac = log(total_counts),
                log_total_counts_wo_highly_expressed_atac = log(total_counts_wo_highly_expressed)
            ) %>%
            tibble::rownames_to_column("barcode") %>%
            dplyr::left_join(covar_df, by = "~{obs_sample_id}") %>%
            dplyr::left_join(atac_lsi_df, by = "barcode") %>%
            tibble::column_to_rownames("barcode")
        obj$set_metadata(atac_obs)

        # --- Get metadata ---
        metadata_cols <- unique(c(base_covariates, offset, lsi_cols_all))
        metadata <- obj$get_metadata(metadata_cols, cell_type = cell_type)

        # --- Get expression matrix for this gene chunk ---
        rna_matrix <- obj$get_rna_matrix(gene_ids, cell_type)

        # --- Fit covariate-only hurdle models ---
        message(sprintf("Testing %d genes x %d K values", length(gene_ids), length(k_values)))

        # Helper: fit covariate-only hurdle for one gene
        fit_gene <- function(gene_id, X, Z, offsetx) {
            y <- as.integer(rna_matrix[, gene_id])
            if (sum(y == 0) == 0 || sum(y > 0) < 5) return(NULL)
            tryCatch({
                model <- fast_negbin_hurdle(X, y, Z = Z, offsetx = offsetx)
                if (!model$converged) return(NULL)
                data.frame(
                    gene_id = gene_id,
                    aic = AIC(model),
                    bic = BIC(model),
                    loglik = as.numeric(logLik(model)),
                    n_cells = length(y),
                    n_expr = sum(y > 0),
                    stringsAsFactors = FALSE
                )
            }, error = function(e) NULL)
        }

        # Filter metadata to complete cases (drops NAs from covar join)
        all_cols <- unique(c(base_covariates, offset, lsi_cols_all))
        complete_rows <- complete.cases(metadata[, base_covariates, drop = FALSE])
        if (sum(!complete_rows) > 0) {
            message(sprintf("Removed %d rows with NA values in covariates", sum(!complete_rows)))
        }
        metadata <- metadata[complete_rows, , drop = FALSE]
        rna_matrix <- rna_matrix[complete_rows, , drop = FALSE]

        results_df <- purrr::map_dfr(k_values, function(K) {
            message(sprintf("  K = %d ...", K))

            # Build covariate set for this K
            if (K > 0) {
                lsi_covariates <- paste0("ATAC_LSI_", seq_len(K))
            } else {
                lsi_covariates <- NULL
            }

            # Build design matrix X once (same for all genes at this K)
            formula_str <- paste("~", paste(base_covariates, collapse = " + "))
            X <- model.matrix(as.formula(formula_str), data = metadata)

            # Build Z for zero model: start from X, optionally add offset, optionally add LSI
            offsetx <- if (!is.null(offset)) metadata[[offset]] else NULL
            Z <- X
            if (!is.null(offset)) {
                Z <- cbind(Z, metadata[[offset]])
                colnames(Z)[ncol(Z)] <- offset
            }
            if (!is.null(lsi_covariates)) {
                lsi_start <- ncol(Z) + 1
                Z <- cbind(Z, as.matrix(metadata[, lsi_covariates, drop = FALSE]))
                colnames(Z)[lsi_start:(lsi_start + length(lsi_covariates) - 1)] <- lsi_covariates
            }

            res <- purrr::map_dfr(gene_ids, fit_gene, X = X, Z = Z, offsetx = offsetx)
            if (nrow(res) > 0) res$K <- K
            message(sprintf("    converged: %d / %d genes", nrow(res), length(gene_ids)))
            res
        })

        results_df$cell_type <- cell_type_full

        write.table(results_df, "optimize_lsi.~{prefix}.results.tsv",
                    row.names = FALSE, quote = FALSE, sep = "\t")
        message("Done!")
        __EOF__

        Rscript script.R && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_results = "optimize_lsi.~{prefix}.results.tsv"
    }

    runtime {
        docker: docker
        cpu: 4
        memory: "64 GB"
        disks: "local-disk 100 HDD"
        zones: zones
        preemptible: 2
    }
}

# ---------------------------------------------------------------------------
# Task 3: Combine chunk results, compute summary, and make scree plot
# ---------------------------------------------------------------------------
task summarize_optimize_lsi {
    input {
        Array[File] results
        String cell_type
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        cat << "__EOF__" > script.R
        library(dplyr)
        library(ggplot2)
        library(data.table)

        cell_type_full <- "~{cell_type}"

        # Combine all chunk results
        result_files <- unlist(stringr::str_split("~{sep=',' results}", ","))
        results_df <- purrr::map_dfr(result_files, function(f) {data.table::fread(f, data.table = FALSE)})
        message(sprintf("Combined %d rows from %d chunks", nrow(results_df), length(result_files)))

        # Derive k_values from the data (allows run_optimize_lsi to override K grid)
        k_values <- sort(unique(results_df$K))
        message(sprintf("K values found in data: %s", paste(k_values, collapse = ", ")))

        # Write combined per-gene results
        write.table(results_df, "optimize_lsi.~{cell_type}.results.tsv",
                    row.names = FALSE, quote = FALSE, sep = "\t")

        # Compute summary: delta-BIC relative to K=0 baseline
        # Restrict to genes that converged at ALL K values
        genes_all_k <- results_df %>%
            dplyr::group_by(gene_id) %>%
            dplyr::filter(dplyr::n() == length(k_values)) %>%
            dplyr::ungroup()

        baseline <- genes_all_k %>%
            dplyr::filter(K == 0) %>%
            dplyr::select(gene_id, bic_0 = bic, aic_0 = aic, loglik_0 = loglik)

        summary_df <- genes_all_k %>%
            dplyr::left_join(baseline, by = "gene_id") %>%
            dplyr::filter(is.finite(bic) & is.finite(bic_0)) %>%
            dplyr::mutate(
                delta_bic = bic - bic_0,
                delta_aic = aic - aic_0,
                delta_loglik = loglik - loglik_0
            ) %>%
            dplyr::group_by(cell_type, K) %>%
            dplyr::summarize(
                n_genes = dplyr::n(),
                mean_bic = mean(bic),
                mean_delta_bic = mean(delta_bic),
                median_delta_bic = median(delta_bic),
                sd_delta_bic = sd(delta_bic),
                mean_aic = mean(aic),
                mean_delta_aic = mean(delta_aic),
                mean_loglik = mean(loglik),
                mean_delta_loglik = mean(delta_loglik),
                .groups = "drop"
            )

        write.table(summary_df, "optimize_lsi.~{cell_type}.summary.tsv",
                    row.names = FALSE, quote = FALSE, sep = "\t")

        # Scree plot
        p <- ggplot(summary_df, aes(x = K, y = mean_delta_bic)) +
            geom_line() +
            geom_point() +
            geom_errorbar(aes(
                ymin = mean_delta_bic - sd_delta_bic / sqrt(n_genes),
                ymax = mean_delta_bic + sd_delta_bic / sqrt(n_genes)
            ), width = 1) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
            labs(
                x = "Number of ATAC LSI components (K)",
                y = "Mean delta-BIC (vs K=0)",
                title = cell_type_full,
                subtitle = sprintf("n = %d genes", summary_df$n_genes[1])
            ) +
            theme_bw()
        ggsave("optimize_lsi.~{cell_type}.scree.pdf", p, width = 6, height = 4)

        message("Done!")
        __EOF__

        Rscript script.R && \
        bgzip optimize_lsi.~{cell_type}.results.tsv && \
        bgzip optimize_lsi.~{cell_type}.summary.tsv && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_results = "optimize_lsi.~{cell_type}.results.tsv.gz"
        File out_summary = "optimize_lsi.~{cell_type}.summary.tsv.gz"
        File out_plot = "optimize_lsi.~{cell_type}.scree.pdf"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "8 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}
