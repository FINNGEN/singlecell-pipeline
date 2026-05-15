version 1.0

task generate_peak_gene_pairs {
    input {
        String chromosome
        File rna_nonzero_prop
        File atac_nonzero_prop
        Float nonzero_prop_threshold = 0.05
        Int chunk_size = 25
        Boolean debug = false
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        cat << "__EOF__" > script.R
        library(multilinkeR)
        library(dplyr)

        nonzero_prop_threshold = ~{nonzero_prop_threshold}
        chunk_size = ~{chunk_size}

        # Read gene/peak names from nonzero_prop files (no h5ad needed)
        rna_nonzero_prop = read.table("~{rna_nonzero_prop}", header = TRUE, sep = "\t")
        atac_nonzero_prop = read.table("~{atac_nonzero_prop}", header = TRUE, sep = "\t") %>%
            dplyr::rename(peak_id = gene_id)

        # Filter by nonzero proportion threshold
        if (nonzero_prop_threshold > 0) {
            rna_nonzero_prop = dplyr::filter(rna_nonzero_prop, nonzero_prop >= nonzero_prop_threshold)
            atac_nonzero_prop = dplyr::filter(atac_nonzero_prop, nonzero_prop >= nonzero_prop_threshold)
        }

        # Filter ATAC peaks to current chromosome (peak_id format: "chr1-1000-2000")
        atac_nonzero_prop = dplyr::filter(atac_nonzero_prop, sub("-.*", "", peak_id) == "~{chromosome}")

        # Generate peak-gene pairs from gene TSS BED and peak names
        gene_info <- multilinkeR::create_gene_granges(subset_genes = rna_nonzero_prop$gene_id)
        peak_info <- multilinkeR::create_peak_granges(atac_nonzero_prop$peak_id)
        peak_gene_pairs <- multilinkeR::generate_peak_gene_pairs(gene_info, peak_info)

        message(sprintf("Generated %d peak-gene pairs for %d genes",
            nrow(peak_gene_pairs), length(unique(peak_gene_pairs$gene_id))))

        # Write full pairs file
        write.table(peak_gene_pairs, "~{chromosome}.peak_gene_pairs.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

        # Chunk by genes, distributing peak-gene pairs evenly across chunks
        gene_ids <- unique(peak_gene_pairs$gene_id)
        peaks_per_gene <- as.integer(table(peak_gene_pairs$gene_id)[gene_ids])
        n_genes <- length(gene_ids)
        total_pairs <- nrow(peak_gene_pairs)

        # In debug mode, limit to one chunk_size worth of genes
        if (~{true='TRUE' false='FALSE' debug}) {
            gene_ids <- head(gene_ids, chunk_size)
            peaks_per_gene <- peaks_per_gene[seq_along(gene_ids)]
            n_genes <- length(gene_ids)
            total_pairs <- sum(peaks_per_gene)
            peak_gene_pairs <- peak_gene_pairs[peak_gene_pairs$gene_id %in% gene_ids, ]
        }

        # Target number of chunks based on chunk_size (genes per chunk)
        n_chunks <- max(1L, ceiling(n_genes / chunk_size))
        target_pairs <- ceiling(total_pairs / n_chunks)

        # Greedy assignment: fill each chunk up to target_pairs, never split a gene
        chunk_assignments <- integer(n_genes)
        current_chunk <- 1L
        current_pairs <- 0L
        for (i in seq_len(n_genes)) {
            # Start new chunk if adding this gene would exceed target and chunk is non-empty
            if (current_pairs > 0 && current_pairs + peaks_per_gene[i] > target_pairs) {
                current_chunk <- current_chunk + 1L
                current_pairs <- 0L
            }
            chunk_assignments[i] <- current_chunk
            current_pairs <- current_pairs + peaks_per_gene[i]
            if (peaks_per_gene[i] > target_pairs) {
                warning(sprintf("gene '%s' has %d peaks (> target %d per chunk)",
                    gene_ids[i], peaks_per_gene[i], target_pairs))
            }
        }
        n_chunks <- max(chunk_assignments)

        # Write chunk files
        for (ch in seq_len(n_chunks)) {
            chunk_genes <- gene_ids[chunk_assignments == ch]
            chunk_pairs <- peak_gene_pairs[peak_gene_pairs$gene_id %in% chunk_genes, c("peak_id", "gene_id")]
            outfile <- sprintf("~{chromosome}.pairs.%d.tsv", ch - 1L)
            write.table(chunk_pairs, outfile, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
            message(sprintf("  Chunk %d: %d genes, %d pairs", ch, length(chunk_genes), nrow(chunk_pairs)))
        }

        message(sprintf("Total: %d chunks (target ~%d pairs each)", n_chunks, target_pairs))
        __EOF__

        Rscript script.R && \
        touch _SUCCESS
    >>>

    output {
        File success = "_SUCCESS"
        File out_peak_gene_pairs = "~{chromosome}.peak_gene_pairs.tsv"
        Array[File] out_chunk_peak_gene_pairs = glob("~{chromosome}.pairs.*.tsv")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}

task run_multilinker {
    input {
        String method = "open4gene"
        String rna_h5ad
        String atac_h5ad
        File peak_gene_pairs
        String cell_type
        String all_cell_annot
        String obs_sample_id
        File covar
        File? atac_lsi
        File? n_lsi
        Array[String] covariates
        String offset = "NULL"
        String test = "score"
        Boolean use_spa = false
        Boolean binarize = false
        String prefix = basename(peak_gene_pairs, ".tsv")
        String docker
        String zones
    }

    String offset_var = if offset != "NULL" then "\"~{offset}\"" else "NULL"

    command <<<
        set -e -o pipefail

        cat << "__EOF__" > script.R
        library(multilinkeR)
        library(dplyr)

        rna_h5ad_path = "~{rna_h5ad}"
        atac_h5ad_path = "~{atac_h5ad}"
        cell_type = "~{cell_type}"
        cell_type_col = stringr::str_replace(cell_type, "\\.[^\\.]+$", "")
        cell_type = if ("~{cell_type}" == "~{all_cell_annot}") {NULL} else {stringr::str_replace(cell_type, "^.*\\.([^\\.]+)$", "\\1")}
        covariates = unlist(stringr::str_split("~{sep=',' covariates}", ","))
        covariates = covariates[nchar(covariates) > 0]
        zero_covariates = NULL
        offset = ~{offset_var}
        spa_cutoff = if (~{true='TRUE' false='FALSE' use_spa}) 2 else NULL

        covar = data.table::fread("~{covar}", data.table = FALSE)
        if ("IID" %in% colnames(covar)) {
            covar <- dplyr::rename(covar, ~{obs_sample_id} = "IID")
        }
        if (any(stringr::str_starts(colnames(covar), "PEER"))) {
            covariates <- c(covariates, colnames(covar)[stringr::str_starts(colnames(covar), "PEER")])
        }

        obj <- multilinkeR::create_multilinker(
            rna = rna_h5ad_path,
            atac = atac_h5ad_path,
            cell_type_col = cell_type_col
        )

        # Add ATAC and sample-level metadata
        atac_obs <- dplyr::transmute(
                obj$atac$obs,
                ~{obs_sample_id},
                log_total_counts_atac = log(total_counts),
                log_total_counts_wo_highly_expressed_atac = log(total_counts_wo_highly_expressed)
            ) %>%
            tibble::rownames_to_column("barcode") %>%
            dplyr::left_join(covar, by = "~{obs_sample_id}")

        if (~{if defined(atac_lsi) then "TRUE" else "FALSE"}) {
            atac_lsi = data.table::fread("~{atac_lsi}", data.table = FALSE)

            if (~{if defined(n_lsi) then "TRUE" else "FALSE"}) {
                K = data.table::fread("~{n_lsi}", data.table = FALSE) %>%
                    dplyr::filter(cell_type == "~{cell_type}") %>%
                    dplyr::pull(n_lsi)
            } else {
                K = 20
            }
            if (K > 0) {
                zero_covariates = paste0("ATAC_LSI_", seq_len(K))
                atac_lsi = atac_lsi[, c("barcode", zero_covariates)]
                atac_obs = dplyr::left_join(atac_obs, atac_lsi, by = "barcode")
            }
        }
        atac_obs = tibble::column_to_rownames(atac_obs, "barcode")
        obj$set_metadata(atac_obs)

        peak_gene_pairs <- read.table("~{peak_gene_pairs}", header = FALSE, sep = "\t") %>%
            dplyr::rename(peak_id = V1, gene_id = V2)

        results <- obj$compute_links(
            method = "~{method}",
            peak_gene_pairs = peak_gene_pairs,
            cell_type = cell_type,
            covariates = covariates,
            zero_covariates = zero_covariates,
            offset = offset,
            binarize = ~{true='TRUE' false='FALSE' binarize},
            test = "~{test}",
            spa_cutoff = spa_cutoff,
            nonzero_prop_threshold = 0,
            n_cores = parallelly::availableCores()
        )
        results$cell_type <- "~{cell_type}"
        if ("~{cell_type}" == "~{all_cell_annot}") {
            results <- dplyr::select(results, peak_id, gene_id, cell_type, dplyr::everything())
        }

        write.table(results, "~{method}.~{prefix}.results.tsv", row.names = FALSE, quote = FALSE, sep = "\t")
        __EOF__

        Rscript script.R && \
        touch _SUCCESS
    >>>

    output {
        File success = "_SUCCESS"
        File out_results = "~{method}.~{prefix}.results.tsv"
    }

    runtime {
        docker: docker
        cpu: 16
        memory: "80 GB"
        disks: "local-disk 20 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task combine_results {
    input {
        Array[File] results
        String prefix
        Boolean bgzip = false
        String docker
        String zones
    }
    String outname = if bgzip then "~{prefix}.results.tsv.gz" else "~{prefix}.results.tsv"
    String out_pipe = if bgzip then "| bgzip -c > ~{outname}" else "> ~{outname}"

    command <<<
        set -e -o pipefail

        awk '
        BEGIN {
            OFS = "\t"
        }
        NR == 1 || FNR > 1 {
            print
        }' ~{sep=" " results} ~{out_pipe}
    >>>

    output {
        File out_results = outname
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
