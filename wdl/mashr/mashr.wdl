version 1.0

workflow mashr {
    input {
        Array[String] chromosomes
        String output_prefix
        File cell_type_list
        String acat_pattern
        String acat_q_col
        String zfile_pattern
        String susie_snp_pattern
        File gex_features
        File indep_variants
        Int n_strong_fallback = 5
        Int n_null = 5
        Int n_random = 5
        Float min_strong_observed_frac = 0.0
        Boolean debug = false
        String docker_suite
        String docker_mashr
        String zones
    }

    Array[String] cell_types = read_lines(cell_type_list)

    scatter (chrom in chromosomes) {
        call preprocess {
            input:
                prefix = "~{output_prefix}.~{chrom}",
                chrom = chrom,
                cell_types = cell_types,
                acat_pattern = acat_pattern,
                acat_q_col = acat_q_col,
                zfile_pattern = zfile_pattern,
                susie_snp_pattern = susie_snp_pattern,
                gex_features = gex_features,
                debug = debug,
                docker = docker_suite,
                zones = zones
        }

        scatter (i in range(length(preprocess.out_chunks))) {
            call munge_input {
                input:
                    cell_types = cell_types,
                    input_txt = preprocess.out_chunks[i],
                    zfiles = read_lines(preprocess.out_zfiles[i]),
                    susie_snps = read_lines(preprocess.out_susie_snps[i]),
                    indep_variants = indep_variants,
                    n_strong_fallback = n_strong_fallback,
                    n_null = n_null,
                    n_random = n_random,
                    min_strong_observed_frac = min_strong_observed_frac,
                    docker = docker_suite,
                    zones = zones
            }
        }
    }

    call combine_inputs {
        input:
            prefix = output_prefix,
            strong = flatten(munge_input.out_strong),
            strong_b = flatten(munge_input.out_strong_b),
            strong_s = flatten(munge_input.out_strong_s),
            null = flatten(munge_input.out_null),
            null_b = flatten(munge_input.out_null_b),
            null_s = flatten(munge_input.out_null_s),
            random = flatten(munge_input.out_random),
            random_b = flatten(munge_input.out_random_b),
            random_s = flatten(munge_input.out_random_s),
            docker = docker_suite,
            zones = zones
    }

    call run_mashr {
        input:
            prefix = output_prefix,
            strong_b = combine_inputs.out_strong_b,
            strong_s = combine_inputs.out_strong_s,
            null_b = combine_inputs.out_null_b,
            null_s = combine_inputs.out_null_s,
            random_b = combine_inputs.out_random_b,
            random_s = combine_inputs.out_random_s,
            docker = docker_mashr,
            zones = zones
    }
}

task preprocess {
    input {
        String prefix
        String chrom
        Array[String] cell_types
        String acat_pattern
        String acat_q_col
        String zfile_pattern
        String susie_snp_pattern
        File gex_features
        Boolean debug = false
        Int chunk_size = 100
        String docker
        String zones
    }

    command <<<
        set -e

        cat << "__EOF__" > script.py
        #!/usr/bin/env python3
        import argparse
        import itertools
        import math
        import numpy as np
        import pandas as pd
        import re
        from google.cloud import storage

        cell_types = "~{sep=',' cell_types}".split(",")
        sig_genes = []
        for cell_type in cell_types:
            print(f"Processing cell type: {cell_type} for ACAT")
            acat_path = "~{acat_pattern}".replace("{CELL_TYPE}", cell_type)
            acat_df = pd.read_csv(acat_path, sep="\t")
            sig_genes += acat_df[acat_df["~{acat_q_col}"] < 0.05]["phenotype_id"].tolist()

        sig_genes = np.unique(sig_genes)
        print(f"Found {len(sig_genes)} unique significant genes")

        df = pd.DataFrame(list(itertools.product(sig_genes, cell_types)), columns=["phenotype_id", "cell_type"])

        if df.phenotype_id[0].startswith("ENSG"):
            features = pd.read_csv(
                "~{gex_features}",
                names=["phenotype_id", "symbol", "type", "chrom", "start", "end", "gene_type"],
                header=None,
                sep="\t",
            )
            df = df.merge(features, on="phenotype_id", how="left")
        else:
            # chr1-12345-67890
            df["chrom"] = df.phenotype_id.str.split("-", expand=True).iloc[:, 0]

        df = df.loc[df.chrom == "~{chrom}", :]
        print(f"Found {df.shape[0]} gene-cell type pairs for chromosome ~{chrom}")

        df["zfile"] = df.apply(lambda x: re.sub(r"\{GENE_ID\}", x["phenotype_id"], re.sub(r"\{CELL_TYPE\}", x["cell_type"], re.sub(r"\{CHR\}", x["chrom"], "~{zfile_pattern}"))), axis=1)
        df["susie_snp"] = df.apply(lambda x: re.sub(r"\{GENE_ID\}", x["phenotype_id"], re.sub(r"\{CELL_TYPE\}", x["cell_type"], re.sub(r"\{CHR\}", x["chrom"], "~{susie_snp_pattern}"))), axis=1)

        storage_client = storage.Client()
        zfile_bucket = storage_client.get_bucket("~{zfile_pattern}".split("/")[2])
        susie_snp_bucket = storage_client.get_bucket("~{susie_snp_pattern}".split("/")[2])
        zfile_exists = [zfile_bucket.blob(x.split("/", 3)[-1]).exists() for x in df.zfile]
        susie_snp_exists = [susie_snp_bucket.blob(x.split("/", 3)[-1]).exists() for x in df.susie_snp]

        # Replace susie_snp with "NA" where files don't exist
        df["susie_snp"] = df["susie_snp"].where(susie_snp_exists, pd.NA)
        df = df.loc[zfile_exists, :]
        df[["phenotype_id", "cell_type", "zfile", "susie_snp"]].to_csv("input.tsv", sep="\t", index=False, header=False)

        # chunk params - ensure genes don't split across chunks
        chunk_size = ~{chunk_size}
        unique_genes = df["phenotype_id"].unique()
        total_genes = len(unique_genes)
        n_chunks = math.ceil(total_genes / chunk_size)
        genes_per_chunk = math.ceil(total_genes / n_chunks)

        if ~{true='True' false='False' debug}:
            n_chunks = 1

        for i in range(n_chunks):
            start_gene_idx = i * genes_per_chunk
            end_gene_idx = min(start_gene_idx + genes_per_chunk, total_genes)
            chunk_genes = unique_genes[start_gene_idx:end_gene_idx]

            # Filter dataframe to include only genes in this chunk
            chunk_df = df[df["phenotype_id"].isin(chunk_genes)]

            chunk_df.to_csv(f"~{prefix}.input.{i}.tsv", sep="\t", index=False, na_rep="NA")
            chunk_df["zfile"].to_csv(f"~{prefix}.zfile.{i}.tsv", sep="\t", index=False, header=False)
            chunk_df["susie_snp"].dropna().to_csv(f"~{prefix}.susie_snp.{i}.tsv", sep="\t", index=False, header=False)
        __EOF__

        python3 script.py && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        Array[File] out_chunks = glob("~{prefix}.input.*.tsv")
        Array[File] out_zfiles = glob("~{prefix}.zfile.*.tsv")
        Array[File] out_susie_snps = glob("~{prefix}.susie_snp.*.tsv")
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

task munge_input {
    input {
        Array[String] cell_types
        File input_txt
        Array[File] zfiles
        Array[File] susie_snps
        File indep_variants
        String prefix = basename(input_txt, ".tsv")
        Int n_strong_fallback = 5
        Int n_null = 5
        Int n_random = 5
        Float strong_pip_threshold = 0.1
        Float strong_z_threshold = 4.891638 # sqrt(qchisq(1e-6, 1, lower.tail = F))
        Float null_z_threshold = 2.0
        Float min_strong_observed_frac = 0.0
        Int seed = 42
        String docker
        String zones
    }

    command <<<
        set -e
        ulimit -s 65536

        mkdir -p zfile susie_snp && \
        ln -s ~{sep=' ' zfiles} zfile/ && \
        ln -s ~{sep=' ' susie_snps} susie_snp/

        cat << "__EOF__" > script.R
        library(dplyr)

        cell_types <- unlist(stringr::str_split("~{sep=',' cell_types}", ","))
        df <- data.table::fread("~{input_txt}", data.table = FALSE)
        indep_variants <- data.table::fread("~{indep_variants}", data.table = FALSE, header = FALSE)[, 1]

        df.z <- purrr::map_dfr(seq_len(nrow(df)), function(i) {
            path = file.path("zfile", basename(df$zfile[i]))
            data.table::fread(path, data.table = FALSE) %>%
                dplyr::transmute(
                    phenotype_id = df$phenotype_id[i],
                    cell_type = df$cell_type[i],
                    variant_id = rsid,
                    beta = beta,
                    se = se
                )
        })
        df.susie <- purrr::map_dfr(seq_len(nrow(df)), function(i) {
            if (is.na(df$susie_snp[i])) {
                return(NULL)
            }
            path <- file.path("susie_snp", basename(df$susie_snp[i]))
            data.table::fread(path, data.table = FALSE) %>%
                dplyr::filter(cs > 0 & prob > ~{strong_pip_threshold} & low_purity != 1) %>%
                dplyr::transmute(
                    phenotype_id = df$phenotype_id[i],
                    cell_type = df$cell_type[i],
                    variant_id = rsid,
                    cs = cs,
                    prob = prob
                )
        })

        mash.input <-
            dplyr::left_join(df.z, df.susie, by = c("phenotype_id", "cell_type", "variant_id")) %>%
            dplyr::mutate(
                id = stringr::str_c(phenotype_id, variant_id, sep = ":"),
                chisq = (beta / se) ** 2,
                independent = variant_id %in% indep_variants
            ) %>%
            dplyr::group_split(phenotype_id) %>%
            purrr::map(function(data) {
                set.seed(~{seed})

                data.chisq <-
                    dplyr::filter(data, independent) %>%
                    dplyr::group_by(variant_id) %>%
                    dplyr::summarize(
                        max_chisq = max(chisq, na.rm = TRUE),
                        min_chisq = min(chisq, na.rm = TRUE)
                    ) %>%
                    dplyr::ungroup()

                if (any(!is.na(data$cs))) {
                    strong.variants = unique(data$variant_id[!is.na(data$cs)])
                } else {
                    strong.variants <-
                        dplyr::filter(data.chisq, max_chisq > (~{strong_z_threshold} ** 2)) %>%
                        dplyr::arrange(dplyr::desc(max_chisq)) %>%
                        dplyr::slice_head(n = ~{n_strong_fallback}) %>%
                        dplyr::pull(variant_id)
                }
                # fallback for no strong variants
                if (length(strong.variants) == 0) {
                    strong.variants <- dplyr::filter(data, chisq == max(chisq, na.rm = TRUE)) %>%
                        dplyr::slice_head(n = 1) %>%
                        dplyr::pull(variant_id)
                }

                null.variants = dplyr::filter(data.chisq, min_chisq < (~{null_z_threshold} ** 2)) %>%
                    dplyr::slice_sample(n = ~{n_null}) %>%
                    dplyr::pull(variant_id)

                random.variants = dplyr::slice_sample(data.chisq, n = ~{n_random}) %>%
                    dplyr::pull(variant_id)

                return(list(strong = dplyr::filter(data, variant_id %in% strong.variants),
                            null = dplyr::filter(data, variant_id %in% null.variants),
                            random = dplyr::filter(data, variant_id %in% random.variants)))
            })

        # Extract and concatenate each type of data from the list of lists
        df.strong <- purrr::map_dfr(mash.input, ~ .x$strong)
        df.null <- purrr::map_dfr(mash.input, ~ .x$null)
        df.random <- purrr::map_dfr(mash.input, ~ .x$random)

        # Function to process beta and se matrices with missing rate filtering
        process_matrices <- function(data, prefix, min_observed_frac = 0) {
            # ensure all cell types are represented
            data <- dplyr::bind_rows(data, data.frame(id = "__DUMMY__", cell_type = cell_types))

            # Beta matrix
            beta_matrix <- tidyr::pivot_wider(data, id_cols = id, names_from = cell_type, values_from = beta)
            beta_matrix <- beta_matrix[, c("id", cell_types)] %>%
                dplyr::filter(id != "__DUMMY__")

            # Filter out rows that have only one non-missing entry
            non_missing_count <- rowSums(!is.na(beta_matrix[, cell_types]))
            n_cell_types <- length(cell_types)
            beta_matrix <- beta_matrix[non_missing_count > 1 & non_missing_count >= (min_observed_frac * n_cell_types), ]
            beta_matrix[is.na(beta_matrix)] <- 0

            # SE matrix - keep same rows as beta matrix
            se_matrix <- tidyr::pivot_wider(data, id_cols = id, names_from = cell_type, values_from = se)
            se_matrix <- se_matrix[, c("id", cell_types)] %>%
                dplyr::filter(id != "__DUMMY__")
            se_matrix <- se_matrix[se_matrix$id %in% beta_matrix$id, ]
            se_matrix[is.na(se_matrix)] <- 0

            # Save matrices
            data.table::fwrite(data, paste0("~{prefix}.", prefix, ".txt.gz"), sep = "\t")
            data.table::fwrite(beta_matrix, paste0("~{prefix}.", prefix, ".b.txt"), sep = "\t")
            data.table::fwrite(se_matrix, paste0("~{prefix}.", prefix, ".s.txt"), sep = "\t")
        }

        # Apply function to each dataset
        process_matrices(df.strong, "strong", min_observed_frac = ~{min_strong_observed_frac})
        process_matrices(df.null, "null")
        process_matrices(df.random, "random")
        __EOF__

        Rscript script.R && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_strong = "~{prefix}.strong.txt.gz"
        File out_null = "~{prefix}.null.txt.gz"
        File out_random = "~{prefix}.random.txt.gz"
        File out_strong_b = "~{prefix}.strong.b.txt"
        File out_strong_s = "~{prefix}.strong.s.txt"
        File out_null_b = "~{prefix}.null.b.txt"
        File out_null_s = "~{prefix}.null.s.txt"
        File out_random_b = "~{prefix}.random.b.txt"
        File out_random_s = "~{prefix}.random.s.txt"
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

task combine_inputs {
    input {
        String prefix
        Array[File] strong
        Array[File] strong_b
        Array[File] strong_s
        Array[File] null
        Array[File] null_b
        Array[File] null_s
        Array[File] random
        Array[File] random_b
        Array[File] random_s
        String docker
        String zones
    }

    command <<<
        set -e
        ulimit -s 65536

        zcat ~{sep=" " strong} | awk 'NR == 1 || $1 != "phenotype_id"' | bgzip -c > ~{prefix}.strong.txt.gz && \
        zcat ~{sep=" " null} | awk 'NR == 1 || $1 != "phenotype_id"' | bgzip -c > ~{prefix}.null.txt.gz && \
        zcat ~{sep=" " random} | awk 'NR == 1 || $1 != "phenotype_id"' | bgzip -c > ~{prefix}.random.txt.gz && \
        awk "NR == 1 || FNR > 1" ~{sep=" " strong_b} | bgzip -c > ~{prefix}.strong.b.txt.gz && \
        awk "NR == 1 || FNR > 1" ~{sep=" " strong_s} | bgzip -c > ~{prefix}.strong.s.txt.gz && \
        awk "NR == 1 || FNR > 1" ~{sep=" " null_b} | bgzip -c > ~{prefix}.null.b.txt.gz && \
        awk "NR == 1 || FNR > 1" ~{sep=" " null_s} | bgzip -c > ~{prefix}.null.s.txt.gz && \
        awk "NR == 1 || FNR > 1" ~{sep=" " random_b} | bgzip -c > ~{prefix}.random.b.txt.gz && \
        awk "NR == 1 || FNR > 1" ~{sep=" " random_s} | bgzip -c > ~{prefix}.random.s.txt.gz && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_strong = "~{prefix}.strong.txt.gz"
        File out_null = "~{prefix}.null.txt.gz"
        File out_random = "~{prefix}.random.txt.gz"
        File out_strong_b = "~{prefix}.strong.b.txt.gz"
        File out_strong_s = "~{prefix}.strong.s.txt.gz"
        File out_null_b = "~{prefix}.null.b.txt.gz"
        File out_null_s = "~{prefix}.null.s.txt.gz"
        File out_random_b = "~{prefix}.random.b.txt.gz"
        File out_random_s = "~{prefix}.random.s.txt.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "8 GB"
        disks: "local-disk 20 HDD"
        zones: zones
        preemptible: 2
    }
}

task run_mashr {
    input {
        String prefix
        File strong_b
        File strong_s
        File null_b
        File null_s
        File random_b
        File random_s
        Int posterior_samples = 0
        String docker
        String zones
    }

    command <<<
        set -e

        cat << "__EOF__" > script.R
        library(dplyr)

        read_matrix <- function(path, id = NULL) {
            mat = data.table::fread(path, data.table = FALSE) %>%
                tibble::column_to_rownames(var = "id") %>%
                as.matrix()
            if (!is.null(id)) {
                mat <- mat[id, , drop = FALSE]
            }
            return(mat)
        }
        strong.b <- read_matrix("~{strong_b}")
        strong.s <- read_matrix("~{strong_s}", id = rownames(strong.b))
        print(sprintf("Loaded %s strong variants", nrow(strong.b)))

        null.b <- read_matrix("~{null_b}")
        null.s <- read_matrix("~{null_s}", id = rownames(null.b))
        print(sprintf("Loaded %s null variants", nrow(null.b)))

        random.b <- read_matrix("~{random_b}")
        random.s <- read_matrix("~{random_s}",  id = rownames(random.b))
        print(sprintf("Loaded %s random variants", nrow(random.b)))

        vhat <- mashr::estimate_null_correlation_simple(
            mashr::mash_set_data(
                Bhat = null.b,
                Shat = null.s,
                zero_Bhat_Shat_reset = 1000
            )
        )
        data.strong <- mashr::mash_set_data(
            Bhat = strong.b,
            Shat = strong.s,
            V = vhat,
            zero_Bhat_Shat_reset = 1000
        )
        data.random <- mashr::mash_set_data(
            Bhat = random.b,
            Shat = random.s,
            V = vhat,
            zero_Bhat_Shat_reset = 1000
        )

        # U.flash <- mashr::cov_flash(data.strong, factors = "default")
        U.flash_nonneg <- mashr::cov_flash(data.strong, factors = "nonneg", tag="non_neg")
        U.pca <- mashr::cov_pca(data.strong, npc = 5)
        U.ed <- mashr::cov_ed(data.strong, c(U.flash_nonneg, U.pca))
        U.can <- mashr::cov_canonical(data.random)
        print("Computed prior covariance matrices")

        model.random <- mashr::mash(
            data = data.random,
            Ulist = c(U.ed, U.can),
            outputlevel = 1
        )
        print("Fitted random model")

        model.strong = mashr::mash(
            data.strong,
            g = ashr::get_fitted_g(model.random),
            fixg = TRUE,
            posterior_samples = ~{posterior_samples}
        )
        print("Fitted strong model")

        df.lfsr = ashr::get_lfsr(model.strong) %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "id")

        fastSave::saveRDS.pigz(model.random, "~{prefix}.mashr.model.random.rds")
        fastSave::saveRDS.pigz(model.strong, "~{prefix}.mashr.model.strong.rds")
        fastSave::saveRDS.pigz(data.strong, "~{prefix}.mashr.data.strong.rds")
        data.table::fwrite(df.lfsr, "~{prefix}.mashr.lfsr.tsv", sep = "\t")
        print("Finished")
        __EOF__

        Rscript script.R && \
        bgzip "~{prefix}.mashr.lfsr.tsv" && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_model_random = "~{prefix}.mashr.model.random.rds"
        File out_model_strong = "~{prefix}.mashr.model.strong.rds"
        File out_data_strong = "~{prefix}.mashr.data.strong.rds"
        File out_df_lfsr = "~{prefix}.mashr.lfsr.tsv.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "64 GB"
        disks: "local-disk 20 HDD"
        zones: zones
        preemptible: 2
    }
}
