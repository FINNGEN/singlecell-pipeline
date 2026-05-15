version 1.0

workflow munge_tensorqtl {
    input {
        File gex_features
        File tss_bed
        String cis_pattern
        String cis_nominal_pattern
        File cell_type_list
        Boolean debug = false
        String docker_suite
        String zones
    }

    Array[String] cell_types = if debug then [read_lines(cell_type_list)[0]] else read_lines(cell_type_list)

    scatter (cell_type in cell_types) {
        File cis = sub(cis_pattern, "\\{CELL_TYPE\\}", cell_type)
        File cis_nominal = sub(cis_nominal_pattern, "\\{CELL_TYPE\\}", cell_type)

        call munge {
            input:
                gex_features = gex_features,
                tss_bed = tss_bed,
                cis = cis,
                cis_nominal = cis_nominal,
                docker = docker_suite,
                zones = zones
        }
    }
}

task munge {
    input {
        File gex_features
        File tss_bed
        File cis
        File cis_nominal
        String prefix = basename(cis_nominal, ".cis_qtl_pairs.tsv.gz")
        String variant_col = "variant_id"
        String gene_col = "phenotype_id"
        String chromosome_col = "#chrom"
        String position_col = "position"
        String beta_col = "slope"
        String se_col = "slope_se"
        String pval_col = "pval_nominal"
        String freq_col = "af"
        Boolean meta_mode = false
        String docker
        String zones
    }

    command <<<
        set -e

        mv ~{cis_nominal} ~{prefix}.gz

        cat << "__EOF__" > script.R
        library(dplyr)

        df.cis = data.table::fread("~{cis}", data.table = FALSE)
        gene_ids = dplyr::filter(df.cis, qval < 0.05)[["~{gene_col}"]]

        df = data.table::fread("~{prefix}.gz", data.table = FALSE) %>%
            dplyr::filter(~{gene_col} %in% gene_ids) %>%
            dplyr::mutate(`~{variant_col}` = stringr::str_replace_all(~{variant_col}, "_", ":"))

        df.out = dplyr::transmute(
                df,
                SNP = ~{variant_col},
                gene = ~{gene_col},
                beta = ~{beta_col},
                `t-stat` = ~{beta_col} / ~{se_col},
                `p-value` = ~{pval_col},
                FDR = 0
            ) %>%
            tidyr::drop_na()

        df.esi = dplyr::transmute(
            df,
            chrom = stringr::str_remove(`~{chromosome_col}`, "^chr"),
            SNP = ~{variant_col},
            cm = 0,
            pos = ~{position_col},
            alt = stringr::str_split_fixed(~{variant_col}, "[:_]", 4)[, 4],
            ref = stringr::str_split_fixed(~{variant_col}, "[:_]", 4)[, 3],
            freq = ~{freq_col}) %>%
            dplyr::distinct()

        if (~{true='TRUE' false='FALSE' meta_mode}) {
            df.esi = dplyr::group_by(df.esi, SNP) %>%
                dplyr::mutate(freq = mean(freq)) %>%
                dplyr::ungroup() %>%
                dplyr::distinct()
        }
        stopifnot(length(unique(df.esi$SNP)) == nrow(df.esi))

        if (stringr::str_starts(gene_ids[1], "ENSG")) {
            df.genes = data.table::fread("~{gex_features}", header = FALSE, data.table = FALSE) %>%
                dplyr::transmute(gene_id = V1, symbol = V2)
            df.epi = data.table::fread("~{tss_bed}", header = FALSE, data.table = FALSE) %>%
                dplyr::transmute(
                    chromosome = stringr::str_remove(V1, "^chr"),
                    gene_id = V4,
                    cM = 0,
                    position = V2,
                ) %>%
                dplyr::left_join(df.genes, by = "gene_id") %>%
                dplyr::mutate(orientation = "+")
        } else {
            parse_peak_id = function(peak_id) {
                x = stringr::str_split_fixed(peak_id, "-", 3)
                return(tibble::tibble(
                    peak_chrom = x[, 1],
                    peak_start = as.numeric(x[, 2]),
                    peak_end = as.numeric(x[, 3])
                ))
            }
            df.epi = dplyr::distinct(df.cis, phenotype_id) %>%
                dplyr::mutate(parse_peak_id(phenotype_id)) %>%
                dplyr::transmute(
                    chromosome = stringr::str_remove(peak_chrom, "^chr"),
                    gene_id = phenotype_id,
                    cM = 0,
                    position = (peak_start + peak_end) %/% 2,
                    symbol = phenotype_id,
                    orientation = "+",
                )
        }

        data.table::fwrite(df.out, "~{prefix}.mateQTL.txt.gz", quote = F, row.names = F, sep = "\t")
        data.table::fwrite(df.esi, "~{prefix}.tmp.esi", quote = F, row.names = F, col.names = F, sep = "\t")
        data.table::fwrite(df.epi, "~{prefix}.tmp.epi", quote = F, row.names = F, col.names = F, sep = "\t")
        __EOF__

        Rscript script.R && \
        smr --eqtl-summary ~{prefix}.mateQTL.txt.gz --matrix-eqtl-format --make-besd --out ~{prefix} && \
        smr --beqtl-summary ~{prefix} --update-esi ~{prefix}.tmp.esi --update-epi ~{prefix}.tmp.epi && \
        touch _SUCCESS
    >>>

    output {
        File out_besd = "~{prefix}.besd"
        File out_esi = "~{prefix}.esi"
        File out_epi = "~{prefix}.epi"
        File success = "_SUCCESS"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "32 GB"
        disks: "local-disk 20 HDD"
        zones: zones
        preemptible: 2
    }
}