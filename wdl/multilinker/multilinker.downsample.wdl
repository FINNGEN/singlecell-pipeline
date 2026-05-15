version 1.0

import "multilinker.tasks.wdl" as tasks

workflow multilinker {
    input {
        String method = "open4gene"
        Array[String] chromosomes
        String cell_type = "predicted.celltype.l1.CD4_T"
        File downsample_condition_list
        String gex_h5ad_pattern
        String atac_h5ad_pattern
        String gex_nonzero_prop_pattern
        String atac_nonzero_prop_pattern
        String all_cell_annot
        String obs_sample_id
        File covar
        Array[String] covariates
        String offset = "NULL"
        Float nonzero_prop_threshold = 0.05
        Int min_chunk_size = 1000
        Boolean binarize = false
        Boolean debug = false
        String docker_suite
        String zones
    }

    Array[String] downsample_conditions = if debug then [read_lines(downsample_condition_list)[0]] else read_lines(downsample_condition_list)

    scatter (condition in downsample_conditions) {
        String gex_h5ad = sub(gex_h5ad_pattern, "\\{CONDITION\\}", condition)

        scatter (chrom in chromosomes) {
            String atac_h5ad = sub(sub(atac_h5ad_pattern, "\\{CHR\\}", chrom), "\\{CONDITION\\}", condition)

            call tasks.generate_peak_gene_pairs {
                input:
                    chromosome = chrom,
                    rna_h5ad = gex_h5ad,
                    atac_h5ad = atac_h5ad,
                    rna_nonzero_prop = sub(sub(gex_nonzero_prop_pattern, "\\{CELL_TYPE\\}", cell_type), "\\{CONDITION\\}", condition),
                    atac_nonzero_prop = sub(sub(sub(atac_nonzero_prop_pattern, "\\{CHR\\}", chrom), "\\{CELL_TYPE\\}", cell_type), "\\{CONDITION\\}", condition),
                    cell_type = cell_type,
                    all_cell_annot = all_cell_annot,
                    nonzero_prop_threshold = nonzero_prop_threshold,
                    min_chunk_size = min_chunk_size,
                    debug = debug,
                    docker = docker_suite,
                    zones = zones
            }

            scatter (chunk_pairs in generate_peak_gene_pairs.out_chunk_peak_gene_pairs) {
                call tasks.run_multilinker {
                    input:
                        method = method,
                        rna_h5ad = gex_h5ad,
                        atac_h5ad = atac_h5ad,
                        peak_gene_pairs = chunk_pairs,
                        cell_type = cell_type,
                        all_cell_annot = all_cell_annot,
                        obs_sample_id = obs_sample_id,
                        covar = covar,
                        covariates = covariates,
                        offset = offset,
                        binarize = binarize,
                        docker = docker_suite,
                        zones = zones
                }
            }

            call tasks.combine_results as combine_results_chunk {
                input:
                    results = run_multilinker.out_results,
                    prefix = "~{method}.~{cell_type}.~{condition}.~{chrom}",
                    docker = docker_suite,
                    zones = zones
            }
        }

        call tasks.combine_results {
            input:
                results = combine_results_chunk.out_results,
                prefix = "~{method}.~{cell_type}.~{condition}",
                bgzip = true,
                docker = docker_suite,
                zones = zones
        }
    }
}
