version 1.0

import "multilinker.tasks.wdl" as tasks

workflow multilinker {
    input {
        String method = "open4gene"
        Array[String] chromosomes
        File cell_type_list
        String gex_h5ad
        String atac_h5ad_pattern
        String gex_nonzero_prop_pattern
        String atac_nonzero_prop_pattern
        String all_cell_annot
        String obs_sample_id
        String covar_pattern
        String atac_lsi_pattern
        File n_lsi
        Array[String] covariates
        String offset = "NULL"
        String test = "score"
        Boolean use_spa = false
        Float nonzero_prop_threshold = 0.05
        Int chunk_size = 25
        Boolean binarize = false
        Boolean debug = false
        String docker_suite
        String zones
    }

    Array[String] cell_types = if debug then [read_lines(cell_type_list)[0]] else read_lines(cell_type_list)

    scatter (cell_type in cell_types) {
        scatter (chrom in chromosomes) {
            String atac_h5ad = sub(atac_h5ad_pattern, "\\{CHR\\}", chrom)

            call tasks.generate_peak_gene_pairs {
                input:
                    chromosome = chrom,
                    rna_nonzero_prop = sub(gex_nonzero_prop_pattern, "\\{CELL_TYPE\\}", cell_type),
                    atac_nonzero_prop = sub(atac_nonzero_prop_pattern, "\\{CELL_TYPE\\}", cell_type),
                    nonzero_prop_threshold = nonzero_prop_threshold,
                    chunk_size = chunk_size,
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
                        covar = sub(covar_pattern, "\\{CELL_TYPE\\}", cell_type),
                        atac_lsi = sub(atac_lsi_pattern, "\\{CELL_TYPE\\}", cell_type),
                        n_lsi = n_lsi,
                        covariates = covariates,
                        offset = offset,
                        test = test,
                        use_spa = use_spa,
                        binarize = binarize,
                        docker = docker_suite,
                        zones = zones
                }
            }

            call tasks.combine_results as combine_results_chunk {
                input:
                    results = run_multilinker.out_results,
                    prefix = "~{method}.~{cell_type}.~{chrom}",
                    docker = docker_suite,
                    zones = zones
            }
        }

        call tasks.combine_results {
            input:
                results = combine_results_chunk.out_results,
                prefix = "~{method}.~{cell_type}",
                bgzip = true,
                docker = docker_suite,
                zones = zones
        }
    }
}
