version development

import "saige_qtl.tasks.wdl" as tasks

workflow saige_qtl_step1 {
    input {
        String out_dir
        File gex_features
        File obs
        String expression_bed_pattern
        String nonzero_prop_pattern
        String h5ad_pattern
        String covar_pattern
        String barcode_pattern
        File cell_type_list
        String grm_bfile
        String obs_sample_id
        String sample_id
        Array[String] covariates
        Array[String] sample_covariates
        Array[String] categorical_covariates = []
        Array[String] dynamic_covariates = []
        String offset = ""
        Boolean add_depth_covariate = false
        Int num_genes_per_chunk
        Boolean overwrite
        Boolean filter_genes = true
        Boolean debug = false
        String docker_saige
        String docker_suite
        String zones
    }

    Array[String] cell_types = if debug then [read_lines(cell_type_list)[0]] else read_lines(cell_type_list)

    scatter (cell_type in cell_types) {
        File expression_bed = sub(expression_bed_pattern, "\\{CELL_TYPE\\}", cell_type)
        File nonzero_prop = sub(nonzero_prop_pattern, "\\{CELL_TYPE\\}", cell_type)
        File covar = sub(covar_pattern, "\\{CELL_TYPE\\}", cell_type)
        File barcode = sub(barcode_pattern, "\\{CELL_TYPE\\}", cell_type)
        String prefix = basename(expression_bed, ".bed.gz")
        String chunk_list = "~{out_dir}/chunk/~{prefix}.chunks.txt"

        call tasks.preprocess_pseudobulk {
            input:
                gex_features = gex_features,
                expression_bed = expression_bed,
                nonzero_prop = nonzero_prop,
                covar = covar,
                sample_id = sample_id,
                num_genes_per_chunk = num_genes_per_chunk,
                filter_genes = filter_genes,
                prefix = prefix,
                chunk_list = chunk_list,
                docker = docker_suite,
                zones = zones
        }

        Array[Array[String]] gene_id_chunks = if debug then [read_tsv(preprocess_pseudobulk.out_chunks)[1]] else read_tsv(preprocess_pseudobulk.out_chunks)

        scatter (zipped in zip(range(length(gene_id_chunks)), gene_id_chunks)) {
            Int i = zipped.left
            Array[String] chunk = zipped.right
            String prefix_chunk = "~{prefix}.chunk~{i}"

            call tasks.preprocess_sc {
                input:
                    h5ad_pattern = h5ad_pattern,
                    chunk = chunk,
                    covar = covar,
                    covariates = if offset != "" then flatten([covariates, [offset]]) else covariates,
                    obs = obs,
                    barcode = barcode,
                    obs_sample_id = obs_sample_id,
                    sample_id = sample_id,
                    prefix = prefix_chunk,
                    docker = docker_suite,
                    zones = zones
            }

            if (!overwrite) {
                call tasks.check_step1_null {
                    input:
                        gene_ids = read_lines(preprocess_sc.out_gene_ids),
                        out_dir = sub(out_dir, "[/\\s]+$", ""),
                        prefix = prefix,
                        docker = docker_suite,
                        zones = zones
                }
            }

            Array[String] gene_ids = read_lines(select_first([check_step1_null.out_gene_ids, preprocess_sc.out_gene_ids]))

            scatter (gene_id in gene_ids) {
                Array[String] step1_covariates = flatten([covariates, preprocess_sc.peer_covariates, dynamic_covariates])

                call tasks.step1_null {
                    input:
                        trait_type = "count",
                        bfile = grm_bfile,
                        pheno = preprocess_sc.out,
                        sample_id = sample_id,
                        gene_id = gene_id,
                        covariates = step1_covariates,
                        sample_covariates = sample_covariates,
                        categorical_covariates = categorical_covariates,
                        dynamic_covariates = dynamic_covariates,
                        offset = offset,
                        n_cells = preprocess_sc.n_cells,
                        out_dir = sub(out_dir, "[/\\s]+$", ""),
                        prefix = "~{prefix}.~{gene_id}",
                        docker = docker_saige,
                        zones = zones
                }
            }
        }
    }
}
