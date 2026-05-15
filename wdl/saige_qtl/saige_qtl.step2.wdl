version development

import "saige_qtl.tasks.wdl" as tasks

workflow saige_qtl_step2 {
    input {
        String out_dir
        String expression_bed_pattern
        String chunk_list_pattern
        String sc_step1_rda_pattern
        String sc_step1_vr_pattern
        File cell_type_list
        File bgen_list
        Boolean skip_non_existing
        Float min_maf = 0.0
        Int min_mac = 20
        Boolean trans = false
        Array[String] dynamic_covariates = []
        Int step2_chunk_size = 0
        String tmp_bucket
        Boolean overwrite
        Boolean debug = false
        String docker_hail
        String docker_saige
        String docker_suite
        String zones
    }

    Boolean is_gxe = length(dynamic_covariates) > 0
    Array[File] bgenfiles = if debug then [read_lines(bgen_list)[0]] else read_lines(bgen_list)
    Array[String] cell_types = if debug then [read_lines(cell_type_list)[0]] else read_lines(cell_type_list)

    scatter (cell_type in cell_types) {
        File expression_bed = sub(expression_bed_pattern, "\\{CELL_TYPE\\}", cell_type)
        File chunk_list_raw = sub(chunk_list_pattern, "\\{CELL_TYPE\\}", cell_type)
        String prefix = basename(expression_bed, ".bed.gz")
        String sc_step1_rda_pattern_cell = sub(sc_step1_rda_pattern, "\\{CELL_TYPE\\}", cell_type)
        String sc_step1_vr_pattern_cell = sub(sc_step1_vr_pattern, "\\{CELL_TYPE\\}", cell_type)

        if (step2_chunk_size > 0) {
            call tasks.rechunk_step2_list {
                input:
                    chunk_list = chunk_list_raw,
                    step2_chunk_size = step2_chunk_size,
                    prefix = prefix,
                    docker = docker_suite,
                    zones = zones
            }
        }

        File chunk_list = select_first([rechunk_step2_list.out_chunk_list, chunk_list_raw])
        Array[Array[String]] gene_id_chunks = if debug then [read_tsv(chunk_list)[0]] else read_tsv(chunk_list)

        scatter (zipped in zip(range(length(gene_id_chunks)), gene_id_chunks)) {
            Int i = zipped.left
            Array[String] chunk = zipped.right
            String prefix_chunk = "~{prefix}.chunk~{i}"

            call tasks.generate_multi_gene_file {
                input:
                    step1_rda_pattern = sc_step1_rda_pattern_cell,
                    step1_vr_pattern = sc_step1_vr_pattern_cell,
                    chunk = chunk,
                    bgenfiles = bgenfiles,
                    expression_bed = expression_bed,
                    trans = trans,
                    prefix = prefix_chunk,
                    tmp_bucket = tmp_bucket,
                    overwrite = overwrite,
                    skip_non_existing = skip_non_existing,
                    docker = docker_suite,
                    zones = zones
            }

            Array[String] gene_ids = read_lines(generate_multi_gene_file.out_gene_ids)

            if (generate_multi_gene_file.out_ht == "") {
                scatter (bgen in generate_multi_gene_file.out_bgenfiles) {
                    call tasks.step2_test {
                        input:
                            gene_ids = gene_ids,
                            modelfiles = read_lines(generate_multi_gene_file.out_modelfiles),
                            varianceratiofiles = read_lines(generate_multi_gene_file.out_varianceratiofiles),
                            bgen = bgen,
                            cis_region = generate_multi_gene_file.out_cis_region,
                            min_maf = min_maf,
                            min_mac = min_mac,
                            is_gxe = is_gxe,
                            prefix = prefix_chunk,
                            docker = docker_saige,
                            zones = zones
                    }
                }

                if (defined(select_first(step2_test.out_results))) {
                    call tasks.write_step2_ht {
                        input:
                            gene_ids = gene_ids,
                            result_paths = flatten(select_all(step2_test.out_results)),
                            expression_bed = expression_bed,
                            cell_type = cell_type,
                            dynamic_covariates = dynamic_covariates,
                            prefix = prefix_chunk,
                            tmp_bucket = tmp_bucket,
                            docker = docker_hail,
                            zones = zones
                    }
                }
            }
        }

        call tasks.combine_cis {
            input:
                out_dir = sub(out_dir, "[/\\s]+$", ""),
                step2_ht_paths = flatten([generate_multi_gene_file.out_ht, select_all(write_step2_ht.out_ht)]),
                prefix = prefix,
                docker = docker_hail,
                zones = zones
        }

        if (!is_gxe) {
            call tasks.step3_gene_pvalue {
                input:
                    cis = combine_cis.out_cis,
                    prefix = prefix,
                    docker = docker_saige,
                    zones = zones
            }
        }

        if (is_gxe) {
            call tasks.step3_gxe_pvalue {
                input:
                    cis = combine_cis.out_cis,
                    prefix = prefix,
                    docker = docker_saige,
                    zones = zones
            }
        }

        File acat_p = select_first([step3_gene_pvalue.out_acat_p, step3_gxe_pvalue.out_acat_p])

        call tasks.postprocess {
            input:
                out_dir = sub(out_dir, "[/\\s]+$", ""),
                acat_p = acat_p,
                step2_ht_list = combine_cis.out_step2_ht_list,
                cis = combine_cis.out_cis,
                skip_gene_ids = generate_multi_gene_file.out_skip_gene_ids,
                prefix = prefix,
                is_gxe = is_gxe,
                docker = docker_suite,
                zones = zones
        }
    }
}
