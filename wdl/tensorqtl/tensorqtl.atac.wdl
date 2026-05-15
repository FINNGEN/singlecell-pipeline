version 1.0

import "tensorqtl.tasks.wdl" as tasks

workflow tensorqtl {
    input {
        String bed_pattern
        String pfile_all
        String pfile_pattern
        Array[String] chromosomes
        File cell_type_list
        Array[String] cell_types = read_lines(cell_type_list)
        File covar
        Array[String] covariates
        String? covar_sample_id
        File? n_peer
        Float maf_threshold = 0.01
        Int mac_threshold = 0
        Int window = 1000000
        Array[String] models = ["additive", "dominant", "recessive"]
        String out_ht_dir
        Boolean overwrite
        Boolean write_ht = true
        String docker_peer
        String docker_qtl
        String docker_suite
        String docker_hail
        String zones
    }

    scatter (cell_type in cell_types) {
        File psam = pfile_all + ".psam"
        File bed = sub(bed_pattern, "\\{CELL_TYPE\\}", cell_type)
        String prefix = basename(bed, ".bed.gz")

        call tasks.compute_loco_expression_pcs as compute_expression_pcs {
            input:
                psam = psam,
                bed = bed,
                covar = covar,
                covariates = covariates,
                covar_sample_id = covar_sample_id,
                n_peer = n_peer,
                cell_type = cell_type,
                prefix = prefix,
                docker = docker_suite,
                zones = zones
        }

        if (compute_expression_pcs.n_samples > 50) {
            scatter (model in models) {
                scatter (chrom in chromosomes) {
                    String pfile = sub(pfile_pattern, "\\{CHR\\}", chrom)
                    call tasks.map as map_cis_nominal {
                        input:
                            mode = "cis_nominal",
                            pfile = pfile,
                            bed = bed,
                            covar = compute_expression_pcs.out_covar,
                            chrom = chrom,
                            maf_threshold = maf_threshold,
                            mac_threshold = mac_threshold,
                            window = window,
                            dominant = model == "dominant",
                            recessive = model == "recessive",
                            docker = docker_qtl,
                            zones = zones
                    }

                    call tasks.map as map_cis {
                        input:
                            mode = "cis",
                            pfile = pfile,
                            bed = bed,
                            covar = compute_expression_pcs.out_covar,
                            chrom = chrom,
                            maf_threshold = maf_threshold,
                            mac_threshold = mac_threshold,
                            window = window,
                            dominant = model == "dominant",
                            recessive = model == "recessive",
                            docker = docker_qtl,
                            zones = zones
                    }

                    call tasks.map as map_trans {
                        input:
                            mode = "trans",
                            pfile = pfile_all,
                            bed = bed,
                            covar = compute_expression_pcs.out_covar,
                            chrom = chrom,
                            maf_threshold = maf_threshold,
                            mac_threshold = mac_threshold,
                            window = window,
                            dominant = model == "dominant",
                            recessive = model == "recessive",
                            docker = docker_qtl,
                            zones = zones
                    }

                }

                call tasks.combine_tsv as combine_tsv_cis {
                    input:
                        mode = "cis_nominal",
                        tsv_paths = map_cis_nominal.out,
                        cell_type = cell_type,
                        prefix = if model == "additive" then "~{prefix}.cis_qtl_pairs" else "~{prefix}.~{model}.cis_qtl_pairs",
                        docker = docker_suite,
                        zones = zones
                }

                call tasks.combine_tsv as combine_tsv_trans {
                    input:
                        mode = "trans",
                        tsv_paths = map_trans.out,
                        cell_type = cell_type,
                        prefix = if model == "additive" then "~{prefix}.trans_qtl_pairs" else "~{prefix}.~{model}.trans_qtl_pairs",
                        docker = docker_suite,
                        zones = zones
                }

                call tasks.calculate_qvalues {
                    input:
                        cis_df = map_cis.out,
                        prefix = if model == "additive" then prefix else "~{prefix}.~{model}",
                        docker = docker_qtl,
                        zones = zones
                }

                call tasks.calculate_qvalues_acat {
                    input:
                        cis_df = map_cis_nominal.out_acat,
                        prefix = if model == "additive" then prefix else "~{prefix}.~{model}",
                        docker = docker_suite,
                        zones = zones
                }
            }
        }
    }

    if (write_ht) {
        scatter (model in models) {
            call tasks.import_hail as import_hail_cis {
                input:
                    tsv_paths = flatten(select_all(combine_tsv_cis.out_tsv)),
                    out_ht = if model == "additive" then "~{out_ht_dir}/cis_nominal.ht" else "~{out_ht_dir}/cis_nominal.~{model}.ht",
                    model = model,
                    overwrite = overwrite,
                    docker = docker_hail,
                    zones = zones
            }

            call tasks.import_hail as import_hail_trans {
                input:
                    tsv_paths = flatten(select_all(combine_tsv_trans.out_tsv)),
                    out_ht = if model == "additive" then "~{out_ht_dir}/trans_nominal.ht" else "~{out_ht_dir}/trans_nominal.~{model}.ht",
                    model = model,
                    overwrite = overwrite,
                    docker = docker_hail,
                    zones = zones
            }
        }
    }
}
