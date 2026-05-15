version 1.0

import "tensorqtl.tasks.wdl" as tasks

workflow tensorqtl {
    input {
        String bed_pattern
        String pfile_all
        String pfile_pattern
        Array[String] chromosomes
        File covar
        Array[String] covariates
        String? covar_sample_id
        String cell_type = "predicted.celltype.l1.CD4_T"
        File downsample_condition_list
        Array[String] downsample_conditions = read_lines(downsample_condition_list)
        Boolean acat = false
        Int min_K = 0
        Int max_K = 60
        Float maf_threshold = 0.01
        Int mac_threshold = 0
        Int window = 1000000
        String docker_peer
        String docker_qtl
        String docker_suite
        String docker_hail
        String zones
    }

    scatter (condition in downsample_conditions) {
        File psam = pfile_all + ".psam"
        File bed = sub(sub(bed_pattern, "\\{CELL_TYPE\\}", cell_type), "\\{CONDITION\\}", condition)
        String prefix = basename(bed, ".bed.gz")

        call tasks.count_n_samples {
            input:
                bed = bed,
                docker = docker_suite,
                zones = zones
        }

        if (count_n_samples.n_samples > 10) {
            Int max_K_capped = if (count_n_samples.n_samples < 2 * max_K) then count_n_samples.n_samples / 2 else max_K
            scatter (i in range(max_K_capped + 1 - min_K)) {
                Int K = i + min_K
                call tasks.compute_loco_expression_pcs as compute_expression_pcs {
                    input:
                        psam = psam,
                        bed = bed,
                        covar = covar,
                        covariates = covariates,
                        covar_sample_id = covar_sample_id,
                        K = K,
                        prefix = "~{prefix}.K~{K}",
                        zones = zones,
                        docker = docker_suite
                }

                if (compute_expression_pcs.n_samples > 10) {
                    if (!acat) {
                        scatter (chrom in chromosomes) {
                            call tasks.map as map_cis {
                                input:
                                    mode = "cis",
                                    pfile = sub(pfile_pattern, "\\{CHR\\}", chrom),
                                    bed = bed,
                                    covar = compute_expression_pcs.out_covar,
                                    chrom = chrom,
                                    maf_threshold = maf_threshold,
                                    mac_threshold = mac_threshold,
                                    window = window,
                                    docker = docker_qtl,
                                    zones = zones
                            }
                        }

                        call tasks.calculate_qvalues {
                            input:
                                cis_df = map_cis.out,
                                prefix = "~{prefix}.K~{K}",
                                docker = docker_qtl,
                                zones = zones
                        }
                    }
                    if (acat) {
                        scatter (chrom in chromosomes) {
                            call tasks.map as map_cis_nominal {
                                input:
                                    mode = "cis_nominal",
                                    pfile = sub(pfile_pattern, "\\{CHR\\}", chrom),
                                    bed = bed,
                                    covar = compute_expression_pcs.out_covar,
                                    chrom = chrom,
                                    maf_threshold = maf_threshold,
                                    mac_threshold = mac_threshold,
                                    window = window,
                                    supress_zfiles = true,
                                    docker = docker_qtl,
                                    zones = zones
                            }
                        }

                        call tasks.calculate_qvalues_acat {
                            input:
                                cis_df = map_cis_nominal.out_acat,
                                prefix = "~{prefix}.K~{K}",
                                docker = docker_suite,
                                zones = zones
                        }
                    }
                }
            }

            call tasks.count_egenes {
                input:
                    result_paths = flatten([select_all(calculate_qvalues.out), select_all(calculate_qvalues_acat.out)]),
                    prefix = prefix,
                    docker = docker_hail,
                    zones = zones
            }
        }
    }

    call tasks.optimize_n_peer {
        input:
            result_paths = select_all(count_egenes.out_txt),
            prefix = sub(basename(bed_pattern, ".inv.bed.gz"), "\\.\\{CELL_TYPE\\}", ""),
            docker = docker_suite,
            zones = zones
    }
}
