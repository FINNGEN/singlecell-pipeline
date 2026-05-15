version 1.0

import "tensorqtl.tasks.wdl" as tasks

workflow tensorqtl {
    input {
        String bed_pattern
        String pfile_all
        String pfile_pattern
        Array[String] chromosomes
        String cell_type = "predicted.celltype.l1.CD4_T"
        File downsample_condition_list
        Array[String] downsample_conditions = read_lines(downsample_condition_list)
        File covar
        Array[String] covariates
        String? covar_sample_id
        File? n_peer
        Float maf_threshold = 0.01
        Int mac_threshold = 0
        Int window = 1000000
        Array[String] models = ["additive"]
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

        call tasks.peer_factor {
            input:
                psam = psam,
                bed = bed,
                covar = covar,
                covariates = covariates,
                covar_sample_id = covar_sample_id,
                n_peer = n_peer,
                cell_type = cell_type,
                prefix = prefix,
                docker = docker_peer,
                zones = zones
        }

        if (peer_factor.n_samples > 10) {
            scatter (model in models) {
                scatter (chrom in chromosomes) {
                    String pfile = sub(pfile_pattern, "\\{CHR\\}", chrom)
                    call tasks.map as map_cis_nominal {
                        input:
                            mode = "cis_nominal",
                            pfile = pfile,
                            bed = bed,
                            covar = peer_factor.out_covar,
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
                            covar = peer_factor.out_covar,
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
                        prefix = if model == "additive" then "~{prefix}.~{condition}.cis_qtl_pairs" else "~{prefix}.~{condition}.~{model}.cis_qtl_pairs",
                        docker = docker_suite,
                        zones = zones
                }

                call tasks.calculate_qvalues_acat {
                    input:
                        cis_df = map_cis_nominal.out_acat,
                        prefix = if model == "additive" then "~{prefix}.~{condition}" else "~{prefix}.~{condition}.~{model}",
                        docker = docker_suite,
                        zones = zones
                }
            }
        }
    }
}
