version 1.0

import "mesc.tasks.wdl" as tasks

workflow mesc {
    input {
        File cell_type_list
        File phenotype_summary
        String expscore_pattern
        String out_prefix
        Int chunk_size = 25
        Boolean debug = false
        String docker_suite
        String docker_ldsc
        String zones
    }

    Array[String] cell_types = if debug then [read_lines(cell_type_list)[0]] else read_lines(cell_type_list)

    call tasks.preprocess_phenotypes {
        input:
            phenotype_summary = phenotype_summary,
            chunk_size = chunk_size,
            docker = docker_suite,
            zones = zones
    }

    scatter (zipped in zip(preprocess_phenotypes.out_phenos, preprocess_phenotypes.out_sumstats)) {
        File pheno_list = zipped.left
        Array[File] sumstats = read_lines(zipped.right)
        String chunk_prefix = basename(pheno_list, ".txt")

        scatter (cell_type in cell_types) {
            File expscore = sub(expscore_pattern, "\\{CELL_TYPE\\}", cell_type)
            String expscore_prefix = basename(expscore, ".tar.gz")
            String prefix = "~{chunk_prefix}.~{expscore_prefix}"

            call tasks.run_mesc {
                input:
                    phenotype_list = pheno_list,
                    cell_type = cell_type,
                    sumstats = sumstats,
                    expscore = expscore,
                    prefix = prefix,
                    expscore_prefix = expscore_prefix,
                    docker = docker_ldsc,
                    zones = zones
            }
        }

        call tasks.combine_results {
            input:
                prefix = chunk_prefix,
                all_h2med = run_mesc.out_all_h2med,
                categories_h2med = run_mesc.out_categories_h2med,
                docker = docker_suite,
                zones = zones
        }
    }

    call tasks.combine_results as combine_results_all {
        input:
            prefix = out_prefix,
            all_h2med = combine_results.out_all_h2med,
            categories_h2med = combine_results.out_categories_h2med,
            bgzip = true,
            docker = docker_suite,
            zones = zones
    }
}
