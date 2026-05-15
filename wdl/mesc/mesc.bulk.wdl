version 1.0

import "mesc.tasks.wdl" as tasks

workflow mesc {
    input {
        File phenotype_summary
        String expscore
        String out_prefix
        Int chunk_size = 25
        Boolean debug = false
        String docker_suite
        String docker_ldsc
        String zones
    }

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
        String expscore_prefix = basename(expscore, ".tar.gz")
        String prefix = "~{chunk_prefix}.~{expscore_prefix}"

        call tasks.run_mesc {
            input:
                phenotype_list = pheno_list,
                cell_type = "bulk",
                sumstats = sumstats,
                expscore = expscore,
                prefix = prefix,
                expscore_prefix = expscore_prefix,
                docker = docker_ldsc,
                zones = zones
        }
    }

    call tasks.combine_results as combine_results_all {
        input:
            prefix = out_prefix,
            all_h2med = run_mesc.out_all_h2med,
            categories_h2med = run_mesc.out_categories_h2med,
            bgzip = true,
            docker = docker_suite,
            zones = zones
    }
}
