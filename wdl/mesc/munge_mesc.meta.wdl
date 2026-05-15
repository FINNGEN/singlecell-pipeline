version 1.0

import "mesc.tasks.wdl" as tasks

workflow munge_mesc {
    input {
        Array[String] chromosomes
        File cell_type_list
        Array[String] sumstats_pattern
        String out_prefix
        Boolean debug = false
        String docker_suite
        String docker_ldsc
        String zones
    }

    Array[String] cell_types = if debug then [read_lines(cell_type_list)[0]] else read_lines(cell_type_list)

    scatter (cell_type in cell_types) {
        String prefix = "~{out_prefix}.~{cell_type}.mesc"

        scatter (chrom in chromosomes) {
            String prefix_chrom = "~{prefix}.~{sub(chrom, '^chr', '')}"
            call concat_sumstats {
                input:
                    cell_type = cell_type,
                    chrom = chrom,
                    sumstats_pattern = sumstats_pattern,
                    prefix_chrom = prefix_chrom,
                    docker = docker_suite,
                    zones = zones
            }

            call tasks.compute_expscore {
                input:
                    sumstats = concat_sumstats.out_sumstats,
                    prefix = prefix,
                    prefix_chrom = prefix_chrom,
                    docker = docker_ldsc,
                    zones = zones
            }
        }

        call tasks.tar_expscore {
            input:
                prefix = prefix,
                expscore_files = flatten([
                    compute_expscore.out_G,
                    compute_expscore.out_ave,
                    compute_expscore.out_expscore,
                    compute_expscore.out_gannot,
                    compute_expscore.out_hsq
                ]),
                docker = docker_suite,
                zones = zones
        }
    }
}

task concat_sumstats {
    input {
        String cell_type
        String chrom
        Array[String] sumstats_pattern
        String prefix_chrom
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        python3 << "__EOF__"
        import pandas as pd
        import re

        sumstats_pattern = "~{sep=',' sumstats_pattern}".split(",")
        df = pd.DataFrame({"sumstats": sumstats_pattern})
        df["sumstats"] = df["sumstats"].apply(lambda x: re.sub(r"\{CELL_TYPE\}", "~{cell_type}", x))
        df["sumstats"] = df["sumstats"].apply(lambda x: re.sub(r"\{CHROM\}", "~{chrom}".replace("chr", ""), x))
        df.to_csv("sumstats.txt", index=False, header=False)
        __EOF__

        mkdir -p sumstats
        for sumstats in $(cat sumstats.txt)
        do
            gcloud storage cp $sumstats sumstats/ && \
            gunzip sumstats/$(basename $sumstats)
        done

        awk '
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            print
        }
        FNR > 1 {
            print | "sort -k2,2n -k1,1 -k5,5n"
        }
        ' sumstats/* | bgzip -c > ~{prefix_chrom}.gz && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_sumstats = "~{prefix_chrom}.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "8 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}
