version 1.0

workflow munge_fg_sumstats {
    input {
        File phenotype_summary
        Float h2_z_threshold = 2
        Int num_gw_significant_threshold = 10
        Boolean debug = false
        String docker_suite
        String zones
    }

    call preprocess_phenotypes {
        input:
            phenotype_summary = phenotype_summary,
            h2_z_threshold = h2_z_threshold,
            num_gw_significant_threshold = num_gw_significant_threshold,
            docker = docker_suite,
            zones = zones
    }

    Array[String] sumstats_list = if debug then [preprocess_phenotypes.out_sumstats[0]] else preprocess_phenotypes.out_sumstats

    scatter (sumstats in sumstats_list) {
        call munge {
            input:
                sumstats = sumstats,
                docker = docker_suite,
                zones = zones
        }
    }
}

task preprocess_phenotypes {
    input {
        File phenotype_summary
        Float h2_z_threshold = 2
        Int num_gw_significant_threshold = 10
        String docker
        String zones
    }

    command <<<
        set -e

        python3 << "__EOF__"
        import pandas as pd

        df = pd.read_csv("~{phenotype_summary}", sep="\t")
        df = df.loc[(df.H2_Z > ~{h2_z_threshold}) & (df.num_gw_significant > ~{num_gw_significant_threshold}), :]

        df[["sumstats"]].to_csv("sumstats.tsv", sep="\t", index=False, header=False)
        __EOF__
    >>>

    output {
        Array[String] out_sumstats = read_lines("sumstats.tsv")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}

task munge {
    input {
        File sumstats
        String prefix = basename(sumstats, ".gz")
        String docker
        String zones
    }

    command <<<
        set -e

        zcat ~{sumstats} | awk '
        BEGIN {
            FS = "\t"
            OFS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
            print "SNP", "A1", "A2", "freq", "b", "se", "p", "n"
        }
        NR > 1 {
            if ($col["#chrom"] == 23) {
                chrom = "X"
            } else {
                chrom = $col["#chrom"]
            }
            variant = "chr"chrom":"$col["pos"]":"$col["ref"]":"$col["alt"]
            print variant, $col["alt"], $col["ref"], $col["af_alt"], $col["beta"], $col["sebeta"], $col["pval"], "NA"
        }
        ' | bgzip -c > ~{prefix}.ma.gz
    >>>

    output {
        File out_ma = "~{prefix}.ma.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}