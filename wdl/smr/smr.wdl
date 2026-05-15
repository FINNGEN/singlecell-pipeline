version 1.0

workflow smr {
    input {
        File phenotype_summary
        File cell_type_list
        String bfile
        String beqtl_pattern
        String ma_sumstats_pattern
        String out_prefix
        Float h2_z_threshold = 2
        Int num_gw_significant_threshold = 10
        Int chunk_size = 10
        Boolean debug = false
        String docker_suite
        String zones
    }

    Array[String] cell_types = if debug then [read_lines(cell_type_list)[0]] else read_lines(cell_type_list)

    call preprocess_phenotypes {
        input:
            phenotype_summary = phenotype_summary,
            ma_sumstats_pattern = ma_sumstats_pattern,
            h2_z_threshold = h2_z_threshold,
            num_gw_significant_threshold = num_gw_significant_threshold,
            chunk_size = chunk_size,
            debug = debug,
            docker = docker_suite,
            zones = zones
    }

    scatter (cell_type in cell_types) {
        String beqtl = sub(beqtl_pattern, "\\{CELL_TYPE\\}", cell_type)
        scatter (zipped in zip(preprocess_phenotypes.out_phenos, preprocess_phenotypes.out_sumstats)) {
            File pheno_list = zipped.left
            Array[File] sumstats = read_lines(zipped.right)

            call run_smr {
                input:
                    cell_type = cell_type,
                    phenotype_list = pheno_list,
                    bfile = bfile,
                    beqtl = beqtl,
                    sumstats = sumstats,
                    docker = docker_suite,
                    zones = zones
            }
        }

        call combine_results {
            input:
                prefix = "~{out_prefix}.~{cell_type}",
                smr_results = run_smr.out_smr,
                docker = docker_suite,
                zones = zones
        }
    }
}

task preprocess_phenotypes {
    input {
        File phenotype_summary
        String ma_sumstats_pattern
        Float h2_z_threshold = 2
        Int num_gw_significant_threshold = 10
        Int chunk_size = 10
        Boolean debug = false
        String docker
        String zones
    }

    command <<<
        set -e

        python3 << "__EOF__"
        import math
        import pandas as pd
        import re

        df = pd.read_csv("~{phenotype_summary}", sep="\t")
        df = df.loc[(df.H2_Z > ~{h2_z_threshold}) & (df.num_gw_significant > ~{num_gw_significant_threshold}), :]
        df["sumstats"] = df["phenocode"].apply(lambda x: re.sub(r"\{PHENO\}", x, "~{ma_sumstats_pattern}"))

        df[["sumstats"]].to_csv("sumstats.tsv", sep="\t", index=False, header=False)

        # chunk params
        chunk_size= ~{chunk_size}
        total_lines = len(df.index)
        n_chunks = math.ceil(total_lines / chunk_size)
        lines_per_chunk = math.ceil(total_lines / n_chunks)

        if ~{true='True' false='False' debug}:
            n_chunks = 1

        for i in range(n_chunks):
            start_line = i * lines_per_chunk
            end_line = min(start_line + lines_per_chunk, total_lines)
            df[["phenocode", "sumstats"]].iloc[start_line:end_line].to_csv(f"pheno.{i}.tsv", sep="\t", index=False, header=False)
            df.sumstats.iloc[start_line:end_line].to_csv(f"sumstats.{i}.tsv", sep="\t", index=False, header=False)
        __EOF__
    >>>

    output {
        Array[File] out_phenos = glob("pheno.*.tsv")
        Array[File] out_sumstats = glob("sumstats.*.tsv")
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


task run_smr {
    input {
        String cell_type
        File phenotype_list
        String prefix = "~{basename(phenotype_list, '.tsv')}.~{cell_type}"
        String bfile
        File bed = "~{bfile}.bed"
        File bim = "~{bfile}.bim"
        File fam = "~{bfile}.fam"
        String beqtl
        File besd = "~{beqtl}.besd"
        File esi = "~{beqtl}.esi"
        File epi = "~{beqtl}.epi"
        Array[File] sumstats
        String docker
        String zones
    }

    command <<<
        set -e

        n_cpu=$(grep -c ^processor /proc/cpuinfo)

        mkdir -p bfile beqtl sumstats out
        ln -s ~{bed} ~{bim} ~{fam} bfile/
        ln -s ~{besd} ~{esi} ~{epi} beqtl/
        ln -s ~{sep=' ' sumstats} sumstats/

        sed -ie 's/_/:/g' ~{bim}

        while IFS=$'\t' read -r pheno sumstats_path
        do
            echo "Processing phenotype: $pheno"
            fname=$(basename $sumstats_path .gz)
            zcat sumstats/$fname.gz > sumstats/$fname

            smr \
            --bfile bfile/$(basename ~{bed} .bed) \
            --gwas-summary sumstats/$fname \
            --beqtl-summary beqtl/$(basename ~{beqtl} .besd) \
            --out $pheno \
            --thread-num $n_cpu
        done < ~{phenotype_list}

        awk -v cell_type="~{cell_type}" '
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            print "cell_type", "phenotype", $0
        }
        FNR > 1 {
            phenotype = FILENAME
            gsub(/\.smr$/, "", phenotype)
            print cell_type, phenotype, $0
        }
        ' *.smr > out/~{prefix}.smr
    >>>

    output {
        File out_smr = "out/~{prefix}.smr"
        Array[File] out_list = glob("*.snp_failed_freq_ck.list")
    }

    runtime {
        docker: docker
        cpu: 16
        memory: "24 GB"
        disks: "local-disk 100 HDD"
        zones: zones
        preemptible: 2
    }
}

task combine_results {
    input {
        String prefix
        Array[File] smr_results
        String docker
        String zones
    }

    command <<<
        set -e

        awk 'NR == 1 || FNR > 1' ~{sep=' ' smr_results} | bgzip -c > ~{prefix}.smr.gz
    >>>

    output {
        File out_smr = "~{prefix}.smr.gz"
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
