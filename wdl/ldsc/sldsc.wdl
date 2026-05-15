version 1.0

import "mesc.tasks.wdl" as tasks

workflow sldsc {
    input {
        File cell_type_list
        File phenotype_summary
        String ldscore_pattern
        File baselineLD
        File weights
        File frqfile
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
            File ldscore = sub(ldscore_pattern, "\\{CELL_TYPE\\}", cell_type)
            String ldscore_prefix = basename(ldscore, ".tar.gz")
            String prefix = "~{chunk_prefix}.~{ldscore_prefix}"

            call run_sldsc {
                input:
                    phenotype_list = pheno_list,
                    cell_type = cell_type,
                    sumstats = sumstats,
                    ldscore = ldscore,
                    baselineLD = baselineLD,
                    weights = weights,
                    frqfile = frqfile,
                    prefix = prefix,
                    ldscore_prefix = ldscore_prefix,
                    docker = docker_ldsc,
                    zones = zones
            }
        }

        call combine_results {
            input:
                prefix = chunk_prefix,
                all_results = run_sldsc.out_all_results,
                docker = docker_suite,
                zones = zones
        }
    }

    call combine_results as combine_results_all {
        input:
            prefix = out_prefix,
            all_results = combine_results.out_all_results,
            bgzip = true,
            docker = docker_suite,
            zones = zones
    }
}

task run_sldsc {
    input {
        File phenotype_list
        String cell_type
        Array[File] sumstats
        File ldscore
        File baselineLD
        File weights
        File frqfile
        String prefix
        String ldscore_prefix
        String docker
        String zones
    }

    command <<<
        set -e

        export PYENV_ROOT="$HOME/.pyenv"
        [[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"
        eval "$(pyenv init - bash)"

        tar xzvf ~{ldscore} --strip-components=1 && \
        tar xzvf ~{baselineLD} --strip-components=1 && \
        tar xzvf ~{weights} --strip-components=1 && \
        tar xzvf ~{frqfile}

        mkdir -p sumstats out
        ln -s ~{sep=' ' sumstats} sumstats/

        cat << "__EOF__" > script.awk
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            print "Phenotype", "Cell_Type", $0
        }
        NR > 1 {
            print phenotype, cell_type, $0
        }
        __EOF__

        while IFS=$'\t' read -r pheno sumstats_path
        do
            echo "Processing phenotype: $pheno"
            prefix=$pheno.~{ldscore_prefix}
            fname=$(basename $sumstats_path)

            python ~/ldsc/ldsc.py \
            --h2 sumstats/$fname \
            --w-ld-chr weights.hm3_noMHC. \
            --ref-ld-chr ~{ldscore_prefix}.chr,baselineLD. \
            --overlap-annot \
            --frqfile-chr 1000G.EUR.hg38. \
            --print-coefficients \
            --out $prefix

            awk -f script.awk -v phenotype=$pheno -v cell_type=~{cell_type} $prefix.results > out/$prefix.all.results
        done < ~{phenotype_list}

        awk 'NR == 1 || FNR > 1' out/*.all.results > ~{prefix}.all.results
    >>>

    output {
        File out_all_results = "~{prefix}.all.results"
        Array[File] out_log = glob("*.log")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "8 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task combine_results {
    input {
        String prefix
        Array[File] all_results
        Boolean bgzip = false
        String docker
        String zones
    }

    command <<<
        set -e

        awk 'NR == 1 || FNR > 1' ~{sep=' ' all_results} > ~{prefix}.all.results

        if [[ ~{bgzip} == "true" ]]
        then
            bgzip ~{prefix}.all.results
        fi
    >>>

    output {
        File out_all_results = if bgzip then "~{prefix}.all.results.gz" else "~{prefix}.all.results"
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