version 1.0

workflow compute_ldscore {
    input {
        Array[String] chromosomes
        File cell_type_list
        String annot_bed_pattern
        String bfile_pattern
        File hm3_snplist
        String out_prefix
        Boolean debug = false
        String docker_suite
        String docker_ldsc
        String zones
    }

    Array[String] cell_types = if debug then [read_lines(cell_type_list)[0]] else read_lines(cell_type_list)

    scatter (cell_type in cell_types) {
        File annot_bed = sub(annot_bed_pattern, "\\{CELL_TYPE\\}", cell_type)

        scatter (chrom in chromosomes) {
            String bfile = sub(bfile_pattern, "\\{CHR\\}", sub(chrom, "^chr", ""))
            String prefix_chrom = "~{out_prefix}.~{cell_type}.~{chrom}"
            call make_annot {
                input:
                    annot_bed = annot_bed,
                    bfile = bfile,
                    prefix = prefix_chrom,
                    docker = docker_ldsc,
                    zones = zones
            }

            call compute_l2 {
                input:
                    bfile = bfile,
                    annot = make_annot.out_annot,
                    snp = hm3_snplist,
                    prefix = prefix_chrom,
                    docker = docker_ldsc,
                    zones = zones
            }
        }

        call tar_ldscore {
            input:
                prefix = "~{out_prefix}.~{cell_type}",
                ldscore_files = flatten([
                    make_annot.out_annot,
                    compute_l2.out_ldscore,
                    compute_l2.out_M,
                    compute_l2.out_M_5_50
                ]),
                docker = docker_suite,
                zones = zones
        }
    }
}

task make_annot {
    input {
        File annot_bed
        String prefix = basename(annot_bed, ".bed")
        String bfile
        File bim = "~{bfile}.bim"
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        export PYENV_ROOT="$HOME/.pyenv"
        [[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"
        eval "$(pyenv init - bash)"

        ~/ldsc/make_annot.py \
        --bed-file ~{annot_bed} \
        --bimfile ~{bim} \
        --annot-file ~{prefix}.annot.gz && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_annot = "~{prefix}.annot.gz"
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

task compute_l2 {
    input {
        String bfile
        File bed = "~{bfile}.bed"
        File bim = "~{bfile}.bim"
        File fam = "~{bfile}.fam"
        String bfile_prefix = basename(bfile)
        File annot
        File snp
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        export PYENV_ROOT="$HOME/.pyenv"
        [[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"
        eval "$(pyenv init - bash)"

        mkdir bfile && \
        ln -s ~{bed} ~{bim} ~{fam} bfile/ && \
        ~/ldsc/ldsc.py \
        --l2 \
        --bfile ./bfile/~{bfile_prefix} \
        --ld-wind-cm 1 \
        --annot ~{annot} \
        --thin-annot \
        --out ~{prefix} \
        --print-snps ~{snp} && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_ldscore = "~{prefix}.l2.ldscore.gz"
        File out_M = "~{prefix}.l2.M"
        File out_M_5_50 = "~{prefix}.l2.M_5_50"
    }

    runtime {
        docker: docker
        cpu: 4
        memory: "16 GB"
        disks: "local-disk 20 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task tar_ldscore {
    input {
        String prefix
        Array[File] ldscore_files
        String docker
        String zones
    }

    command <<<
        set -e

        mkdir -p ~{prefix}
        cp ~{sep=' ' ldscore_files} ~{prefix}/
        tar czvf ~{prefix}.tar.gz ~{prefix}
    >>>

    output {
        File out_tar = "~{prefix}.tar.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}