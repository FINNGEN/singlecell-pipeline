version 1.0

workflow subset_multiplex_vcf {
    input {
        Array[String] chromosomes
        File lexsorted_vcf
        File lexsorted_vcf_tbi = lexsorted_vcf + ".tbi"
        File mapping
        File input_sample_sheet
        String docker
        String zones
    }

    call config_samples {
        input:
            input_sample_sheet = input_sample_sheet,
            docker = docker,
            zones = zones
    }

    scatter (samples in config_samples.multiplexed_samples) {
        call subset_vcf {
            input:
                chromosomes = chromosomes,
                vcf = lexsorted_vcf,
                vcf_tbi = lexsorted_vcf_tbi,
                mapping = mapping,
                samples = samples,
                docker = docker,
                zones = zones
        }
    }

    output {
        Array[File] out_subset_vcfs = subset_vcf.out_vcf
    }
}

task config_samples {
    input {
        File input_sample_sheet
        String docker
        String zones
    }

    command <<<
        set -e

        # prepare multiplexed sample config
        cut -f1 -d, ~{input_sample_sheet} | cut -f1 -d "-" | tail -n+2 | sort | uniq > multiplexed_samples.txt

    >>>

    output {
        Array[String] multiplexed_samples = read_lines("multiplexed_samples.txt")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 1 HDD"
        zones: zones
        preemptible: 2
    }
}

task subset_vcf {
    input {
        Array[String] chromosomes
        File vcf
        File vcf_tbi
        File mapping
        String samples

        String docker
        String zones
    }

    command <<<
        set -e

        # remap back to FFIDs for demultiplexing
        awk 'NR > 1{print $1, $2}' ~{mapping} > update_ids.txt

        bcftools reheader -s update_ids.txt ~{vcf} |
        bcftools view --types 'snps' --min-ac 1:minor -s $(echo ~{samples} | tr '_' ',') | awk '
        BEGIN {
            OFS = "\t"
        }
        $0 ~ /^#/ {
            print
            next
        }
        {
            for (i = 10; i <= NF; i++) {
                if ($i !~ /:\.$/ && $9 != "GT") {
                    continue
                }
                gt = substr($i, 1, 1) + substr($i, 1, 3)
                if (gt == 0) {
                    $i = substr($i, 1, 3)":1,0,0"
                } else if (gt == 1) {
                    $i = substr($i, 1, 3)":0,1,0"
                } else if (gt == 2) {
                    $i = substr($i, 1, 3)":0,0,1"
                }
            }
            if ($9 == "GT") {
                $9 = "GT:GP"
            }
            print
        }' | bgzip -c > ~{samples}.vcf.gz
        tabix ~{samples}.vcf.gz

    >>>

    output {
        File out_vcf = samples + ".vcf.gz"
        File out_vcf_tbi = samples + ".vcf.gz.tbi"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "7 GB"
        disks: "local-disk 100 HDD"
        zones: zones
        preemptible: 2
    }
}
