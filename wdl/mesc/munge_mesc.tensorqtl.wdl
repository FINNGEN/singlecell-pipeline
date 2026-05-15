version 1.0

import "mesc.tasks.wdl" as tasks

workflow munge_mesc {
    input {
        Array[String] chromosomes
        File cell_type_list
        String sumstats_pattern
        String bed_pattern
        File tss_full_bed
        File rsid_annot
        Boolean debug = false
        String docker_suite
        String docker_ldsc
        String zones
    }

    Array[String] cell_types = if debug then [read_lines(cell_type_list)[0]] else read_lines(cell_type_list)

    scatter (cell_type in cell_types) {
        File sumstats = sub(sumstats_pattern, "\\{CELL_TYPE\\}", cell_type)
        File bed = sub(bed_pattern, "\\{CELL_TYPE\\}", cell_type)
        String prefix = "~{basename(sumstats, '.cis_qtl_pairs.tsv.gz')}.mesc"

        scatter (chrom in chromosomes) {
            String prefix_chrom = "~{prefix}.~{sub(chrom, '^chr', '')}"
            call munge_tensorqtl {
                input:
                    chrom = chrom,
                    sumstats = sumstats,
                    bed = bed,
                    tss_full_bed = tss_full_bed,
                    rsid_annot = rsid_annot,
                    prefix_chrom = prefix_chrom,
                    docker = docker_suite,
                    zones = zones
            }

            call tasks.compute_expscore {
                input:
                    sumstats = munge_tensorqtl.out_sumstats,
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

task munge_tensorqtl {
    input {
        String chrom
        File sumstats
        File sumstats_tbi = "~{sumstats}.tbi"
        File bed
        File tss_full_bed
        File rsid_annot
        String prefix_chrom
        String docker
        String zones
    }

    command <<<
        set -e

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        n_samples=$(zcat ~{bed} | head -n1 | tr '\t' '\n' | tail -n+5 | wc -l)

        tabix -h ~{sumstats} ~{chrom} | awk -v n=$n_samples '
        BEGIN {
            OFS = "\t"
            filenum = 1
            print "GENE", "GENE_COORD", "SNP", "CHR", "SNP_COORD", "N", "Z"
        }
        FNR == NR {
            gene_coord[$4] = int(($2 + $3) / 2)
            next
        }
        FNR == 1 && FNR < NR {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
            filenum++
        }
        FNR > 1 && filenum == 2 {
            rsid[$col["variant_id"]] = $col["rsid"]
        }
        FNR > 1 && filenum == 3 && $col["slope_se"] != "NA" {
            if ($col["variant_id"] in rsid) {
                snp = rsid[$col["variant_id"]]
            } else {
                snp = $col["variant_id"]
            }
            chrom = $col["#chrom"]
            sub(/^chr/, "", chrom)

            gene = $col["phenotype_id"]
            if (gene in gene_coord) {
                gene_mid_pos = gene_coord[gene]
            } else if (gene ~ /^chr/) {
                peak_coord = split(gene, a, "-")
                gene_mid_pos = int((a[2] + a[3]) / 2)
            } else {
                print "Error: Gene not found in TSS bed:", gene > "/dev/stderr"
                exit 1
            }
            z = $col["slope"] / $col["slope_se"]
            print gene, gene_mid_pos, snp, chrom, $col["position"], n, z | "sort -k2,2n -k5,5n"
        }
        ' <(cat ~{tss_full_bed}) <(zcat ~{rsid_annot}) - | bgzip -c > ~{prefix_chrom}.gz && \
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
