version 1.0

import "mesc.tasks.wdl" as tasks

workflow munge_mesc {
    input {
        File probe_info
        String sumstats_pattern
        File rsid_annot
        String out_prefix
        Boolean debug = false
        String docker_suite
        String docker_ldsc
        String zones
    }

    call preprocess_probes {
        input:
            probe_info = probe_info,
            sumstats_pattern = sumstats_pattern,
            window = 1000000,
            docker = docker_suite,
            zones = zones
    }

    Array[File] phenotype_list = if debug then [preprocess_probes.out_phenos[0]] else preprocess_probes.out_phenos

    scatter (pheno_list in phenotype_list) {

        String prefix_chrom = '~{out_prefix}.~{sub(basename(pheno_list, ".tsv"), "^pheno\\.", "")}'

        call munge_pqtl {
            input:
                phenotype_list = pheno_list,
                rsid_annot = rsid_annot,
                prefix_chrom = prefix_chrom,
                docker = docker_suite,
                zones = zones
        }

        call tasks.compute_expscore {
            input:
                sumstats = munge_pqtl.out_sumstats,
                prefix = out_prefix,
                prefix_chrom = prefix_chrom,
                docker = docker_ldsc,
                zones = zones
        }
    }

    call tasks.tar_expscore {
        input:
            prefix = out_prefix,
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

task preprocess_probes {
    input {
        File probe_info
        String sumstats_pattern
        Int window = 1000000
        String docker
        String zones
    }

    command <<<
        set -e

        python3 << "__EOF__"
        import math
        import pandas as pd
        import re

        window = ~{window}
        autosomes = [str(i) for i in range(1, 23)]

        df = pd.read_csv("~{probe_info}", sep="\t")
        # make sure to use chromosome names without "chr" prefix
        df["chrom"] = df.chrom.astype(str).str.replace("chr", "")
        df = df.loc[df.chrom.isin(autosomes), :]        
        df["cis_start"] = (df.start - window).astype(int)
        df["cis_end"] = (df.end + window).astype(int)
        df["pheno_coord"] = ((df.start + df.end) // 2).astype(int)
        df["pheno"] = df.apply(lambda x: f"{x.Probe}.{x.phenotype_id}", axis=1)
        df["region"] = df.apply(lambda x: f"{x.chrom}:{x.cis_start}-{x.cis_end}", axis=1)
        df["sumstats"] = df["Probe"].apply(lambda x: re.sub(r"\{PHENO\}", x, "~{sumstats_pattern}"))

        # chunk by chromosomes
        for chrom in autosomes:
            df_chrom = df.loc[df.chrom == chrom, :][["pheno", "pheno_coord", "region", "sumstats"]]
            df_chrom.to_csv(f"pheno.{chrom}.tsv", sep="\t", index=False, header=False)
        __EOF__
    >>>

    output {
        Array[File] out_phenos = glob("pheno.*.tsv")
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

task munge_pqtl {
    input {
        File phenotype_list
        File rsid_annot
        String prefix_chrom
        String chrom_col = "CHR"
        String position_col = "POS"
        String id_col = "ID"
        String beta_col = "BETA"
        String se_col = "SE"
        String n_col = "N"
        String docker
        String zones
    }

    command <<<
        set -e

        mkdir -p sumstats

        while IFS=$'\t' read -r pheno pheno_coord region sumstats_path
        do
            export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
            echo "Processing phenotype: $pheno"

            # if header.txt doesn't exist, create it
            if [[ ! -f header.txt ]]
            then
                gcloud storage cat $sumstats_path | zcat | head -n1 | awk '
                BEGIN {
                    OFS = "\t"
                }
                {
                    print "phenotype_id", "pheno_coord", $0
                }' > header.txt
            fi

            tabix $sumstats_path $region | awk -v pheno=$pheno -v coord=$pheno_coord '
            BEGIN {
                OFS = "\t"
            }
            {
                print pheno, coord, $0
            }' > sumstats/$pheno.txt
        done < ~{phenotype_list}

        cat header.txt sumstats/*.txt | awk '
        BEGIN {
            OFS = "\t"
            print "GENE", "GENE_COORD", "SNP", "CHR", "SNP_COORD", "N", "Z"
        }
        FNR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
        }
        FNR > 1 && FNR == NR {
            rsid[$col["variant_id"]] = $col["rsid"]
        }
        FNR > 1 && FNR < NR && $col["~{se_col}"] != "NA" {
            gsub(/:/, "_", $col["~{id_col}"])
            if ($col["~{id_col}"] !~ /^chr/) {
                $col["~{id_col}"] = "chr"$col["~{id_col}"]
            }
            if ($col["~{id_col}"] in rsid) {
                snp = rsid[$col["~{id_col}"]]
            } else {
                snp = $col["~{id_col}"]
            }
            chrom = $col["~{chrom_col}"]
            sub(/^chr/, "", chrom)

            gene = $col["phenotype_id"]
            gene_mid_pos = $col["pheno_coord"]
            z = $col["~{beta_col}"] / $col["~{se_col}"]
            print gene, gene_mid_pos, snp, chrom, $col["~{position_col}"], $col["~{n_col}"], z | "sort -k2,2n -k1,1 -k5,5n"
        }
        ' <(zcat ~{rsid_annot}) - | bgzip -c > ~{prefix_chrom}.gz && \
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
