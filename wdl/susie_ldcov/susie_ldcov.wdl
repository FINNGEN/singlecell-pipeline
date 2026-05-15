version 1.0

workflow susie_ldcov {
    input {
        String output_directory
        String bcor_output_directory_pattern
        Array[String] chromosomes
        File cell_type_list
        Array[String] cell_types = read_lines(cell_type_list)
        String bgen_pattern
        String covar_pattern
        String zfile_pattern
        String acat_pattern
        String qval_col
        String variant_id_col
        Int n_causal_snps
        File tss_bed
        Float var_y = 1.0
        Float min_cs_corr = 0
        Int chunk_size = 100
        Int window = 1000000
        String? set_variant_id_map_chr
        Boolean copy_z_files = true
        Boolean debug = false
        String zones
        String docker_colocboost
        String docker_suite

    }

    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

    scatter (cell_type in cell_types) {
        String zfile_cell_type = sub(zfile_pattern, "\\{CELL_TYPE\\}", cell_type)
        File covar = sub(covar_pattern, "\\{CELL_TYPE\\}", cell_type)
        File acat = sub(acat_pattern, "\\{CELL_TYPE\\}", cell_type)

        call precompute_projection {
            input:
                prefix = basename(zfile_cell_type, ".{CHR}.{GENE_ID}.z"),
                sample = sub(sub(bgen_pattern, "\\{CHR\\}", "chr1"), "\\.bgen$", ".sample"),
                covar = covar,
                zones = zones,
                docker = docker_colocboost
        }

        scatter (chrom in chromosomes) {
            File bgen = sub(bgen_pattern, "\\{CHR\\}", chrom)
            File bgi = "~{bgen}.bgi"
            File sample = sub(bgen, "\\.bgen$", ".sample")

            call preprocess {
                input:
                    chrom = chrom,
                    cell_type = cell_type,
                    zfile_pattern = sub(zfile_cell_type, "\\{CHR\\}", chrom),
                    acat = acat,
                    qval_col = qval_col,
                    variant_id_col = variant_id_col,
                    covar = covar,
                    chunk_size = chunk_size,
                    debug = debug,
                    docker = docker_suite,
                    zones = zones
            }

            scatter (zfile_list in preprocess.zfiles) {
                Array[File] zfiles = read_lines(zfile_list)

                call run_ldcov_susie {
                    input:
                        bgen = bgen,
                        bgi = bgi,
                        sample = sample,
                        projection = precompute_projection.out_projection,
                        zfiles = zfiles,
                        n_samples=preprocess.n_samples,
                        n_causal_snps=n_causal_snps,
                        var_y=var_y,
                        min_cs_corr=min_cs_corr,
                        tss_bed=tss_bed,
                        window=window,
                        docker=docker_colocboost,
                        zones=zones
                }

                if (copy_z_files) {
                    call copy as copy_zfiles {
                        input:
                            output_dir = "~{output_directory_stripped}/~{basename(zfile_cell_type, '.{CHR}.{GENE_ID}.z')}",
                            files = flatten([run_ldcov_susie.snp_txt, run_ldcov_susie.snp_txt_tbi]),
                            docker = docker_suite,
                            zones = zones
                    }

                    call copy as copy_bcor {
                        input:
                            output_dir = sub(sub(bcor_output_directory_pattern, "\\{CELL_TYPE\\}", cell_type), "[/\\s]+$", ""),
                            files = run_ldcov_susie.bcor,
                            docker = docker_suite,
                            zones = zones
                    }
                }

                call combine {
                    input:
                        susie_snp=run_ldcov_susie.snp,
                        susie_cred=run_ldcov_susie.cred,
                        susie_cred_99=run_ldcov_susie.cred_99,
                        pheno_chrom=basename(sub(zfile_cell_type, "\\{CHR\\}", chrom), ".{GENE_ID}.z"),
                        zones=zones,
                        docker=docker_suite
                }
            }

            call combine_all {
                input:
                    susie_snp=combine.out_susie_snp,
                    susie_cred=combine.out_susie_cred,
                    susie_cred_99=combine.out_susie_cred_99,
                    failed_zfiles=run_ldcov_susie.failed_zfiles,
                    pheno_chrom=basename(sub(zfile_cell_type, "\\{CHR\\}", chrom), ".{GENE_ID}.z"),
                    zones=zones,
                    docker=docker_suite
            }
        }
    }
}

task preprocess {
    input {
        String chrom
        String cell_type
        String zfile_pattern
        File acat
        String qval_col
        String variant_id_col
        File covar
        Int chunk_size = 100
        Boolean debug = false
        String docker
        String zones
    }

    command <<<
        set -e

        python3 << "__EOF__"
        import math
        import pandas as pd
        import sys

        acat_df = pd.read_csv("~{acat}", sep="\t")
        acat_df["chrom"] = acat_df["~{variant_id_col}"].str.split("_").str[0]
        sig_genes = acat_df[(acat_df["~{qval_col}"] < 0.05) & (acat_df["chrom"] == "~{chrom}")]["phenotype_id"].unique()

        if len(sig_genes) == 0:
            print("No significant genes found.")
            with open("n_samples.txt", "w") as f:
                f.write("0")
            with open("~{cell_type}.incl", "w") as f:
                f.write("")
            sys.exit(0)

        chunk_size = ~{chunk_size}
        zfile_paths = ["~{zfile_pattern}".replace("{GENE_ID}", x) for x in sig_genes]
        df = pd.DataFrame({'zfile_path': zfile_paths})
        total_lines = len(df.index)
        n_chunks = math.ceil(total_lines / chunk_size)
        lines_per_chunk = math.ceil(total_lines / n_chunks)

        if ~{true='True' false='False' debug}:
            n_chunks = 1

        for i in range(n_chunks):
            start_line = i * lines_per_chunk
            end_line = min(start_line + lines_per_chunk, total_lines)
            df.iloc[start_line:end_line].to_csv(f"~{chrom}.zfiles.{i}.txt", sep="\t", index=False, header=False)
        __EOF__

        zcat ~{covar} | tail -n+2 | cut -f1 > ~{cell_type}.incl
        wc -l ~{cell_type}.incl | cut -f1 -d' ' > n_samples.txt
    >>>

    output {
        Int n_samples = read_int("n_samples.txt")
        File incl = cell_type + ".incl"
        Array[File] zfiles = glob("~{chrom}.zfiles.*.txt")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "7 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task precompute_projection {
    input {
        String prefix
        File sample
        File covar
        String zones
        String docker
    }

    command <<<
        #!/usr/bin/env bash

        set -e

        ldcov \
        --sample ~{sample} \
        --covariates ~{covar} \
        --precompute-projection \
        --out ~{prefix}
    >>>

    output {
        File out_projection = "~{prefix}.proj.npz"
    }

    runtime {

        docker: "~{docker}"
        cpu: 1
        memory: "8 GB"
        disks: "local-disk 10 HDD"
        zones: "~{zones}"
        preemptible: 2
        noAddress: true
    }
}

task run_ldcov_susie {
    input {
        File bgen
        File bgi
        File sample
        File projection
        Array[File] zfiles
        Int n_samples
        Int n_causal_snps
        Float var_y
        File tss_bed
        Int window
        String zones
        String docker
        Float min_cs_corr
    }

    command <<<
        #!/usr/bin/env bash

        set -eu

        n_cpu=$(grep -c ^processor /proc/cpuinfo)
        touch failed_zfiles.txt

        for zfile in ~{sep=" " zfiles}
        do
            prefix=$(basename "$zfile" .z)
            prefix_wo_chr=$(echo "$prefix" | awk -F. 'BEGIN{OFS="."}{$((NF-1))=""; gsub(/\.+/, "."); gsub(/^\.|\.$/,""); print}')

            trait=$(echo "$prefix" | awk -F'.' '{print $NF}')
            if [[ "$trait" == ENSG* ]]
            then
                region=$(awk -v t="$trait" -v w=~{window} '$4==t {s = $2 - w + 1; if (s < 0) {s = 0}; print $1":"s"-"$3 + w}' ~{tss_bed})
            else
                region=$(echo "$trait" | awk -v w=~{window} -F'-' '{s = $2 - w + 1; if (s < 0) {s = 0}; print $1":"s"-"$3 + w}')
            fi

            echo "Processing SNP file for trait: $trait, region: $region"
            {
                ldcov \
                --bgen ~{bgen} \
                --bgi ~{bgi} \
                --sample ~{sample} \
                --projection ~{projection} \
                --z "$zfile" \
                --compute-ld \
                --output-format bcor \
                --nan-action mean \
                --out "$prefix" && \
                run_susieR.R \
                --z "$zfile" \
                --ld "$prefix.bcor" \
                -n ~{n_samples} \
                --L ~{n_causal_snps} \
                --var-y ~{var_y} \
                --snp "$prefix.susie.snp" \
                --cred "$prefix.susie.cred" \
                --log "$prefix.susie.log" \
                --susie-obj "$prefix.susie.rds" \
                --save-susie-obj \
                --write-alpha \
                --write-single-effect \
                --write-lbf-variable \
                --min-cs-corr ~{min_cs_corr} && \
                awk -v trait="$trait" -v region="$region" '
                BEGIN {
                    OFS = "\t"
                }
                NR == 1 {
                    for (i = 1; i <= NF; i++) {
                        col[$i] = i
                    }
                    gsub(" ", "\t")
                    print "trait", "region", \
                        "rsid", "chromosome", "position", \
                        "cs", "low_purity", "prob", \
                        "lbf_variable1", "lbf_variable2", "lbf_variable3", \
                        "lbf_variable4", "lbf_variable5", "lbf_variable6", \
                        "lbf_variable7", "lbf_variable8", "lbf_variable9", \
                        "lbf_variable10", "beta", "se", "p"
                }
                NR > 1 {
                    print trait, region, \
                        $col["rsid"], $col["chromosome"], $col["position"], \
                        $col["cs"], $col["low_purity"], $col["prob"], \
                        $col["lbf_variable1"], $col["lbf_variable2"], $col["lbf_variable3"], \
                        $col["lbf_variable4"], $col["lbf_variable5"], $col["lbf_variable6"], \
                        $col["lbf_variable7"], $col["lbf_variable8"], $col["lbf_variable9"], \
                        $col["lbf_variable10"], $col["beta"], $col["se"], $col["p"]
                }
                ' $prefix.susie.snp | bgzip -c -@ $n_cpu > $prefix_wo_chr.txt.gz && \
                tabix -s 4 -b 5 -e 5 -S 1 $prefix_wo_chr.txt.gz
            } || {
                exit_status=$?
                # Handle OOM/Killed (exit codes 137=KILL, 143=TERM)
                if [[ $exit_status -eq 137 || $exit_status -eq 143 ]]; then
                    echo "FATAL ERROR (OOM/Killed) processing $zfile. Exiting."
                    echo "$zfile" >> failed_zfiles.log
                    exit 1
                # Handle specific Susie error
                elif grep -q "Estimating residual variance failed: the estimated value is negative" "$prefix.susie.log"
                then
                    echo "Susie residual variance error for $zfile (logged, continuing)"
                    echo "$zfile" >> failed_zfiles.txt
                elif grep -q "The number of variants ([0-9]*) is less than current L ([0-9]*)" "$prefix.susie.log"
                then
                    echo "Variant count error for $zfile (logged, continuing)"
                    echo "$zfile" >> failed_zfiles.txt
                else
                    echo "Fatal error processing $zfile (logged, exiting)"
                    echo "$zfile" >> failed_zfiles.txt
                    exit 1
                fi
            }
        done

        echo "Processing complete."
        echo "Failed files count: $(wc -l < failed_zfiles.txt)"
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        Array[File] log = glob("*.susie.log")
        Array[File] snp = glob("*.susie.snp")
        Array[File] cred = glob("*.susie.cred")
        Array[File] cred_99 = glob("*.susie.cred_99")
        Array[File] rds = glob("*.susie.rds")
        File failed_zfiles = "failed_zfiles.txt"
        Array[File] snp_txt = glob("*.txt.gz")
        Array[File] snp_txt_tbi = glob("*.txt.gz.tbi")
        Array[File] bcor = glob("*.bcor")
    }

    runtime {

        docker: "~{docker}"
        cpu: 4
        memory: "36 GB"
        disks: "local-disk 100 HDD"
        zones: "~{zones}"
        maxRetries: 1
        preemptible: 2
        noAddress: true
    }
}

task combine {
    input {
        String pheno_chrom
        String pheno = sub(pheno_chrom, "\\.chr[0-9X]+$", "")
        Array[File] susie_snp
        Array[File] susie_cred
        Array[File] susie_cred_99
        String zones
        String docker
    }

    command <<<
        set -eux

        n_cpu=$(grep -c ^processor /proc/cpuinfo)

        cat << "__EOF__" > combine_snp.awk
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
            gsub(" ", "\t")
            print "trait", "region", "v", $0
        }
        FNR == 1 {
            match(FILENAME, /(chr[0-9X]+)\.(ENSG[0-9]+|chr[0-9X]+-[0-9]+-[0-9]+)\./, a)
            region = a[2]
        }
        FNR > 1 {
            chrom = substr($col["chromosome"], 4)
            sub(/^0/, "", chrom)
            v = sprintf( \
                "%s:%s:%s:%s", \
                chrom, \
                $col["position"], \
                $col["allele1"], \
                $col["allele2"] \
            )
            gsub(" ", "\t")
            print pheno, region, v, $0 | "sort -V -k3"
        }
        __EOF__

        # Combine susie .snp files
        awk -f combine_snp.awk -v pheno=~{pheno} ~{sep=" " susie_snp} | bgzip -c -@ $n_cpu > ~{pheno_chrom}.SUSIE.snp.bgz
        tabix -s 5 -b 6 -e 6 -S 1 ~{pheno_chrom}.SUSIE.snp.bgz

        # Combine susie .cred files
        awk -v pheno=~{pheno} '
        BEGIN {
            OFS = "\t"
            header_printed=0
        }
        FNR == 1 {
            gsub(/[ \t]+$/, "", $0)
            if(header_printed==0 && $0!="") {
                print "trait", "region", $0
                header_printed=1
            };
            match(FILENAME, /(chr[0-9X]+)\.(ENSG[0-9]+|chr[0-9X]+-[0-9]+-[0-9]+)\./, a)
            region = a[2]
        }
        FNR > 1 {
            print pheno, region, $0
        }
        ' ~{sep=" " susie_cred} | bgzip -c -@ $n_cpu > ~{pheno_chrom}.SUSIE.cred.bgz

        awk -v pheno=~{pheno} '
        BEGIN {
            OFS = "\t"
            header_printed=0
        }
        FNR == 1 {
            if(header_printed==0) {
                print "trait", "region", $0
                header_printed=1
            };
            match(FILENAME, /(chr[0-9X]+)\.(ENSG[0-9]+|chr[0-9X]+-[0-9]+-[0-9]+)\./, a)
            region = a[2]
        }
        FNR > 1 {
            print pheno, region, $0
        }
        ' ~{sep=" " susie_cred_99} | bgzip -c -@ $n_cpu > ~{pheno_chrom}.SUSIE.cred_99.bgz
    >>>

    output {
        File out_susie_snp = "~{pheno_chrom}.SUSIE.snp.bgz"
        File out_susie_snp_tbi = "~{pheno_chrom}.SUSIE.snp.bgz.tbi"
        File out_susie_cred = "~{pheno_chrom}.SUSIE.cred.bgz"
        File out_susie_cred_99 = "~{pheno_chrom}.SUSIE.cred_99.bgz"
    }

    runtime {

        docker: "~{docker}"
        cpu: "1"
        memory: "7 GB"
        disks: "local-disk 100 HDD"
        zones: "~{zones}"
        preemptible: 2
        noAddress: true
    }
}

task combine_all {
    input {
        String pheno_chrom
        String pheno = sub(pheno_chrom, "\\.chr[0-9X]+$", "")
        Array[File] susie_snp
        Array[File] susie_cred
        Array[File] susie_cred_99
        Array[File] failed_zfiles
        String zones
        String docker
    }

    command <<<
        set -eux

        n_cpu=$(grep -c ^processor /proc/cpuinfo)

        # Combine susie .snp files
        zcat ~{sep=" " susie_snp} | awk 'NR == 1 {print} $1 != "trait" {print | "sort -V -k3"}' | bgzip -c -@ $n_cpu > ~{pheno_chrom}.SUSIE.snp.bgz
        tabix -s 5 -b 6 -e 6 -S 1 ~{pheno_chrom}.SUSIE.snp.bgz

        # Filter in_cs
        zcat ~{pheno_chrom}.SUSIE.snp.bgz | awk '
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
            print
        }
        NR > 1 && $col["cs"] > 0 && $col["low_purity"] == 0 {
            print
        }
        ' | bgzip -c -@ $n_cpu > ~{pheno_chrom}.SUSIE.in_cs.snp.bgz
        tabix -s 5 -b 6 -e 6 -S 1 ~{pheno_chrom}.SUSIE.in_cs.snp.bgz

        # Combine susie .cred files
        zcat ~{sep=" " susie_cred} | awk 'NR == 1 || $1 != "trait"' | bgzip -c -@ $n_cpu > ~{pheno_chrom}.SUSIE.cred.bgz
        zcat ~{sep=" " susie_cred_99} | awk 'NR == 1 || $1 != "trait"' | bgzip -c -@ $n_cpu > ~{pheno_chrom}.SUSIE.cred_99.bgz

        cat ~{sep=" " failed_zfiles} | sort -V > ~{pheno_chrom}.failed_zfiles.txt
    >>>

    output {
        File out_susie_snp = "~{pheno_chrom}.SUSIE.snp.bgz"
        File out_susie_snp_tbi = "~{pheno_chrom}.SUSIE.snp.bgz.tbi"
        File out_susie_in_cs_snp = "~{pheno_chrom}.SUSIE.in_cs.snp.bgz"
        File out_susie_in_cs_snp_tbi = "~{pheno_chrom}.SUSIE.in_cs.snp.bgz.tbi"
        File out_susie_cred = "~{pheno_chrom}.SUSIE.cred.bgz"
        File out_susie_cred_99 = "~{pheno_chrom}.SUSIE.cred_99.bgz"
        File out_failed_zfiles = "~{pheno_chrom}.failed_zfiles.txt"
    }

    runtime {

        docker: "~{docker}"
        cpu: "1"
        memory: "7 GB"
        disks: "local-disk 100 HDD"
        zones: "~{zones}"
        preemptible: 2
        noAddress: true
    }
}

task copy {
    input {
        String output_dir
        Array[String] files
        String docker
        String zones
    }

    command <<<
        set -e

        python3 << "__EOF__"
        #!/usr/bin/env python3
        files = "~{sep=',' files}".split(",")

        with open("files.txt", "w") as f:
            for file in files:
                f.write(f"{file}\n")
        __EOF__

        cat files.txt | gcloud storage cp -I ~{output_dir}/ && \
        touch _SUCCESS
    >>>

    meta {
        volatile: true
    }

    output {
        File out_success = "_SUCCESS"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "7 GB"
        disks: "local-disk 100 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}
