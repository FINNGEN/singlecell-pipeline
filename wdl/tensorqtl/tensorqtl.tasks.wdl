version 1.0

task peer_factor {
    input {
        File psam
        File bed
        String prefix
        File covar
        Array[String] covariates
        String? covar_sample_id
        File? subset_samples
        Int? K
        File? n_peer
        String? cell_type
        Int max_iterations = 10000
        String docker
        String zones
    }

    command <<<
        set -e

        K=""
        if ~{true='true' false='false' defined(n_peer) && defined(cell_type)}
        then
            K=$(awk -v t="~{cell_type}" '$1 == t{if ($2 != "NA"){print "--K "$2} else {print ""}}' ~{n_peer})
            echo "n_peer file sets $K"
        fi
        if ~{true='true' false='false' defined(K)}
        then
            K="--K ~{K}"
            echo "K input sets $K"
        fi

        compute_peer_factor.R \
        --psam ~{psam} \
        --expression-bed ~{bed} \
        --covar ~{covar} \
        --covariates ~{sep=" " covariates} \
        ~{true='--covar-sample-id ' false='' defined(covar_sample_id)}~{covar_sample_id} \
        ~{true='--subset-samples ' false='' defined(subset_samples)}~{subset_samples} \
        $K \
        --max-iterations ~{max_iterations} \
        --prefix ~{prefix}

        zcat ~{prefix}.peer_alpha.tsv.gz | tail -n+2 | wc -l > ~{prefix}.K.txt
        zcat ~{prefix}.peer_factors_covariates.tsv.gz | tail -n+2 | wc -l > n_samples.txt
    >>>

    output {
        File out_K = prefix + ".K.txt"
        File out_alpha = prefix + ".peer_alpha.tsv.gz"
        File out_covar = prefix + ".peer_factors_covariates.tsv.gz"
        File out_bed = prefix + ".peer_residuals.bed.gz"
        Int n_samples = read_int("n_samples.txt")
    }

    runtime {
        docker: docker
        cpu: 8
        memory: "30 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task compute_loco_expression_pcs {
    input {
        String? chrom
        File psam
        File bed
        String prefix
        File covar
        Array[String] covariates
        String? covar_sample_id
        Int? K
        File? n_peer
        String? cell_type
        String docker
        String zones
    }

    Int task_size = ceil(20 * size(bed, 'GB')) + 36
    Int memory = if (task_size < 8) then 8 else task_size

    command <<<
        set -e

        K=""
        if ~{true='true' false='false' defined(n_peer) && defined(cell_type)}
        then
            K=$(awk -v t="~{cell_type}" '$1 == t{if ($2 != "NA"){print "--K "$2} else {print ""}}' ~{n_peer})
            echo "n_peer file sets $K"
        fi
        if ~{true='true' false='false' defined(K)}
        then
            K="--K ~{K}"
            echo "K input sets $K"
        fi

        compute_loco_expression_pcs.R \
        --psam ~{psam} \
        ~{true='--chromosome ' false='' defined(chrom)}~{chrom} \
        --expression-bed ~{bed} \
        --covar ~{covar} \
        --covariates ~{sep=" " covariates} \
        ~{true='--covar-sample-id ' false='' defined(covar_sample_id)}~{covar_sample_id} \
        $K \
        --prefix ~{prefix}

        zcat ~{prefix}.*gex_pcs_covariates.tsv.gz | tail -n+2 | wc -l > n_samples.txt
    >>>

    output {
        File out_covar = if defined(chrom) then "~{prefix}.loco_~{chrom}.gex_pcs_covariates.tsv.gz" else "~{prefix}.gex_pcs_covariates.tsv.gz"
        Int n_samples = read_int("n_samples.txt")
    }

    runtime {
        docker: docker
        cpu: 8
        memory: "~{memory} GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task map {
    input {
        String mode
        String pfile
        File psam = pfile + ".psam"
        File pgen = pfile + ".pgen"
        File pvar = pfile + ".pvar"
        File bed
        File covar
        File? interact
        String chrom
        Int window = 1000000
        Float maf_threshold = 0.01
        Int mac_threshold = 0
        Boolean dominant = false
        Boolean recessive = false
        Boolean supress_zfiles = false
        String docker
        String zones
    }

    String prefix = if mode == "cis" then basename(bed, ".bed.gz") else basename(bed, ".bed.gz") + "." + chrom
    String write_zfiles = if (mode == "cis_nominal" && !supress_zfiles) then "--write-z-files" else ""
    String acat = if mode == "cis_nominal" then "--acat" else ""
    String interaction = if mode == "cis_interaction" then "--interaction ~{interact}" else ""

    command <<<
        set -e

        if [[ "~{sub(pgen, '.pgen$', '.psam')}" != "~{psam}" ]]
        then
            mv "~{psam}" "~{sub(pgen, '.pgen$', '.psam')}"
        fi
        if [[ "~{sub(pgen, '.pgen$', '.pvar')}" != "~{pvar}" ]]
        then
            mv "~{pvar}" "~{sub(pgen, '.pgen$', '.pvar')}"
        fi

        run_tensorqtl.py \
        --mode ~{mode} \
        --pfile ~{sub(pgen, ".pgen$", "")} \
        --expression-bed ~{bed} \
        --covar ~{covar} \
        ~{interaction} \
        --chrom ~{chrom} \
        --window ~{window} \
        --maf-threshold ~{maf_threshold} \
        --mac-threshold ~{mac_threshold} \
        ~{true='--dominant ' false='' dominant} \
        ~{true='--recessive ' false='' recessive} \
        --write-tsv \
        ~{write_zfiles} \
        ~{acat} \
        --prefix ~{prefix}

        if ~{true='true' false='false' defined(interact)}
        then
            zcat ~{prefix}.cis_qtl_top_assoc.txt.gz | awk '
            NR == 1 {
                print
            }
            NR > 1 {
                print | "sort -V -k2,2"
            }
            ' | gzip -c > ~{prefix}.cis_qtl_top_assoc.tsv.gz
        fi

        if ~{true='true' false='false' mode != "cis_nominal"}
        then
            touch ~{prefix}.cis_qtl_acat.tsv.gz
        fi

    >>>

    output {
        File out = prefix + (if mode == "cis_nominal" || mode == "cis_interaction" then ".cis_qtl_pairs.tsv.gz" else (if mode == "cis" then ".cis_qtl_egenes.tsv.gz" else ".trans_qtl_pairs.tsv.gz"))
        File out_acat = "~{prefix}.cis_qtl_acat.tsv.gz"
        Array[File] out_interaction_top = glob("*.cis_qtl_top_assoc.tsv.gz")
        Array[File] zfiles = glob("*.z")
    }

    runtime {
        docker: docker
        cpu: 16
        memory: (ceil(0.08 * size(bed, 'MB')) + (if mode == "trans" then 96 else 24)) + " GB"
        disks: "local-disk 20 HDD"
        zones: zones
        preemptible: 2
        maxRetries: 1
        noAddress: true
    }
}

task calculate_qvalues {
    input {
        Array[File] cis_df
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e

        calculate_qvalues.py \
        --cis-df ~{sep=" " cis_df} \
        --prefix ~{prefix}

    >>>

    output {
        File out = "~{prefix}.cis_qtl_egenes.tsv.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "16 GB"
        disks: "local-disk ~{2 * ceil(size(cis_df, 'GB')) + 10} HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task calculate_qvalues_acat {
    input {
        Array[File] cis_df
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e

        cat << "__EOF__" > qvalue.R
        library(dplyr)

        paths = unlist(stringr::str_split("~{sep="," cis_df}", ","))
        df = purrr::map_dfr(paths, function(path) {
            data.table::fread(path, data.table=FALSE)
        }) %>%
            dplyr::mutate(qval = qvalue::qvalue(pval_acat)$qvalue)

        write.table(df, "~{prefix}.cis_qtl_acat.tsv", quote = F, row.names = F, sep = "\t")
        __EOF__

        Rscript qvalue.R && \
        bgzip ~{prefix}.cis_qtl_acat.tsv
    >>>

    output {
        File out = "~{prefix}.cis_qtl_acat.tsv.gz"
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

task combine_tsv {
    input {
        String mode
        Array[File] tsv_paths
        String cell_type
        String sort_cmd = if mode == "trans" then '| "sort -V -k1,2"' else ""
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e

        n_cpu=$(grep -c ^processor /proc/cpuinfo)

        zcat ~{sep=" " tsv_paths} | awk -v cell_type="~{cell_type}" '
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
            print "#chrom", "position", "cell_type", $0
            next
        }
        $0 !~ /variant_id/ {
            split($col["variant_id"], a, /[_:]/)
            chrom = a[1]
            pos = a[2]
            print chrom, pos, cell_type, $0 ~{sort_cmd}
        }
        ' | bgzip -@${n_cpu} -c > ~{prefix}.tsv.gz

        tabix -@${n_cpu} -s 1 -b 2 -e 2 -S 1 ~{prefix}.tsv.gz
    >>>

    output {
        File out_tsv = prefix + ".tsv.gz"
        File out_tsv_tbi = prefix + ".tsv.gz.tbi"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "7 GB"
        disks: "local-disk 50 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task import_hail {
    input {
        Array[String] tsv_paths
        String out_ht
        Boolean overwrite
        String model = "additive"
        Int min_partitions = 20000
        String docker
        String zones
    }

    Boolean is_additive = model == "additive"

    command <<<
        set -e

        mem=$(grep MemTotal /proc/meminfo | awk '{printf("%dg\n", $2 / 1024 / 1024)}')

        echo $mem
        df -h

        if ~{true='true' false='false' is_additive}
        then
            echo '~{sep=" " tsv_paths}' | tr " " "\n" | grep -vF ".dominant." | grep -vF ".recessive." > tsv_paths.txt
        else
            echo '~{sep=" " tsv_paths}' | tr " " "\n" | grep -F ".~{model}." > tsv_paths.txt
        fi

        cat << "__EOF__" > script.py
        import pandas as pd
        import uuid
        import hail as hl
        from hail.utils import new_temp_file

        paths = pd.read_csv("tsv_paths.txt", header=None).iloc[:, 0].to_list()
        print(f"Importing {len(paths)} tsv files")
        ht = hl.import_table(paths, force_bgz=True, impute=True, min_partitions=~{min_partitions})
        ht = ht.annotate(**hl.parse_variant(ht.variant_id.replace("_", ":"), reference_genome="GRCh38"))
        ht = ht.key_by("locus", "alleles", "phenotype_id", "cell_type")
        ht = ht.drop("variant_id", "#chrom", "position")
        ht = ht.checkpoint("~{out_ht}", overwrite=~{true='True' false='False' overwrite})
        __EOF__

        PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=$mem pyspark-shell" \
        python3 script.py && \
        touch _SUCCESS

        df -h
    >>>

    output {
        File success = "_SUCCESS"
    }

    runtime {
        docker: docker
        cpu: 64
        memory: "512 GB"
        disks: "local-disk 400 HDD"
        zones: zones
        preemptible: 2
    }
}

task count_n_samples {
    input {
        File bed
        String docker
        String zones
    }

    command <<<
        set -e

        zcat ~{bed} | head -n1 | tr '\t' '\n' | tail -n+5 | wc -l > n_samples.txt
    >>>

    output {
        Int n_samples = read_int("n_samples.txt")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}

task count_egenes {
    input {
        Array[String] result_paths
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e

        mem=$(grep MemTotal /proc/meminfo | awk '{printf("%dg\n", $2 / 1024 / 1024)}')
        echo "$mem"

        cat << "__EOF__" > script.py
        import pandas as pd
        import hail as hl
        from hail.methods.statgen import _lambda_gc_agg
        from hail.utils import new_temp_file

        paths = "~{sep=',' result_paths}".split(",")
        ht = hl.import_table(paths, force_bgz=True, source_file_field="source", impute=True)
        ht = ht.annotate(
            cell_type=ht.source.replace(r"^.*\.fgid(?:\.sct\.counts)?(?:\.qc)?\.(.*)\.[^\.]*\.inv.*$", "$1"),
            K=hl.int32(ht.source.replace(r"^.*\.inv\.K([0-9]+)\..*$", "$1")),
        )
        ht = ht.select("phenotype_id", "cell_type", "K", "qval")
        ht = ht.checkpoint(new_temp_file())

        ht_agg = ht.group_by("cell_type", "K").aggregate(
            n_genes=hl.agg.count(),
            n_egenes1=hl.agg.count_where(ht["qval"] < 0.5),
            n_egenes2=hl.agg.count_where(ht["qval"] < 0.1),
            n_egenes3=hl.agg.count_where(ht["qval"] < 0.05),
            n_egenes4=hl.agg.count_where(ht["qval"] < 0.01),
            n_egenes5=hl.agg.count_where(ht["qval"] < 0.001),
            n_egenes6=hl.agg.count_where(ht["qval"] < 0.0001),
            n_egenes7=hl.agg.count_where(ht["qval"] < 0.00001),
        )
        ht_agg.export("~{prefix}.n_egenes.txt")
        __EOF__

        PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=$mem pyspark-shell" \
        python3 script.py && \
        touch _SUCCESS
    >>>

    output {
        File out_txt = "~{prefix}.n_egenes.txt"
        File success = "_SUCCESS"
    }

    runtime {
        docker: docker
        cpu: 16
        memory: "32 GB"
        disks: "local-disk 50 HDD"
        zones: zones
        preemptible: 2
    }
}

task optimize_n_peer {
    input {
        Array[File] result_paths
        String n_egenes_column = "n_egenes3"
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e

        cat << "__EOF__" > script.R
        library(dplyr)
        # https://github.com/powellgenomicslab/PEER_factors/blob/main/8-Supp_Figures_Tables/Sup_Figure8_local_greedy_PFs.R

        find_local_greedy_n_peer <- function(data, n_egenes_column = "~{n_egenes_column}") {
            if (nrow(data) <= 1) {
                return(data)
            }

            # Calculate percent change
            data$percent <-
                c(NA, diff(data[[n_egenes_column]]) / data[[n_egenes_column]][-nrow(data)])

            # Skip smoothing if insufficient data points
            if (nrow(data) < 4) {
                data$smooth_y2 <- NA
                data$local_greedy <- NA
                return(data)
            }

            # Perform smoothing using smooth.spline with cross-validation
            data$smooth_y2 <- NA

            # Subset data to remove NA values (first row has NA percent)
            valid_data <- data[2:nrow(data), ]
            valid_data <- valid_data[!is.na(valid_data$percent), ]

            if (nrow(valid_data) < 4) {
                # Not enough valid data points for smoothing
                data$local_greedy <- NA
                return(data)
            }

            # Try to fit smooth.spline with cross-validation
            spline_fit <- tryCatch({
                # cv=TRUE enables cross-validation for automatic smoothing parameter selection
                smooth.spline(x = valid_data$K, y = valid_data$percent, cv = TRUE)
            }, error = function(e) {
                warning(sprintf(
                "Smoothing failed for cell type %s: %s",
                data$cell_type[1],
                e$message
                ))
                return(NULL)
            })

            if (!is.null(spline_fit)) {
                # Predict smoothed values for all K values (except first row)
                pred_values <- predict(spline_fit, x = data$K[2:nrow(data)])
                data$smooth_y2[2:nrow(data)] <- pred_values$y

                # Find local greedy point
                negative_idx <- which(data$smooth_y2 < 0)
                data$local_greedy <- if (length(negative_idx) > 0) {
                    data$K[min(negative_idx) - 1]
                } else {
                    NA
                }
            } else {
                data$local_greedy <- NA
            }

            return(data)
        }

        paths <- unlist(stringr::str_split("~{sep=',' result_paths}", ","))

        df <- purrr::map_dfr(paths, function(x) {data.table::fread(x, data.table = FALSE)}) %>%
            dplyr::arrange(cell_type, K) %>%
            dplyr::group_split(cell_type) %>%
            purrr::map_dfr(find_local_greedy_n_peer)

        df.n_peer <- dplyr::distinct(df, cell_type, local_greedy)

        write.table(df, "~{prefix}.n_egenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
        write.table(df.n_peer, "~{prefix}.n_peer.txt", sep = "\t", quote = FALSE, row.names = FALSE)
        __EOF__

        Rscript script.R && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_n_egenes = "~{prefix}.n_egenes.txt"
        File out_n_peer = "~{prefix}.n_peer.txt"
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