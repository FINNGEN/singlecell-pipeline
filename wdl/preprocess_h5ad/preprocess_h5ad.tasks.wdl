version 1.0

task update_h5ad {
    input {
        File h5ad
        String type
        String prefix = basename(h5ad, ".h5ad") + ".fgid"
        File? features
        Array[File]? obs_annot
        Array[String]? obs_key
        Array[String]? obs_annot_key
        String obs_sample_id
        String? obs_pool_id
        Array[String]? drop_obs_columns
        Array[String]? drop_var_columns
        Array[String]? sanitize_cell_type_annots
        Boolean update_barcodes
        String? resequenced_suffix
        File? keep_samples
        Float? downsample
        Array[String]? vars_to_regress
        Boolean filter_cells = true
        Int min_genes = 200
        Int? min_counts
        Int? max_counts
        Int min_cells = 400
        Boolean filter_samples = true
        Int min_cells_per_sample = 3000
        Boolean calculate_qc_metrics = true
        String docker
        String zones
    }

    command <<<
        set -e

        update_h5ad.py \
        --prefix "~{prefix}" \
        --h5ad "~{h5ad}" \
        --type "~{type}" \
        ~{true='--features ' false='' defined(features)}~{features} \
        ~{true='--obs-annot ' false='' defined(obs_annot)}~{sep=" " obs_annot} \
        ~{true='--obs-key ' false='' defined(obs_key)}~{sep=" " obs_key} \
        ~{true='--obs-annot-key ' false='' defined(obs_annot_key)}~{sep=" " obs_annot_key} \
        --obs-sample-id ~{obs_sample_id} \
        ~{true='--obs-pool-id ' false='' defined(obs_pool_id)}~{sep=" " obs_pool_id} \
        ~{true='--drop-obs-columns ' false='' defined(drop_obs_columns)}~{sep=" " drop_obs_columns} \
        ~{true='--drop-var-columns ' false='' defined(drop_var_columns)}~{sep=" " drop_var_columns} \
        ~{true='--sanitize-cell-type-annots ' false='' defined(sanitize_cell_type_annots)}~{sep=" " sanitize_cell_type_annots} \
        ~{true='--update-barcodes' false='' update_barcodes} \
        ~{true='--resequenced-suffix ' false='' defined(resequenced_suffix)}~{resequenced_suffix} \
        ~{true='--keep-samples ' false='' defined(keep_samples)}~{keep_samples} \
        ~{true='--downsample ' false='' defined(downsample)}~{downsample} \
        ~{true='--filter-cells ' false='' filter_cells} \
        --min-genes ~{min_genes} \
        ~{true='--min-counts ' false='' defined(min_counts)}~{min_counts} \
        ~{true='--max-counts ' false='' defined(max_counts)}~{max_counts} \
        --min-cells ~{min_cells} \
        ~{true='--filter-samples ' false='' filter_samples} \
        --min-cells-per-sample ~{min_cells_per_sample} \
        ~{true='--vars-to-regress ' false='' defined(vars_to_regress)}~{sep=" " vars_to_regress} \
        ~{true='--calculate-qc-metrics' false='' calculate_qc_metrics}

    >>>

    output {
        File out_h5ad = prefix + ".h5ad"
        File out_obs = prefix + ".obs.tsv.gz"
        File out_var = prefix + ".var.tsv.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "624 GB"
        disks: "local-disk " + ceil(2.5 * size(h5ad, 'GB')) + " HDD"
        zones: zones
        preemptible: 2
    }
}

task convert_to_bpcells {
    input {
        File h5ad
        File obs
        String batch_var
        String prefix = basename(h5ad, ".h5ad")
        String docker
        String zones
    }

    command <<<
        set -e

        echo $TMPDIR

        cat << "__EOF__" > script.R
        library(dplyr)

        # convert to BPCells matrix
        mat <- BPCells::open_matrix_anndata_hdf5("~{h5ad}")
        # ensure unsigned integer type
        mat <- BPCells::convert_matrix_type(mat, "uint32_t")
        mat <- BPCells::write_matrix_dir(mat, "~{prefix}.matrix_dir")

        BPCells::write_matrix_hdf5(mat, "~{prefix}.h5", "counts")

        # convert to Seurat object
        metadata <- tibble::tibble(barcode = colnames(mat)) %>%
            dplyr::left_join(data.table::fread("~{obs}", data.table = F), by = "barcode") %>%
            tibble::column_to_rownames("barcode")
        stopifnot(all(colnames(mat) == metadata$barcode))

        median_umi <- median(metadata$total_counts)
        write.table(median_umi, "median_umi.txt", quote = F, col.names = F, row.names = F)

        batches <- unique(metadata[["~{batch_var}"]])
        write.table(batches, "batches.txt", quote = F, col.names = F, row.names = F)
        __EOF__

        Rscript script.R && \
        touch _SUCCESS
    >>>

    output {
        File success = "_SUCCESS"
        File out_h5 = "~{prefix}.h5"
        Float out_median_umi = read_float("median_umi.txt")
        Array[String] out_batches = read_lines("batches.txt")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "32 GB"
        disks: "local-disk " + ceil(2 * size(h5ad, 'GB')) + " HDD"
        zones: zones
        preemptible: if size(h5ad, 'GB') > 100 then 0 else 2
    }
}

task sctransform_by_batch_theta {
    input {
        File h5
        File obs
        String batch_var
        String batch
        Int min_cells = 5
        String malt1_id = "ENSG00000251562"
        String prefix = "~{basename(h5, '.h5')}.~{batch}"
        String docker
        String zones
    }

    command <<<
        set -e

        echo $TMPDIR

        cat << "__EOF__" > sctransform.R
        library(dplyr)
        options(future.globals.maxSize = 6000 * 1024^2)

        align_rows <- function(mat, new_rownames) {
            missing_rownames <- setdiff(new_rownames, rownames(mat))
            if (length(missing_rownames) > 0) {
                mat <- rbind(
                    mat,
                    Matrix::Matrix(0, nrow = length(missing_rownames), ncol = NCOL(mat), dimnames = list(missing_rownames, colnames(mat)))
                )
            }
            mat <- mat[new_rownames, , drop = FALSE]
            return(mat)
        }

        mat <- BPCells::open_matrix_hdf5("~{h5}", "counts")

        # convert to Seurat object
        metadata <- tibble::tibble(barcode = colnames(mat)) %>%
            dplyr::left_join(data.table::fread("~{obs}", data.table = F), by = "barcode") %>%
            dplyr::filter(!!rlang::sym("~{batch_var}") == "~{batch}") %>%
            tibble::column_to_rownames("barcode")

        # filter to batch
        mat <- BPCells::write_matrix_dir(mat[, rownames(metadata), drop=FALSE], "~{prefix}.matrix_dir")
        stopifnot(all(colnames(mat) == metadata$barcode))
        obj <- Seurat::CreateSeuratObject(count = mat, meta.data = metadata)

        ncells <- dim(obj)[2]

        # run SCTransform
        obj <- Seurat::SCTransform(
            obj,
            ncells = ncells,
            min_cells = ~{min_cells},
            return.only.var.genes = FALSE
        )

        # save VST parameters
        vst.out <- Seurat:::SCTModel_to_vst(obj[["SCT"]]@SCTModel.list[[1]])

        log_umi <- log(obj$nCount_RNA)
        theta <- vst.out$model_pars_fit[, "theta"]
        alpha <- vst.out$model_pars_fit[, "(Intercept)"]
        mu <- exp(alpha + mean(log_umi))
        sd <- sqrt(mu + mu^2 / theta)

        model_pars_fit <- tibble::tibble(
            `~{batch_var}` = "~{batch}",
            gene_id = rownames(vst.out$model_pars_fit),
            theta = theta,
            intercept = alpha,
            mu = mu,
            sd = sd
        )
        write.table(model_pars_fit, "~{prefix}.model_pars_fit.txt", quote = F, row.names = F, sep = "\t")

        # MALAT1 thresholding
        library(ggplot2)
        source("https://raw.githubusercontent.com/BaderLab/MALAT1_threshold/880623e42e684f30972863999e338384a8b3e381/malat1_function.R")
        norm_counts <- obj[["SCT"]]$data["~{malt1_id}", ]
        res <- define_malat1_threshold(norm_counts, return_plots = TRUE)
        idx <- norm_counts > res$threshold
        malt1_qc_barcodes <- names(norm_counts)[idx]

        write.table(malt1_qc_barcodes, "~{prefix}.malat1_qc.barcodes.txt", quote = F, row.names = F, col.names = F)
        ggplot2::ggsave(
            filename = "~{prefix}.malat1_qc.png",
            plot = res$plots,
            device = "png",
            width = 7.2,
            height = 7.2,
            dpi = 300
        )

        # write a temporary h5 file for merge
        counts <- align_rows(obj[["SCT"]]$counts[, idx, drop=FALSE], rownames(mat))
        BPCells::write_matrix_hdf5(
            mat = BPCells::convert_matrix_type(counts, "uint32_t"),
            "~{prefix}.sct.counts.h5",
            "counts"
        )

        scale.data <- align_rows(obj[["SCT"]]$scale.data[, idx, drop=FALSE], rownames(mat))
        BPCells::write_matrix_hdf5(
            mat = BPCells::convert_matrix_type(scale.data, "float"),
            "~{prefix}.sct.scaled.h5",
            "scale.data"
        )
        __EOF__

        Rscript sctransform.R && \
        touch _SUCCESS
    >>>

    output {
        File success = "_SUCCESS"
        File out_model_pars_fit = "~{prefix}.model_pars_fit.txt"
        File out_malat1_qc_barcodes = "~{prefix}.malat1_qc.barcodes.txt"
        File out_malat1_qc_plot = "~{prefix}.malat1_qc.png"
        File out_sct_counts_h5 = "~{prefix}.sct.counts.h5"
        File out_sct_scaled_h5 = "~{prefix}.sct.scaled.h5"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: ceil(5 * size(h5, 'GB') + 10) + " GB"
        disks: "local-disk " + ceil(2 * size(h5, 'GB') + 10) + " HDD"
        zones: zones
        preemptible: 2
    }
}

task chunk_chrom {
    input {
        File h5ad
        String type
        String chrom
        String prefix = basename(h5ad, ".h5ad")
        String docker
        String zones
    }

    command <<<
        set -e

        chunk_h5ad.py \
        --h5ad "~{h5ad}" \
        --type "~{type}" \
        --chunk-by "chromosome" \
        --chromosomes "~{chrom}" \
        --prefix "~{prefix}"
    >>>

    output {
        File out_h5ad = prefix + "." + chrom + ".h5ad"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{ceil(2 * size(h5ad, 'GB')) + 10} GB"
        disks: "local-disk " + ceil(1.5 * size(h5ad, 'GB')) + " HDD"
        zones: zones
        preemptible: 2
    }
}

task apply_malat1_qc {
    input {
        File h5ad
        Array[File] malt1_qc_barcodes
        String prefix = basename(h5ad, ".h5ad")
        String docker
        String zones
    }

    command <<<
        set -e

        echo $TMPDIR

        cat ~{sep=" " malt1_qc_barcodes} > malat1_qc_barcodes.txt

        cat << "__EOF__" > save_anndata.py
        import anndata as ad
        import bgzip
        import pandas as pd


        def write_bgz(df, out_file, **kwargs):
            with open(out_file, "wb") as raw:
                with bgzip.BGZipWriter(raw) as fh:
                    fh.write(df.to_csv(**kwargs).encode("utf-8"))

        adata = ad.read_h5ad("~{h5ad}", backed="r")
        malat1_qc_barcodes = pd.read_csv("malat1_qc_barcodes.txt", header=None).iloc[:, 0].tolist()

        adata = adata[adata.obs.index.isin(malat1_qc_barcodes), :]
        adata.write_h5ad("~{prefix}.qc.h5ad")
        write_bgz(adata.obs, "~{prefix}.qc.obs.tsv.gz", index_label="barcode", sep="\t", na_rep="NA")
        write_bgz(adata.var, "~{prefix}.qc.var.tsv.gz", index_label="gene_id", sep="\t", na_rep="NA")
        __EOF__

        python3 save_anndata.py && \
        touch _SUCCESS
    >>>

    output {
        File success = "_SUCCESS"
        File out_h5ad = "~{prefix}.qc.h5ad"
        File out_obs = "~{prefix}.qc.obs.tsv.gz"
        File out_var = "~{prefix}.qc.var.tsv.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "384 GB"
        disks: "local-disk ~{ceil(2 * size(h5ad, 'GB')) + 10} HDD"
        zones: zones
        preemptible: 2
    }
}

task combine_sctransform_outputs {
    input {
        File obs
        Array[File] model_pars_fit
        String batch_var
        String prefix = basename(obs, ".qc.obs.tsv.gz")
        String docker
        String zones
    }

    command <<<
        set -e

        echo $TMPDIR

        awk 'NR == 1 || FNR > 1' ~{sep=" " model_pars_fit} | bgzip -c > ~{prefix}.sct.model_pars_fit.txt.gz

        cat << "__EOF__" > script.R
        library(dplyr)
        df.obs = data.table::fread("~{obs}", data.table = FALSE)
        df.model = data.table::fread("~{prefix}.sct.model_pars_fit.txt.gz", data.table = FALSE)

        df.lm =
            dplyr::left_join(
                df.model,
                dplyr::group_by(df.obs, pool_name) %>%
                dplyr::summarize(mean_log_umi = mean(log(total_counts))),
                by = "~{batch_var}"
            ) %>%
            dplyr::group_split(gene_id) %>%
            purrr::map_dfr(function(data) {
                model = lm(sd ~ mean_log_umi, data = data)
                coefs = coef(summary(model))
                tibble::tibble(
                    gene_id = data$gene_id[1],
                    slope = coefs["mean_log_umi", 1],
                    intercept = coefs["(Intercept)", 1],
                    p = coefs["mean_log_umi", 4],
                )
            })
        write.table(df.lm, "~{prefix}.sct.model_pars_fit.lm.txt", quote = F, row.names = F, sep = "\t")
        __EOF__

        Rscript script.R && \
        touch _SUCCESS
    >>>

    output {
        File success = "_SUCCESS"
        File out_model_pars_fit = "~{prefix}.sct.model_pars_fit.txt.gz"
        File out_model_pars_fit_lm = "~{prefix}.sct.model_pars_fit.lm.txt"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "16 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}

task chunk_gene_id {
    input {
        File h5ad
        String type
        Array[String] gene_ids
        File? tss_bed
        Boolean omit_barcode = false
        String prefix = basename(h5ad, ".h5ad")
        String docker
        String zones
    }

    command <<<
        set -e

        chunk_h5ad.py \
        --h5ad "~{h5ad}" \
        --type "~{type}" \
        --chunk-by "gene_id" \
        --gene-ids ~{sep=" " gene_ids} \
        ~{true='--tss-bed ' false='' defined(tss_bed)}~{tss_bed} \
        ~{true='--omit-barcode ' false='' omit_barcode} \
        --prefix "~{prefix}"
    >>>

    output {
        Array[File] out_tsv_gz = glob("*.tsv.gz")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{ceil(4 * size(h5ad, 'GB')) + 10} GB"
        disks: "local-disk " + ceil(2 * size(h5ad, 'GB') + 10) + " HDD"
        zones: zones
        preemptible: 2
    }
}

task merge_atac_obs {
    input {
        Array[File] obs
        String prefix
        String docker
        String zones
    }

    Int disk_space = 2 * ceil(size(obs, 'GB')) + 10

    command <<<
        set -e

        touch obs.txt
        mkdir -p obs
        for i in ~{sep=" " obs}
        do
            zcat $i | awk '
            BEGIN {
                FS = "\t"
                OFS = "\t"
            }
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    col[$i] = i
                }
                print "barcode", "n_genes", "n_counts", "n_genes_by_counts", "total_counts", "total_counts_wo_highly_expressed"
            }
            NR > 1 {
                print $col["barcode"], $col["n_genes"], $col["n_counts"], $col["n_genes_by_counts"], $col["total_counts"], $col["total_counts_wo_highly_expressed"]
            }
            ' > obs/$(basename $i)
            echo obs/$(basename $i) >> obs.txt
        done

        cat << "__EOF__" > script.R
        library(dplyr)

        obs_paths <- read.table("obs.txt", header = F, stringsAsFactors = F)[,1]
        obs <- purrr::map_dfr(obs_paths, function(x) {data.table::fread(x, data.table = F)}) %>%
            dplyr::group_by(barcode) %>%
            dplyr::summarize(n_genes = sum(n_genes),
                             n_counts = sum(n_counts),
                             n_genes_by_counts = sum(n_genes_by_counts),
                             total_counts = sum(total_counts),
                             total_counts_wo_highly_expressed = sum(total_counts_wo_highly_expressed))
        df <- data.table::fread("~{obs[1]}", data.table = F) %>%
            dplyr::select(-n_genes, -n_counts, -n_genes_by_counts, -total_counts, -total_counts_wo_highly_expressed) %>%
            dplyr::left_join(obs, by = "barcode")
        write.table(df, "~{prefix}.obs.tsv", sep = "\t", quote = F, row.names = F)
        __EOF__

        Rscript script.R && \
        bgzip ~{prefix}.obs.tsv && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_obs = "~{prefix}.obs.tsv.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "32 GB"
        disks: "local-disk ~{disk_space} HDD"
        zones: zones
        preemptible: 2
    }
}

task update_atac_h5ad {
    input {
        File h5ad
        File obs
        String prefix = basename(h5ad, ".h5ad")
        String docker
        String zones
    }

    command <<<
        set -e

        cat << "__EOF__" > script.py
        import anndata as ad
        import bgzip
        import pandas as pd

        def write_bgz(df, out_file, **kwargs):
            with open(out_file, "wb") as raw:
                with bgzip.BGZipWriter(raw) as fh:
                    fh.write(df.to_csv(**kwargs).encode("utf-8"))

        adata = ad.read_h5ad("~{h5ad}", backed="r")
        obs = pd.read_csv("~{obs}", sep="\t", index_col=0)
        keep_columns = list(set(adata.obs.columns) - set(obs.columns))
        adata.obs = adata.obs[keep_columns].join(obs, how="left")

        adata.strings_to_categoricals()
        adata.write_h5ad("~{prefix}.h5ad")
        write_bgz(adata.obs, "~{prefix}.obs.tsv.gz", index_label="barcode", sep="\t", na_rep="NA")
        write_bgz(adata.var, "~{prefix}.var.tsv.gz", index_label="gene_id", sep="\t", na_rep="NA")
        __EOF__

        python3 script.py && \
        touch _SUCCESS
    >>>

    output {
        File success = "_SUCCESS"
        File out_h5ad = prefix + ".h5ad"
        File out_obs = prefix + ".obs.tsv.gz"
        File out_var = prefix + ".var.tsv.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "512 GB"
        disks: "local-disk 600 HDD"
        zones: zones
        preemptible: 2
    }
}

task subset_h5ad {
    input {
        File h5ad
        Array[File] barcode_files
        String prefix = basename(h5ad, '.h5ad')
        Float high_expression_max_frac = 0.05
        Boolean recompute_metrics = true
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        cat << "__EOF__" > script.py
        import anndata as ad
        import bgzip
        import numpy as np
        import os
        import os.path
        import pandas as pd
        import scanpy as sc
        import shutil
        from numba import njit


        def write_bgz(df, out_file, **kwargs):
            with open(out_file, "wb") as raw:
                with bgzip.BGZipWriter(raw) as fh:
                    fh.write(df.to_csv(**kwargs).encode("utf-8"))


        @njit
        def calculate_highly_expressed_genes(data, indices, indptr, total_counts, max_fraction, n_genes):
            gene_exceed_counts = np.zeros(n_genes, dtype=np.int32)

            for gene_idx in range(n_genes):
                start = indptr[gene_idx]
                end = indptr[gene_idx + 1]

                for data_idx in range(start, end):
                    cell_idx = indices[data_idx]
                    cell_total = total_counts[cell_idx]
                    if cell_total == 0:
                        continue

                    threshold = cell_total * max_fraction
                    if data[data_idx] > threshold:
                        gene_exceed_counts[gene_idx] += 1
                        break

            return gene_exceed_counts == 0


        barcode_paths = "~{sep=',' barcode_files}".split(",")
        for path in barcode_paths:
            print(f"Processing {path}...")
            barcodes = pd.read_csv(path, header=None).iloc[:, 0].tolist()
            barcode_prefix = os.path.basename(path).removesuffix(".txt.gz")
            fname = f"~{prefix}.{barcode_prefix}"

            adata = ad.read_h5ad("~{h5ad}", backed="r")
            adata = adata[adata.obs.index.isin(barcodes), :]
            adata.write_h5ad(f"tmp/{fname}.h5ad")

            adata = ad.read_h5ad(f"tmp/{fname}.h5ad")
            if ~{if recompute_metrics then "True" else "False"}:
                qc_vars = ["mito"] if "_gex_" in "~{prefix}" else []
                sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=None, log1p=False, inplace=True)

                gene_subset = calculate_highly_expressed_genes(
                    adata.X.data,
                    adata.X.indices,
                    adata.X.indptr,
                    adata.obs.total_counts.values,
                    ~{high_expression_max_frac},
                    adata.X.shape[1],
                )
                print(f"Excluded genes during normalization:\n{adata.var_names[~gene_subset].tolist()}")
                adata.obs["total_counts_wo_highly_expressed"] = np.ravel(adata.X[:, gene_subset].sum(axis=1))
                adata.write_h5ad(f"{fname}.h5ad")
            else:
                shutil.copy(f"tmp/{fname}.h5ad", f"{fname}.h5ad")

            write_bgz(adata.obs, f"{fname}.obs.tsv.gz", index_label="barcode", sep="\t", na_rep="NA")
            write_bgz(adata.var, f"{fname}.var.tsv.gz", index_label="gene_id", sep="\t", na_rep="NA")

            os.remove(f"tmp/{fname}.h5ad")
        __EOF__

        mkdir -p tmp && \
        python3 script.py && \
        touch _SUCCESS
    >>>

    output {
        File success = "_SUCCESS"
        Array[File] out_h5ad = glob("~{prefix}.*.h5ad")
        Array[File] out_obs = glob("~{prefix}.*.obs.tsv.gz")
        Array[File] out_var = glob("~{prefix}.*.var.tsv.gz")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "384 GB"
        disks: "local-disk 512 HDD"
        zones: zones
        preemptible: 2
    }
}

task pseudobulk {
    input {
        File h5ad
        File tss_bed
        String sample_id
        String agg_method
        String? obs_key
        File? obs_annot
        String? obs_annot_key
        String cell_type_annot
        String? all_cell_annot
        Array[String]? extra_conditions
        File? keep_sample_id
        File? keep_obs_key
        Int? min_cells
        Float? min_sample_prop
        Float? min_cell_prop
        Int scale_factor = 1000000
        Boolean normalize_by_peak_length = false
        String prefix = basename(h5ad, ".h5ad") + "." + cell_type_annot
        String docker
        String zones
    }

    Int memory = ceil(5.5 * size(h5ad, 'GB')) + 10

    command <<<
        set -e

        pseudobulk.py \
        --h5ad "~{h5ad}" \
        --bed "~{tss_bed}" \
        --sample-id "~{sample_id}" \
        --agg-method "~{agg_method}" \
        ~{true='--obs-key ' false='' defined(obs_key)}~{obs_key} \
        ~{true='--obs-annot ' false='' defined(obs_annot)}~{obs_annot} \
        ~{true='--obs-annot-key ' false='' defined(obs_annot_key)}~{obs_annot_key} \
        --cell-type-annot "~{cell_type_annot}" \
        ~{true='--all-cell-annot ' false='' defined(all_cell_annot)}~{all_cell_annot} \
        ~{true='--extra-conditions ' false='' defined(extra_conditions)}~{sep=' ' extra_conditions} \
        ~{true='--keep-sample-id ' false='' defined(keep_sample_id)}~{keep_sample_id} \
        ~{true='--keep-obs-key ' false='' defined(keep_obs_key)}~{keep_obs_key} \
        ~{true='--min-cells ' false='' defined(min_cells)}~{min_cells} \
        ~{true='--min-sample-prop ' false='' defined(min_sample_prop)}~{min_sample_prop} \
        ~{true='--min-cell-prop ' false='' defined(min_cell_prop)}~{min_cell_prop} \
        ~{true='--normalize-by-peak-length ' false='' normalize_by_peak_length} \
        --scale-factor ~{scale_factor} \
        --export-cell-counts \
        --write-raw \
        --prefix "~{prefix}"
    >>>

    output {
        Array[File] out_inv_bed = glob("*.inv.bed.gz")
        Array[File] out_raw_bed = glob("*.raw.bed.gz")
        Array[File] out_nonzero_prop = glob("*.nonzero_prop.txt")
        Array[File] out_barcode = glob("*.barcode.txt.gz")
        File out_cell_count = prefix + ".cell_counts.txt"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{memory} GB"
        disks: "local-disk ~{ceil(2 * size(h5ad, 'GB')) + 10} HDD"
        zones: zones
        preemptible: 2
    }
}

task split_ensgene {
    input {
        File gex_features
        String chrom
        Int num_genes
        String docker
        String zones
    }

    command <<<
        set -e

        zcat ~{gex_features} | awk -v chrom="~{chrom}" '
        BEGIN {
            FS = "\t"
        }
        $4 == chrom {
            print $1
        }
        ' > ensgene_ids.txt
        split -l ~{num_genes} -d --verbose ensgene_ids.txt ensgene_list.
    >>>

    output {
        Array[File] out_lists = glob("ensgene_list.*")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "7 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}

task combine_pseudobulk {
    input {
        Array[File] inv_bed
        Array[File] raw_bed
        Array[File] nonzero_prop
        Array[File] barcode
        String agg_method
        String prefix = sub(basename(inv_bed[0], ".~{agg_method}.inv.bed.gz"), "\\.chr[0-9X]+", "")
        String docker
        String zones
    }

    command <<<
        set -e

        zcat ~{sep=" " inv_bed} | awk 'NR == 1 || $2 != "start"' | bgzip -c > ~{prefix}.~{agg_method}.inv.bed.gz
        zcat ~{sep=" " raw_bed} | awk 'NR == 1 || $2 != "start"' | bgzip -c > ~{prefix}.~{agg_method}.raw.bed.gz
        awk 'NR == 1 || FNR > 1' ~{sep=" " nonzero_prop} > ~{prefix}.nonzero_prop.txt
        zcat ~{sep=" " barcode} | bgzip -c > ~{prefix}.barcode.txt.gz

    >>>

    output {
        File out_inv_bed = "~{prefix}.~{agg_method}.inv.bed.gz"
        File out_raw_bed = "~{prefix}.~{agg_method}.raw.bed.gz"
        File out_nonzero_prop = "~{prefix}.nonzero_prop.txt"
        File out_barcode = "~{prefix}.barcode.txt.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "7 GB"
        disks: "local-disk 20 HDD"
        zones: zones
        preemptible: 2
    }
}

task generate_peak_bed {
    input {
        File var
        String docker
        String zones
    }

    command <<<
        set -e

        zcat ~{var} | awk '
        BEGIN {
            OFS = "\t"
        }
        NR > 1 {
            split($1, a, "-")
            print a[1], a[2], a[3], $1
        }
        ' > peak.bed

    >>>

    output {
        File out_bed = "peak.bed"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "7 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}