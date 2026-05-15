version 1.0

# iLISI and cLISI on integrated embeddings (Korsunsky et al., 2019, Nat. Methods).
#
# iLISI quantifies batch mixing, cLISI cell-type preservation, as the
# effective number of neighbors per cell in a local neighborhood (perplexity).
# Higher iLISI = better batch mixing (range 1..n_batches). Lower cLISI
# (closer to 1) = better cell-type preservation (range 1..n_cell_types).
#
# Runs on scanpy_umap or snapatac2_umap output h5ads, both of which store a
# low-dimensional embedding in obsm (X_pca for RNA, X_spectral for ATAC) and
# batch/cell-type columns in obs. X need not be present.
#
# Three tasks: (1) extract_embedding writes TSVs of the embedding and obs
# once; (2) run_lisi computes LISI in R (lisi package), scattered over
# seeds for subsample-robustness; (3) combine_seeds aggregates per-seed
# summaries into median/IQR across seeds. For subsample = 0 (all cells)
# seeds are no-ops and a single-element seeds array suffices.

workflow compute_lisi {
    input {
        File h5ad
        String embedding_key           # "X_pca" for RNA, "X_spectral" for ATAC
        Array[String] batch_cols       # obs columns for iLISI, e.g. ["pool_name", "batch"]
        String celltype_col = ""       # obs column for cLISI; "" to skip
        String label = ""              # written into summary (e.g. "RNA" / "ATAC")
        Float perplexity = 30.0
        Int subsample = 0              # 0 = all cells
        Array[Int] seeds = [42]        # one scatter shard per seed; use >1 with subsampling
        String out_prefix = basename(h5ad, ".h5ad") + ".lisi"
        String docker_suite            # Python + anndata + pandas; also R + lisi
        String zones
    }

    call extract_embedding {
        input:
            h5ad = h5ad,
            embedding_key = embedding_key,
            prefix = basename(h5ad, ".h5ad"),
            docker = docker_suite,
            zones = zones
    }

    scatter (seed in seeds) {
        call run_lisi {
            input:
                embedding = extract_embedding.out_embedding,
                metadata = extract_embedding.out_obs,
                batch_cols = batch_cols,
                celltype_col = celltype_col,
                label = label,
                perplexity = perplexity,
                subsample = subsample,
                seed = seed,
                prefix = "~{out_prefix}.seed~{seed}",
                docker = docker_suite,
                zones = zones
        }
    }

    call combine_seeds {
        input:
            summaries = run_lisi.out_summary,
            prefix = out_prefix,
            docker = docker_suite,
            zones = zones
    }

    output {
        File out_embedding = extract_embedding.out_embedding
        File out_obs = extract_embedding.out_obs
        Array[File] out_per_cell = run_lisi.out_per_cell
        Array[File] out_summary = run_lisi.out_summary
        File out_all_seeds = combine_seeds.out_all_seeds
        File out_combined = combine_seeds.out_combined
    }
}


task extract_embedding {
    input {
        File h5ad
        String embedding_key
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        python3 << "__EOF__"
        import anndata as ad
        import pandas as pd

        adata = ad.read_h5ad("~{h5ad}")
        embedding_key = "~{embedding_key}"

        if embedding_key not in adata.obsm:
            raise KeyError(
                f"{embedding_key} not in obsm. Available: {list(adata.obsm.keys())}"
            )

        X = adata.obsm[embedding_key]
        n_dims = X.shape[1]
        cols = [f"dim_{i+1}" for i in range(n_dims)]
        df_emb = pd.DataFrame(X, index=adata.obs_names, columns=cols)
        df_emb.to_csv(
            "~{prefix}.embedding.tsv.gz",
            sep="\t", index_label="barcode", compression="gzip",
        )
        print(f"Wrote embedding: {X.shape[0]} cells x {n_dims} dims")

        # Always use adata.obs_names as the authoritative barcode. Drop any
        # pre-existing "barcode" column (e.g. the snapatac2 atlas carries one)
        # to avoid a duplicate-column collision on insert.
        obs = adata.obs.copy()
        if "barcode" in obs.columns:
            obs = obs.drop(columns=["barcode"])
        obs.insert(0, "barcode", adata.obs_names)
        obs = obs.reset_index(drop=True)
        obs.to_csv(
            "~{prefix}.obs.tsv.gz",
            sep="\t", index=False, compression="gzip",
        )
        print(f"Wrote obs: {obs.shape[0]} rows x {obs.shape[1]} cols")
        __EOF__
    >>>

    output {
        File out_embedding = "~{prefix}.embedding.tsv.gz"
        File out_obs = "~{prefix}.obs.tsv.gz"
    }

    runtime {
        docker: docker
        cpu: 2
        memory: "~{ceil(4 * size(h5ad, 'GB')) + 16} GB"
        disks: "local-disk ~{ceil(3 * size(h5ad, 'GB')) + 30} HDD"
        zones: zones
        preemptible: 2
    }
}


task run_lisi {
    input {
        File embedding
        File metadata
        Array[String] batch_cols
        String celltype_col
        String label
        Float perplexity
        Int subsample
        Int seed
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        cat << "__EOF__" > compute_lisi.R
        suppressPackageStartupMessages({
            library(data.table)
            library(lisi)
        })

        batch_cols   <- strsplit("~{sep=',' batch_cols}", ",", fixed = TRUE)[[1]]
        celltype_col <- "~{celltype_col}"
        perplexity   <- ~{perplexity}
        subsample    <- ~{subsample}
        seed         <- ~{seed}
        label        <- "~{label}"
        prefix       <- "~{prefix}"

        message("Loading embedding: ~{embedding}")
        emb <- fread("~{embedding}")
        stopifnot("barcode" %in% names(emb))
        non_bc <- setdiff(names(emb), "barcode")
        is_num <- vapply(non_bc, function(c) is.numeric(emb[[c]]), logical(1))
        dim_cols <- non_bc[is_num]
        if (length(dim_cols) < 2) {
            stop("Fewer than 2 numeric embedding columns detected.")
        }
        message(sprintf("  %d cells x %d dims", nrow(emb), length(dim_cols)))

        message("Loading metadata: ~{metadata}")
        meta <- fread("~{metadata}")
        stopifnot("barcode" %in% names(meta))

        use_clisi <- nzchar(celltype_col)
        stopifnot(all(batch_cols %in% names(meta)))
        if (use_clisi) {
            stopifnot(celltype_col %in% names(meta))
        }
        label_cols <- if (use_clisi) c(batch_cols, celltype_col) else batch_cols

        setkey(emb, barcode)
        setkey(meta, barcode)
        common <- intersect(emb[["barcode"]], meta[["barcode"]])
        if (length(common) < 100) {
            stop(sprintf("Only %d barcodes common between embedding and metadata",
                         length(common)))
        }
        message(sprintf("  %d barcodes common", length(common)))
        emb  <- emb[.(common),  on = "barcode"]
        meta <- meta[.(common), on = "barcode"]

        if (subsample > 0L && subsample < nrow(emb)) {
            set.seed(seed)
            idx <- sort(sample.int(nrow(emb), subsample))
            message(sprintf("Subsampling %d cells (of %d)", subsample, nrow(emb)))
            emb  <- emb[idx, ]
            meta <- meta[idx, ]
        }

        x_mat <- as.matrix(emb[, ..dim_cols])
        rownames(x_mat) <- emb[["barcode"]]
        meta_df <- as.data.frame(meta)
        rownames(meta_df) <- meta_df[["barcode"]]

        ok <- Reduce(`&`, lapply(label_cols, function(c) !is.na(meta_df[[c]])))
        if (!all(ok)) {
            message(sprintf("Dropping %d cells with NA in a label column", sum(!ok)))
            x_mat   <- x_mat[ok, , drop = FALSE]
            meta_df <- meta_df[ok, , drop = FALSE]
        }

        message(sprintf("Computing LISI (perplexity = %g) on %d cells for: %s",
                        perplexity, nrow(x_mat), paste(label_cols, collapse = ", ")))
        lisi_scores <- lisi::compute_lisi(
            X              = x_mat,
            meta_data      = meta_df,
            label_colnames = label_cols,
            perplexity     = perplexity
        )

        per_cell <- data.table(
            barcode = rownames(x_mat),
            as.data.table(lisi_scores)
        )
        fwrite(per_cell, sprintf("%s.per_cell.tsv.gz", prefix),
               sep = "\t", quote = FALSE)

        qs <- c(0.25, 0.5, 0.75)
        summary_rows <- lapply(label_cols, function(col) {
            vals <- lisi_scores[[col]]
            kind <- if (use_clisi && col == celltype_col) "cLISI" else "iLISI"
            q <- quantile(vals, qs, na.rm = TRUE)
            data.frame(
                label    = label,
                metric   = kind,
                column   = col,
                n_levels = length(unique(meta_df[[col]])),
                q25      = q[[1]],
                median   = q[[2]],
                q75      = q[[3]],
                mean     = mean(vals, na.rm = TRUE)
            )
        })
        summary_df <- do.call(rbind, summary_rows)
        fwrite(summary_df, sprintf("%s.summary.tsv", prefix),
               sep = "\t", quote = FALSE)

        message("Done.")
        __EOF__

        Rscript compute_lisi.R && touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_per_cell = "~{prefix}.per_cell.tsv.gz"
        File out_summary = "~{prefix}.summary.tsv"
    }

    runtime {
        docker: docker
        cpu: 4
        memory: "~{ceil(20 * size(embedding, 'GB')) + 32} GB"
        disks: "local-disk ~{ceil(5 * size(embedding, 'GB')) + 30} HDD"
        zones: zones
        preemptible: 2
    }
}


task combine_seeds {
    input {
        Array[File] summaries
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        python3 << "__EOF__"
        import os
        import re
        import pandas as pd

        paths = "~{sep=',' summaries}".split(",")
        seed_re = re.compile(r"\.seed(-?\d+)\.summary\.tsv$")

        frames = []
        for p in paths:
            df = pd.read_csv(p, sep="\t")
            m = seed_re.search(os.path.basename(p))
            df["seed"] = int(m.group(1)) if m else -1
            frames.append(df)
        all_df = pd.concat(frames, ignore_index=True)
        all_df.to_csv("~{prefix}.all_seeds.tsv", sep="\t", index=False)
        print(f"Wrote per-seed summary: {all_df.shape[0]} rows from {len(paths)} seeds")

        # Aggregate per-seed medians/means across seeds.
        group_keys = ["label", "metric", "column", "n_levels"]
        agg = (
            all_df
            .groupby(group_keys, dropna=False)
            .agg(
                n_seeds=("seed", "count"),
                seed_median_q25=("median", lambda x: x.quantile(0.25)),
                seed_median_median=("median", "median"),
                seed_median_q75=("median", lambda x: x.quantile(0.75)),
                seed_mean_q25=("mean", lambda x: x.quantile(0.25)),
                seed_mean_median=("mean", "median"),
                seed_mean_q75=("mean", lambda x: x.quantile(0.75)),
            )
            .reset_index()
        )
        agg.to_csv("~{prefix}.combined.tsv", sep="\t", index=False)
        print(f"Wrote combined summary: {agg.shape[0]} rows")
        __EOF__
    >>>

    output {
        File out_all_seeds = "~{prefix}.all_seeds.tsv"
        File out_combined = "~{prefix}.combined.tsv"
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
