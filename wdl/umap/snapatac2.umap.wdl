version 1.0

import "snapatac2.tasks.wdl" as tasks

workflow snapatac2_umap {
    input {
        String obs
        String h5ad_pattern
        String out_prefix
        File cell_type_list
        String all_cell_annot
        Int n_features = 50000
        Int n_components = 50
        Boolean write_peaks_h5ad = false
        Boolean debug = false
        String docker_snapatac2
        String docker_suite
        String zones
    }

    call tasks.prepare_h5ad_list {
        input:
            obs = obs,
            h5ad_pattern = h5ad_pattern,
            docker = docker_suite,
            zones = zones
    }

    scatter (h5ad in prepare_h5ad_list.out_h5ad) {
        call tasks.convert_to_csr {
            input:
                h5ad = h5ad,
                docker = docker_snapatac2,
                zones = zones
        }
    }

    Array[String] cell_types = if debug then [read_lines(cell_type_list)[0]] else read_lines(cell_type_list)

    scatter (cell_type in cell_types) {
        call run_umap {
            input:
                obs = obs,
                h5ad = convert_to_csr.out_h5ad,
                cell_type = cell_type,
                all_cell_annot = all_cell_annot,
                prefix = "~{out_prefix}.~{cell_type}",
                n_features = n_features,
                n_components = n_components,
                write_peaks_h5ad = write_peaks_h5ad,
                docker = docker_snapatac2,
                zones = zones
        }
    }

    output {
        Array[File] out_lsi = run_umap.out_lsi
        Array[File] out_h5ad = run_umap.out_h5ad
        Array[File] out_h5ads = run_umap.out_h5ads
        Array[File?] out_peaks_h5ad = run_umap.out_peaks_h5ad
    }
}

task run_umap {
    input {
        String obs
        Array[File] h5ad
        String cell_type
        String all_cell_annot
        String prefix
        Int n_features = 50000
        Int n_components = 50
        Boolean write_peaks_h5ad = false
        String docker
        String zones
    }

    command <<<
        set -e

        export n_cpu=$(grep -c ^processor /proc/cpuinfo)
        ulimit -n 16384

        cat << "__EOF__" > script.py
        import anndata as ad
        import numpy as np
        import os
        import os.path
        import pandas as pd
        import scanpy as sc
        import snapatac2 as snap

        n_jobs = int(os.environ['n_cpu'])
        cell_type = "~{cell_type}"

        # Parse cell_type: "predicted.celltype.l1.CD4_T" -> col="predicted.celltype.l1", val="CD4_T"
        cell_type_col = ".".join(cell_type.rsplit(".", 1)[:-1])
        cell_type_val = cell_type.rsplit(".", 1)[-1]

        # Build ffid_barcode -> barcode mapping from obs
        obs = pd.read_csv("~{obs}", sep="\t", index_col=0)
        suffix = obs.index.to_series().str.split('-', expand=True).iloc[:, 1]
        mask = obs['pool_name'].str.startswith('pool')
        obs['ffid_barcode'] = obs['ffid_2'] + '-' + suffix
        obs.loc[mask, 'ffid_barcode'] = (
            obs.loc[mask, 'pool_name'] + '_' + obs.loc[mask, 'ffid'] + '-' + suffix[mask]
        )
        obs = obs.reset_index(names="barcode").set_index('ffid_barcode')

        # Create AnnDataSet
        h5ad_paths = "~{sep=',' h5ad}".split(",")
        adatas = [(os.path.basename(x).split(".")[0], x) for x in h5ad_paths]
        adataset = snap.AnnDataSet(adatas=adatas, filename="~{prefix}.h5ads", add_key="key")

        # Subset to cell type using obs file.
        # Two supported configurations:
        #   (A) cell_type != all_cell_annot: obs is multi-cell-type with a
        #       cell_type_col. We intersect adataset with obs (defensively
        #       dropping adataset cells absent from obs, e.g. when obs was
        #       pre-filtered) and then keep cells labeled cell_type_val.
        #   (B) cell_type == all_cell_annot: treat obs as the definitive
        #       whitelist. If obs is full (every adataset cell present),
        #       this is a no-op; if obs is pre-filtered to a subset, we
        #       intersect adataset with obs. Either way, no cell_type_col
        #       filter is applied here.
        adataset_names = list(adataset.obs_names)
        obs_set = set(obs.index)
        if "~{cell_type}" != "~{all_cell_annot}":
            in_obs_idx = [i for i, b in enumerate(adataset_names) if b in obs_set]
            n_dropped = len(adataset_names) - len(in_obs_idx)
            if n_dropped > 0:
                print(f"Pre-filtering adataset: dropping {n_dropped} barcodes not in obs")
            in_obs_names = [adataset_names[i] for i in in_obs_idx]
            obs_ct = obs.loc[in_obs_names, cell_type_col].str.replace(" ", "_")
            keep_idx = [in_obs_idx[i] for i in np.where(obs_ct.to_numpy() == cell_type_val)[0]]
            print(f"Cell type {cell_type}: {len(keep_idx)} / {len(adataset_names)} cells")
            adataset_sub, reorder = adataset.subset(obs_indices=keep_idx, out="~{prefix}.subset.h5ads")
        else:
            in_obs_idx = [i for i, b in enumerate(adataset_names) if b in obs_set]
            n_dropped = len(adataset_names) - len(in_obs_idx)
            if n_dropped > 0:
                print(f"Pre-filtering adataset to {len(in_obs_idx)} cells present in obs (dropped {n_dropped})")
                adataset_sub, reorder = adataset.subset(obs_indices=in_obs_idx, out="~{prefix}.subset.h5ads")
            else:
                print(f"All {len(adataset_names)} adataset cells present in obs; no pre-filter needed")
                adataset_sub = adataset

        # Compute LSI + UMAP
        snap.pp.select_features(adataset_sub, n_features=~{n_features}, max_iter=1, n_jobs=n_jobs)
        snap.tl.spectral(adataset_sub, n_comps=~{n_components}, distance_metric="cosine", feature_weights=None)
        snap.tl.umap(adataset_sub, use_rep="X_spectral", random_state=None)

        # Export LSI components with remapped barcodes
        spectral = adataset_sub.obsm["X_spectral"]
        ffid_barcodes = list(adataset_sub.obs_names)
        barcodes = obs.loc[ffid_barcodes, "barcode"].values

        cols = [f"ATAC_LSI_{i+1}" for i in range(spectral.shape[1])]
        df_lsi = pd.DataFrame(spectral, index=barcodes, columns=cols)
        df_lsi.index.name = "barcode"
        df_lsi.to_csv("~{prefix}.lsi.tsv", index_label="barcode", sep="\t")
        print(f"Exported {df_lsi.shape[0]} cells x {df_lsi.shape[1]} LSI components")

        # Per-cell total fragments across all peaks, for downstream:
        # - accessibility-rate normalization in KNN smoothing
        # - library-size covariate in raw-count models (peak-gene hurdle
        #   negbin, etc.)
        # Prefer an existing obs column when available (free). Otherwise
        # compute from the peak matrix in chunks, to avoid materializing
        # the full sparse matrix for a sum operation at PBMC scale.
        obs_out = obs.loc[ffid_barcodes, :].set_index('barcode').copy()
        _tf_col = next(
            (c for c in ("total_fragments", "n_fragment", "n_fragments", "frag_count") if c in obs_out.columns),
            None,
        )
        if _tf_col is not None:
            print(f"Using existing obs column '{_tf_col}' for total_fragments")
            total_fragments = obs_out[_tf_col].to_numpy().astype(np.int64)
        else:
            print("Computing total_fragments from peak matrix (chunked)...")
            n_cells = adataset_sub.n_obs
            total_fragments = np.zeros(n_cells, dtype=np.int64)
            chunk_size = 50000
            for start in range(0, n_cells, chunk_size):
                end = min(start + chunk_size, n_cells)
                X_chunk = adataset_sub.X[start:end, :]
                total_fragments[start:end] = np.asarray(X_chunk.sum(axis=1)).flatten()
        obs_out['total_fragments'] = total_fragments

        # Export embedding-only AnnData (small; obs + obsm, no X)
        adata = ad.AnnData(
            obs=obs_out,
            obsm={
                "X_spectral": adataset_sub.obsm["X_spectral"],
                "X_umap": adataset_sub.obsm["X_umap"],
            },
        )

        if len(adata.obs_names) != len(adataset_sub.obs_names):
            raise ValueError("Mismatch in number of obs names between adataset and adata")

        adata.strings_to_categoricals()
        adata.write_h5ad("~{prefix}.umap.h5ad")
        print(f"Exported {adata.shape[0]} cells to umap.h5ad")

        # Optional: standard h5ad with the full peak matrix, for downstream
        # chromVAR (full-matrix use) and KNN smoothing (column-wise access via
        # gcs_anndata). Matrix is stored as CSC so gcs_anndata.get_columns()
        # can efficiently fetch a subset of peaks without downloading the
        # full file. Compression is intentionally disabled: GCS disk is cheap
        # relative to the high-memory VMs needed for materialization, and
        # uncompressed HDF5 offsets are random-access friendly for
        # partial-read clients.
        if ~{if write_peaks_h5ad then "True" else "False"}:
            from scipy.sparse import csc_matrix, vstack, issparse
            print("Materializing peak matrix for peaks.h5ad export (chunked CSR, then CSC)...")

            # Chunk the read to avoid holding two full-matrix copies at once
            # during the eventual CSR->CSC transpose. Each chunk is held as
            # CSR (snapatac2's native layout), then vstack'd, then tocsc().
            n_cells = adataset_sub.n_obs
            chunk_size = 100000
            chunks = []
            for start in range(0, n_cells, chunk_size):
                end = min(start + chunk_size, n_cells)
                X_chunk = adataset_sub.X[start:end, :]
                if not issparse(X_chunk):
                    X_chunk = csc_matrix(X_chunk)
                chunks.append(X_chunk.tocsr())
                print(f"  read chunk {start}-{end} ({end-start} cells)")
            X_peaks_csr = vstack(chunks, format="csr") if len(chunks) > 1 else chunks[0]
            del chunks

            # Final transpose to CSC for gcs_anndata column-access efficiency.
            # Peak memory: roughly 2x the sparse size during tocsc().
            print(f"Converting CSR -> CSC: nnz = {X_peaks_csr.nnz}")
            X_peaks = X_peaks_csr.tocsc()
            del X_peaks_csr
            print(f"Peak matrix: {X_peaks.shape[0]} cells x {X_peaks.shape[1]} peaks, nnz = {X_peaks.nnz}")

            var_df = pd.DataFrame(index=pd.Index(adataset_sub.var_names, name=None))

            adata_peaks = ad.AnnData(
                X=X_peaks,
                obs=obs_out,
                var=var_df,
                obsm={
                    "X_spectral": adataset_sub.obsm["X_spectral"],
                    "X_umap": adataset_sub.obsm["X_umap"],
                },
            )
            adata_peaks.strings_to_categoricals()
            # compression=None explicitly: the anndata default is already no
            # compression, but being explicit makes the intent obvious and
            # aligns with the design choice (disk is cheap on GCS relative
            # to the high-memory VMs needed here; uncompressed HDF5 offsets
            # are also random-access friendly for gcs_anndata clients).
            adata_peaks.write_h5ad("~{prefix}.peaks.h5ad", compression=None)
            print(f"Exported peaks.h5ad: {adata_peaks.shape[0]} cells x {adata_peaks.shape[1]} peaks (CSC, uncompressed)")
        __EOF__

        python3 script.py && \
        bgzip "~{prefix}.lsi.tsv" && \
        touch _SUCCESS

    >>>

    output {
        File out_success = "_SUCCESS"
        File out_lsi = "~{prefix}.lsi.tsv.gz"
        File out_h5ad = "~{prefix}.umap.h5ad"
        File out_h5ads = "~{prefix}.h5ads"
        File? out_peaks_h5ad = "~{prefix}.peaks.h5ad"
    }

    runtime {
        docker: docker
        cpu: 96
        memory: "624 GB"
        disks: "local-disk 1200 HDD"
        zones: zones
        preemptible: 2
    }
}
