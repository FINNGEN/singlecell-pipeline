#!/usr/bin/env python3

import anndata as ad
import argparse
import bgzip
import numpy as np
import pandas as pd
import re
import scanpy as sc
import scipy.sparse as sp
from numba import njit


def write_bgz(df, out_file, **kwargs):
    with open(out_file, "wb") as raw:
        with bgzip.BGZipWriter(raw) as fh:
            fh.write(df.to_csv(**kwargs).encode("utf-8"))


def sanitize_filename(filename):
    filename = re.sub(r"[^a-zA-Z0-9_.-]", "_", filename)
    filename = re.sub(r"[_-]+", "_", filename)
    return filename


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


def main(args):
    adata = ad.read_h5ad(args.h5ad)
    adata.raw = None

    if args.type == "gex" and args.features is not None:
        features = pd.read_csv(
            args.features,
            names=["gene_id", "symbol", "type", "chrom", "start", "end", "gene_type"],
            header=None,
            sep="\t",
        )

        features_index = "gene_id" if np.all(adata.var_names.str.startswith("ENSG")) else "symbol"
        features = features.set_index(features_index).drop("type", axis=1)

        adata.var["idx"] = np.arange(len(adata.var.index))
        adata.var = (
            adata.var.merge(features, left_index=True, right_index=True).sort_values("idx").drop(columns=["idx"])
        )

        if features_index == "symbol":
            adata.var = adata.var.reset_index(names=["symbol"]).set_index("gene_id")

        if np.any(adata.var.index.duplicated()):
            print(adata.var)
            raise ValueError("var_names are not unique")

    if args.type == "atac":
        # peak_id = gene_id for atac
        adata.var.index = adata.var.index.rename("gene_id")

    if args.obs_annot is not None:
        if "barcode" in adata.obs.columns:
            adata.obs = adata.obs.drop(columns=["barcode"])
        adata.obs.index = adata.obs.index.rename("barcode")
        adata.obs["idx"] = np.arange(len(adata.obs.index))

        # primary_obs_key (e.g. sample_id) should be the first one
        primary_obs_key = args.obs_key[0]
        if primary_obs_key not in adata.obs.columns:
            print(f"{primary_obs_key} not found: imputing from barcodes")
            adata.obs[primary_obs_key] = adata.obs.index.to_series().str.slice(
                stop=-(args.barcode_length + len(args.barcode_sep))
            )

        # remove _2 suffix
        if args.resequenced_suffix is not None:
            resequenced_key = primary_obs_key + args.resequenced_suffix
            adata.obs[resequenced_key] = adata.obs[primary_obs_key]
            adata.obs[primary_obs_key] = adata.obs[primary_obs_key].str.replace(f"_[0-9]$", "", regex=True)
            adata.obs["is_pool"] = (
                adata.obs[args.obs_pool_id].str.startswith("pool") if args.obs_pool_id is not None else False
            )

            # keep resequened suffix for samples with multiple libraries
            suffix = (
                adata.obs[[primary_obs_key, resequenced_key, "is_pool"]]
                .sort_values([primary_obs_key, resequenced_key, "is_pool"])
                .drop_duplicates()
                .assign(
                    idx=lambda x: x.groupby(primary_obs_key).cumcount().add(1),
                    # f-string doesn't work here
                    resequenced_suffix=lambda x: np.where(
                        x.is_pool,
                        x[resequenced_key].str.replace("^.*(_[0-9])$", "\\1", regex=True),
                        np.where(
                            x.idx > 1,
                            "_" + x.idx.astype(str),
                            "",
                        ),
                    ),
                )
                .set_index(primary_obs_key)
                .drop(columns=["idx", "is_pool"])
            )
            # add suffix for annotating barcodes from samples with multiple libraries
            adata.obs = (
                adata.obs.reset_index()
                .merge(suffix, left_on=resequenced_key, right_on=resequenced_key, how="left")
                .drop(columns=["is_pool"])
                .set_index("barcode")
            )
        else:
            adata.obs["resequenced_suffix"] = ""
        print(adata.obs)

        obs = adata.obs.copy().reset_index()
        for i in range(len(args.obs_annot)):
            obs_annot = pd.read_csv(args.obs_annot[i], sep="\t")
            obs = obs.merge(obs_annot, left_on=args.obs_key[i], right_on=args.obs_annot_key[i], how="left")
        adata.obs = obs.sort_values("idx").drop(columns=["idx"]).set_index("barcode")
        print(adata.obs)

        if args.update_barcodes:
            adata.obs = (
                adata.obs.reset_index()
                .assign(
                    s=lambda x: x[args.obs_sample_id].str.cat(x.resequenced_suffix, sep=""),
                    barcode=lambda x: x.s.str.cat(
                        x.barcode.str.slice(start=-args.barcode_length), sep=args.barcode_sep
                    ),
                )
                .drop(columns=["s", "resequenced_suffix"])
                .set_index("barcode")
            )
            print(adata.obs)

        if np.any(adata.obs.index.duplicated()):
            print(adata.obs.loc[adata.obs.index.duplicated(), :])
            raise ValueError("obs_names are not unique")

    obs_idx = np.array([True] * len(adata.obs.index), dtype=bool)
    if args.vars_to_regress is not None:
        vars_to_regress = list(set(args.vars_to_regress) - {"pct_counts_mito"})
        obs_idx = adata.obs[vars_to_regress].notnull().all(axis=1)

    if args.keep_samples is not None:
        keep = pd.read_csv(args.keep_samples, header=None).iloc[:, 0]
        obs_idx = np.logical_and(obs_idx, adata.obs[args.obs_sample_id].isin(keep))
    if not np.all(obs_idx):
        adata._inplace_subset_obs(obs_idx)

    if args.downsample is not None:
        adata.strings_to_categoricals()
        rng = np.random.default_rng(seed=args.seed)
        rnd_idx = (rng.uniform(0, 1, adata.obs.shape[0]) + adata.obs[args.obs_sample_id].cat.codes) > (
            adata.obs[args.obs_sample_id].cat.codes + (1 - args.downsample)
        )
        adata._inplace_subset_obs(rnd_idx)

    if args.filter_cells:
        sc.pp.filter_cells(adata, min_genes=args.min_genes)
        if args.min_counts is not None:
            sc.pp.filter_cells(adata, min_counts=args.min_counts)
        if args.max_counts is not None:
            sc.pp.filter_cells(adata, max_counts=args.max_counts)
        sc.pp.filter_genes(adata, min_cells=args.min_cells)
    if args.filter_samples:
        keep = adata.obs[args.obs_sample_id].value_counts() > args.min_cells_per_sample
        obs_idx = adata.obs[args.obs_sample_id].isin(keep[keep].index)
        adata._inplace_subset_obs(obs_idx)

    if args.calculate_qc_metrics:
        qc_vars = []
        if args.type == "gex":
            adata.var["mito"] = adata.var.symbol.str.startswith("MT-")
            qc_vars = ["mito"]
        sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=None, log1p=False, inplace=True)

        # Calculate total counts without highly expressed genes
        # https://github.com/scverse/scanpy/blob/1.11.1/src/scanpy/preprocessing/_normalization.py#L59-L288
        gene_subset = calculate_highly_expressed_genes(
            adata.X.data,
            adata.X.indices,
            adata.X.indptr,
            adata.obs.total_counts.values,
            args.high_expression_max_frac,
            adata.X.shape[1],
        )

        print(f"Excluded genes during normalization:\n{adata.var_names[~gene_subset].tolist()}")
        adata.obs["total_counts_wo_highly_expressed"] = np.ravel(adata.X[:, gene_subset].sum(axis=1))

    if args.sanitize_cell_type_annots is not None:
        for col in args.sanitize_cell_type_annots:
            adata.obs[col] = adata.obs[col].apply(sanitize_filename)

    if args.drop_obs_columns is not None:
        adata.obs = adata.obs.drop(columns=args.drop_obs_columns)
    if args.drop_var_columns is not None:
        adata.var = adata.var.drop(columns=args.drop_var_columns)

    # make sure that string columns are categorical
    adata.strings_to_categoricals()
    # Check if adata.X is a sparse matrix and convert to CSC if needed
    if not sp.isspmatrix_csc(adata.X):
        adata.X = sp.csc_matrix(adata.X)
    # make sure that X is uint32
    # https://github.com/bnprks/BPCells/issues/49#issuecomment-1932869295
    adata.X = adata.X.astype(np.uint32)
    adata.write(f"{args.prefix}.h5ad")
    write_bgz(adata.obs, f"{args.prefix}.obs.tsv.gz", index_label="barcode", sep="\t", na_rep="NA")
    write_bgz(adata.var, f"{args.prefix}.var.tsv.gz", index_label="gene_id", sep="\t", na_rep="NA")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", type=str, required=True)
    parser.add_argument("--type", type=str, choices=["gex", "atac"], default="gex")
    parser.add_argument("--h5ad", type=str, required=True)
    parser.add_argument("--features", type=str)
    parser.add_argument("--obs-annot", type=str, nargs="+")
    parser.add_argument("--obs-key", type=str, nargs="+")
    parser.add_argument("--obs-annot-key", type=str, nargs="+")
    parser.add_argument("--obs-sample-id", type=str)
    parser.add_argument("--obs-pool-id", type=str)
    parser.add_argument("--drop-obs-columns", type=str, nargs="+")
    parser.add_argument("--drop-var-columns", type=str, nargs="+")
    parser.add_argument("--update-barcodes", action="store_true")
    parser.add_argument("--barcode-length", type=int, default=16)
    parser.add_argument("--barcode-sep", type=str, default="-")
    parser.add_argument("--sanitize-cell-type-annots", type=str, nargs="+")
    parser.add_argument("--keep-samples", type=str, help="Keep only samples")
    parser.add_argument("--calculate-qc-metrics", action="store_true")
    parser.add_argument("--high-expression-max-frac", type=float, default=0.05)
    parser.add_argument("--resequenced-suffix", type=str, help='Suffix for resequenced samples (e.g. "_2")')
    parser.add_argument("--filter-cells", action="store_true")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-counts", type=int, default=200)
    parser.add_argument("--max-counts", type=int)
    parser.add_argument("--min-cells", type=int, default=10)
    parser.add_argument("--filter-samples", action="store_true")
    parser.add_argument("--min-cells-per-sample", type=int, default=1000)
    parser.add_argument("--vars-to-regress", type=str, nargs="+")
    parser.add_argument("--downsample", type=float, help="Fraction of downsampling")
    parser.add_argument("--seed", type=int, default=123456, help="Seed of random generator for downsampling")
    args = parser.parse_args()

    if (args.obs_annot is not None) and ((args.obs_key is None) or (args.obs_annot_key is None)):
        raise ValueError("--obs-annot requires --obs-key and --obs-annot-key")
    if args.update_barcodes and (args.obs_sample_id is None):
        raise ValueError("--update-barcodes requires --obs-sample-id")
    if (args.downsample is not None) and ((args.downsample < 0) or (args.downsample > 1)):
        raise ValueError("--downsample value should be between 0 and 1")
    if len(args.obs_annot) != len(args.obs_key):
        raise ValueError("Number of --obs-annot should match number of --obs-key")
    if len(args.obs_annot) != len(args.obs_annot_key):
        raise ValueError("Number of --obs-annot should match number of --obs-annot-key")

    main(args)
