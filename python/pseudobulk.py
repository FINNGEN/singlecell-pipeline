#!/usr/bin/env python3

import anndata as ad
import argparse
import itertools
import numpy as np
import pandas as pd
import qtl.io
import re
import scanpy as sc
from scipy import sparse
from scipy.stats import rankdata, norm


def cpm(counts, scale_factor=1e6):
    return scale_factor * counts / np.sum(counts, axis=0)


def inv_normlize(x):
    return norm.ppf((rankdata(x) - 0.5) / np.sum(~np.isnan(x)))


def agg_by(adata: ad.AnnData, cols: list, method: str) -> ad.AnnData:
    """
    Aggregate `.X` entries for each unique combination of values in the provided columns,
    according to the specified method (mean or sum).
    Modified from this post: https://www.kaggle.com/code/awater1223/op2-02-pseudobulk-v2

    Parameters:
    adata: AnnData object to process.
    cols: List of column names to combine for aggregation.
    method: Aggregation method, either 'sum' or 'mean'.
    """
    assert method in ["sum", "mean"], "Method must be 'sum' or 'mean'"

    # Validate columns and create a composite column
    for col in cols:
        assert isinstance(adata.obs[col].dtype, pd.CategoricalDtype), f"{col} is not a categorical dtype"

    adata.obs["composite"] = adata.obs[cols].astype(str).agg("-".join, axis=1)

    # Aggregate based on the composite column
    composite_cat = adata.obs["composite"].astype("category")
    indicator = sparse.coo_matrix(
        (np.broadcast_to(True, adata.n_obs), (composite_cat.cat.codes, np.arange(adata.n_obs))),
        shape=(len(composite_cat.cat.categories), adata.n_obs),
    )

    agg_matrix = indicator @ adata.X
    if method == "mean":
        agg_matrix = agg_matrix / indicator.sum(axis=1)

    if sparse.issparse(agg_matrix):
        agg_matrix = agg_matrix.toarray()

    agg_adata = ad.AnnData(
        agg_matrix,
        var=adata.var,
        obs=pd.DataFrame(index=composite_cat.cat.categories),
    )

    # Handle `.obs` attributes and retain original columns
    joining_df = adata.obs[cols + ["composite"]].drop_duplicates().set_index("composite")
    agg_adata.obs = agg_adata.obs.join(joining_df)
    agg_adata.obs.index.name = "composite"

    # Reset the index and change index type to string for consistency
    agg_adata.obs = agg_adata.obs.reset_index()
    agg_adata.obs.index = agg_adata.obs.index.astype("str")

    return agg_adata


def sanitize_filename(filename):
    filename = re.sub(r"[^a-zA-Z0-9_.-]", "_", filename)
    filename = re.sub(r"[_-]+", "_", filename)
    return filename


def indexer(df, cols, combo, all_cell_annot=None):
    if all_cell_annot is not None:
        if len(cols) == 1:
            return np.array([True] * len(df))
        cols = cols[1:]
        combo = combo[1:]

    condition = True
    for col, val in zip(cols, combo):
        condition = np.logical_and(condition, (df[col] == val).values)
    return condition


def main(args):
    if args.raw_h5ad is not None:
        adata = ad.read_h5ad(args.h5ad, backed="r")
        adata_raw = ad.read_h5ad(args.raw_h5ad)
        if not (
            np.all(adata.obs.index.isin(adata_raw.obs.index)) and np.all(adata.var.index.isin(adata_raw.var.index))
        ):
            raise ValueError("adata and adata_raw must have the same shape")
        obs_idx = adata_raw.obs.index.get_indexer(adata.obs.index)
        var_idx = adata_raw.var.index.get_indexer(adata.var.index)
        nonzero_X = adata_raw.X[np.ix_(obs_idx, var_idx)] > 0
        del adata_raw
    else:
        nonzero_X = None

    adata = ad.read_h5ad(args.h5ad)

    # annotate cell types
    if args.obs_annot is not None:
        obs_annot = pd.read_csv(args.obs_annot, sep="\t")
        obs_annot_cols = [
            x
            for x in [args.obs_annot_key, args.cell_type_annot, args.all_cell_annot] + args.extra_conditions
            if x in obs_annot.columns
        ]
        obs_annot = obs_annot[obs_annot_cols]

        adata.obs.index = adata.obs.index.rename("barcode")
        adata.obs["idx"] = np.arange(len(adata.obs.index))

        adata.obs = (
            adata.obs.reset_index()
            .merge(obs_annot, left_on=args.obs_key, right_on=args.obs_annot_key, how="left")
            .sort_values("idx")
            .drop(columns=["idx"])
            .set_index("barcode")
        )
        print(adata.obs)

    # check all_cell_annot
    cell_types = adata.obs[args.cell_type_annot].unique().tolist()
    if args.all_cell_annot in cell_types:
        raise ValueError(f"--all-cell-label should not exist in annotated cell types")

    # filter to valid sample_id
    obs_idx = ~pd.isnull(adata.obs[args.sample_id])
    if args.keep_sample_id is not None:
        obs_idx = np.logical_and(
            obs_idx, adata.obs[args.sample_id].isin(pd.read_csv(args.keep_sample_id, header=None)[0])
        )
    if args.keep_obs_key is not None:
        obs_idx = np.logical_and(obs_idx, adata.obs[args.obs_key].isin(pd.read_csv(args.keep_obs_key, header=None)[0]))
    adata._inplace_subset_obs(obs_idx)
    assert len(adata.obs.index) > 0, f"No valid samples found after --keep-sample-id and --keep-obs-key filtering"

    if args.keep_feature_id is not None:
        var_idx = adata.var.index.isin(pd.read_csv(args.keep_feature_id, header=None)[0])
        adata._inplace_subset_var(var_idx)
    assert len(adata.var.index) > 0, f"No valid features found after --keep-feature-id filtering"

    # get all combinations of conditions
    conditions = [args.cell_type_annot] + args.extra_conditions
    condition_values = [adata.obs[x].unique().tolist() for x in conditions]
    combinations = list(
        [combo for combo in itertools.product(*condition_values) if np.sum(indexer(adata.obs, conditions, combo)) > 0]
    )
    if args.all_cell_annot is not None:
        condition_values[0] = [args.all_cell_annot]
        combinations += list(
            [
                combo
                for combo in itertools.product(*condition_values)
                if np.sum(indexer(adata.obs, conditions, combo, all_cell_annot=args.all_cell_annot)) > 0
            ]
        )
    print(combinations)

    if nonzero_X is None:
        nonzero_X = adata.X > 0

    # compute nonzero proportions before normalization
    nonzero_prop = {
        combo: np.array(
            np.mean(
                nonzero_X[
                    indexer(
                        adata.obs,
                        conditions,
                        combo,
                        all_cell_annot=(None if combo[0] != args.all_cell_annot else args.all_cell_annot),
                    ),
                    :,
                ],
                axis=0,
            )
        ).ravel()
        for combo in combinations
    }

    # CPM normalization per cell
    if args.agg_method == "mean":
        sc.pp.normalize_total(adata, target_sum=args.scale_factor)

    # psedubulk all cells
    agg_method = args.agg_method.replace("_sct_pearson", "")
    if args.all_cell_annot is not None:
        if args.all_cell_annot not in adata.obs.columns:
            adata.obs[args.all_cell_annot] = args.all_cell_annot
            adata.strings_to_categoricals()
        padata_all = agg_by(
            adata, cols=[args.sample_id, args.all_cell_annot] + args.extra_conditions, method=agg_method
        )

    # filter to cells with valid annot - this happens after all pseudobulk
    nonmissing_obs = adata.obs.index.isin(adata.obs.dropna(subset=conditions).index)
    adata._inplace_subset_obs(nonmissing_obs)
    padata = agg_by(adata, cols=[args.sample_id] + conditions, method=agg_method)

    if args.bed is not None:
        bed_template_df = pd.read_csv(args.bed, names=["chr", "start", "end", "gene_id"], header=None, sep="\t")
        bed_template_df.index = bed_template_df.gene_id.rename(None)
    elif args.gtf is not None:
        bed_template_df = qtl.io.gtf_to_tss_bed(args.gtf, feature="transcript")

    # compute cell counts per condition
    count = (
        adata.obs.groupby([args.sample_id] + conditions, observed=False)
        .count()
        .iloc[:, 0]
        .rename("count")
        .reset_index()
        .pivot(index=args.sample_id, columns=conditions, values="count")
    )
    if args.all_cell_annot is not None:
        count_all = (
            adata.obs.groupby([args.sample_id, args.all_cell_annot] + args.extra_conditions, observed=False)
            .count()
            .iloc[:, 0]
            .rename("count")
            .reset_index()
            .pivot(index=args.sample_id, columns=[args.all_cell_annot] + args.extra_conditions, values="count")
        )
        count = pd.concat([count, count_all], axis=1)

    if args.export_cell_counts:
        frac = count.div(count.sum(axis=1), axis=0)
        # rename columns for export
        frac.columns = [sanitize_filename(".".join(map(str, x))) for x in count.columns.to_flat_index()]
        frac = pd.DataFrame(
            np.apply_along_axis(inv_normlize, axis=0, arr=frac.values), index=frac.index, columns=frac.columns
        )
        frac.to_csv(f"{args.prefix}.cell_counts.txt", index=True, sep="\t")

    if args.normalize_by_peak_length:
        split_df = pd.DataFrame(adata.var_names.str.split("-", expand=True).tolist(), columns=["chr", "start", "end"])
        split_df[["start", "end"]] = split_df[["start", "end"]].astype(int)
        peak_length_kb = ((split_df["end"] - split_df["start"] + 1.0) / 1000).values

    for combo in combinations:
        cell_type = combo[0]
        barcodes = adata.obs.index[
            indexer(
                adata.obs,
                conditions,
                combo,
                all_cell_annot=None if cell_type != args.all_cell_annot else args.all_cell_annot,
            )
        ].to_series()
        if cell_type == args.all_cell_annot:
            subset = padata_all[indexer(padata_all.obs, conditions, combo, all_cell_annot=args.all_cell_annot), :]
        else:
            subset = padata[indexer(padata.obs, conditions, combo), :]

        fname = sanitize_filename(f"{args.prefix}.{'.'.join(combo)}")
        # export barcodes
        barcodes.to_csv(f"{fname}.barcode.txt.gz", index=False, header=False)
        # export nonzero prop
        adata.var["nonzero_prop"] = nonzero_prop[combo]
        adata.var["pseudo_nonzero_prop"] = np.array(np.mean(subset.X > 0, axis=0)).ravel()
        cols = [x for x in ["symbol", "gene_type", "nonzero_prop", "pseudo_nonzero_prop"] if x in adata.var.columns]
        adata.var[cols].to_csv(f"{fname}.nonzero_prop.txt", index_label="gene_id", sep="\t")

        X = subset.X.T
        # keep samples with >10 cells per sample-cell type pairs
        keep_samples = subset.obs[args.sample_id].isin(
            count.index[count[combo if len(combo) > 1 else combo[0]] > args.min_cells]
        )
        # keep features expressed in >50% samples & >10% cells
        keep_features = np.logical_and(
            np.asarray((np.mean(X != 0, axis=1) > args.min_sample_prop)).squeeze(),
            nonzero_prop[combo] > args.min_cell_prop,
        )

        print(combo)
        print(keep_features.shape)
        print(X[np.ix_(keep_features, keep_samples)].shape)

        if args.write_raw:
            bed = bed_template_df.merge(
                pd.DataFrame(subset.X.T, index=subset.var_names, columns=subset.obs[args.sample_id]),
                on="gene_id",
                how="inner",
            )
            bed.to_csv(f"{fname}.{args.agg_method}.raw.bed.gz", index=False, sep="\t", na_rep="NA")

        if np.sum(keep_features) == 0:
            print("No features left, skipping")
            columns = bed_template_df.columns.tolist() + subset.obs[args.sample_id][keep_samples].tolist()
            pd.DataFrame(columns=columns).to_csv(
                f"{fname}.{args.agg_method}.inv.bed.gz", index=False, sep="\t", na_rep="NA"
            )
            continue

        if args.normalize_by_peak_length:
            if sparse.issparse(X):
                X = X.multiply(1.0 / peak_length_kb.rehaspe(-1, 1))
            else:
                X = X * (1.0 / peak_length_kb.reshape(-1, 1))

        if args.agg_method == "sum":
            # inv-normlize log2(CPM + 1)
            arr = np.log2(cpm(X[np.ix_(keep_features, keep_samples)], scale_factor=args.scale_factor) + 1)
        elif args.agg_method == "mean":
            arr = np.log2(X[np.ix_(keep_features, keep_samples)] + 1)
        elif args.agg_method.endswith("_sct_pearson"):
            arr = X[np.ix_(keep_features, keep_samples)]
        X = np.apply_along_axis(inv_normlize, axis=1, arr=arr)
        bed = bed_template_df.merge(
            pd.DataFrame(X, index=subset.var_names[keep_features], columns=subset.obs[args.sample_id][keep_samples]),
            on="gene_id",
            how="inner",
        )
        bed.to_csv(f"{fname}.{args.agg_method}.inv.bed.gz", index=False, sep="\t", na_rep="NA")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", type=str, required=True)
    parser.add_argument(
        "--agg-method", type=str, choices=["sum", "mean", "sum_sct_pearson", "mean_sct_pearson"], default="sum"
    )
    parser.add_argument("--h5ad", type=str, required=True)
    parser.add_argument("--raw-h5ad", type=str)
    parser.add_argument("--gtf", type=str)
    parser.add_argument("--bed", type=str)
    parser.add_argument("--sample-id", type=str, required=True)
    parser.add_argument("--obs-key", type=str, default="barcode")
    parser.add_argument("--obs-annot", type=str)
    parser.add_argument("--obs-annot-key", type=str, default="barcode")
    parser.add_argument("--cell-type-annot", type=str, required=True)
    parser.add_argument("--all-cell-annot", type=str)
    parser.add_argument(
        "--min-cells", type=int, default=10, help="Minimum cell counts per sample-cell type pair to keep samples"
    )
    parser.add_argument("--extra-conditions", type=str, nargs="*")
    parser.add_argument("--keep-sample-id", type=str, help="File with sample IDs to keep")
    parser.add_argument("--keep-feature-id", type=str, help="File with gene IDs to keep")
    parser.add_argument("--keep-obs-key", type=str, help="File with IDs in --obs-key to keep")
    parser.add_argument("--min-sample-prop", type=float, default=0.5, help="Minimum sample proportion to keep features")
    parser.add_argument("--min-cell-prop", type=float, default=0.1, help="Minimum cell proportion to keep features")
    parser.add_argument("--scale-factor", type=int, default=10000, help="Scale factor for normalization")
    parser.add_argument("--export-cell-counts", action="store_true")
    parser.add_argument("--normalize-by-peak-length", action="store_true")
    parser.add_argument("--write-raw", action="store_true")
    args = parser.parse_args()

    if args.gtf is None and args.bed is None:
        raise ValueError("Either --gtf or --bed is required")

    if args.extra_conditions is None:
        args.extra_conditions = []

    main(args)
