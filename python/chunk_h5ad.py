#!/usr/bin/env python3

import argparse
import bgzip
import anndata as ad
import numpy as np
import pandas as pd

CHROMOSOMES = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y", "M"]]


def write_bgz(df, out_file, **kwargs):
    with open(out_file, "wb") as raw:
        with bgzip.BGZipWriter(raw) as fh:
            fh.write(df.to_csv(**kwargs).encode("utf-8"))


def chunk_by_var_names(adata, var_names, fname, omit_barcode=False):
    idx = adata.var_names.isin(var_names)
    if omit_barcode:
        kwargs = {"index": False}
    else:
        kwargs = {"index_label": "barcode"}
    if np.sum(idx) > 0:
        write_bgz(adata[:, idx].to_df(), out_file=fname, sep="\t", **kwargs)


def main(args):
    adata = ad.read_h5ad(args.h5ad, backed="r")

    if args.chunk_by == "chromosome":
        for chrom in args.chromosomes:
            fname = f"{args.prefix}.{chrom}.h5ad"
            if args.type == "gex":
                adata[:, adata.var_names[adata.var.chrom == chrom]].copy(fname)
            elif args.type == "atac":
                adata[:, adata.var_names.str.startswith(f"{chrom}-")].copy(fname)
    elif args.chunk_by == "gene_id" and args.type == "gex":
        for gene_id in args.gene_ids:
            chunk_by_var_names(
                adata, var_names=[gene_id], fname=f"{gene_id}.{args.type}.tsv.gz", omit_barcode=args.omit_barcode
            )
    elif args.chunk_by == "gene_id" and args.type == "atac":
        bed = pd.read_csv(args.tss_bed, sep="\t")
        bed.columns = ["chr", "start", "end", "gene_id"]
        bed = bed.loc[bed.gene_id.isin(args.gene_ids), :]

        peaks = adata.var_names.to_series().str.split("-", expand=True)
        peaks.columns = ["chr", "start", "end"]
        peaks["gene_id"] = peaks.index
        peaks["pos"] = (peaks.start.astype(int) + peaks.end.astype(int)) // 2

        for _, row in bed.iterrows():
            chrom, start, end = row.chr, row.start - args.window, row.end + args.window
            gene_id = row.gene_id
            peak_ids = peaks.loc[(peaks.chr == chrom) & (start < peaks.pos) & (peaks.pos < end)].gene_id

            chunk_by_var_names(
                adata, var_names=peak_ids, fname=f"{gene_id}.{args.type}.tsv.gz", omit_barcode=args.omit_barcode
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", type=str, required=True)
    parser.add_argument("--type", type=str, choices=["gex", "atac"], default="gex")
    parser.add_argument("--chunk-by", type=str, choices=["chromosome", "gene_id"], default="chromosome")
    parser.add_argument("--h5ad", type=str, required=True)
    parser.add_argument("--chromosomes", type=str, choices=CHROMOSOMES, default=CHROMOSOMES, nargs="+")
    parser.add_argument("--gene-ids", type=str, nargs="*")
    parser.add_argument("--tss-bed", type=str)
    parser.add_argument("--window", type=int, default=500000)
    parser.add_argument("--omit-barcode", action="store_true")
    args = parser.parse_args()

    main(args)
