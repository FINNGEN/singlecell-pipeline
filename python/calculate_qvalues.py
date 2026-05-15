#!/usr/bin/env python3

import argparse
import bgzip
import io
import numpy as np
import pandas as pd
import tensorqtl


def write_bgz(df, out_file, **kwargs):
    with open(out_file, "wb") as raw:
        with bgzip.BGZipWriter(raw) as fh:
            fh.write(df.to_csv(**kwargs).encode("utf-8"))


def main(args):
    cis_df = pd.concat([pd.read_csv(x, sep="\t", index_col=0, compression="gzip") for x in args.cis_df], axis=0)
    tensorqtl.calculate_qvalues(cis_df)

    if "pval_nominal_threshold" not in cis_df.columns:
        cis_df["pval_nominal_threshold"] = np.nan

    write_bgz(cis_df, f"{args.prefix}.cis_qtl_egenes.tsv.gz", index_label="phenotype_id", sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", type=str, required=True)
    parser.add_argument("--cis-df", type=str, nargs="+", required=True)
    args = parser.parse_args()

    main(args)
