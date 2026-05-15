#!/usr/bin/env python3

import argparse
import bgzip
import numpy as np
import pandas as pd
import polars as pl
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
from pycct import cct

print(f"PyTorch {torch.__version__}")
print(f"Pandas {pd.__version__}")
print(f"Polars {pl.__version__}")

CHROMOSOMES = [f"chr{i}" for i in list(range(1, 23)) + ["X"]]


def write_bgz(df, out_file, **kwargs):
    with open(out_file, "wb") as raw:
        with bgzip.BGZipWriter(raw) as fh:
            if isinstance(df, pl.DataFrame):
                # Convert kwargs to be compatible with polars
                pl_kwargs = {}
                if "sep" in kwargs:
                    pl_kwargs["separator"] = kwargs["sep"]
                if "na_rep" in kwargs:
                    pl_kwargs["null_value"] = kwargs["na_rep"]
                fh.write(df.write_csv(**pl_kwargs).encode("utf-8"))
            else:
                fh.write(df.to_csv(**kwargs).encode("utf-8"))


def read_cis_qtl_pairs(prefix, variant_df):
    # read back is more straightforward since cis.map_nominal(write_stats=False) doesn't return df
    chromosomes = variant_df.chrom.unique()
    dfs = []

    for chrom in chromosomes:
        pl_df = pl.read_parquet(f"{prefix}.cis_qtl_pairs.{chrom}.parquet")

        if len(pl_df) == 0:
            raise ValueError(f"No cis QTL pairs found for chromosome {chrom}.")

        first_variant = pl_df["variant_id"][0]
        # Determine separator (either ':' or '_')
        separator = ":" if ":" in first_variant else "_"

        # Split using the determined separator and extract position
        pl_df = pl_df.with_columns(pl.col("variant_id").str.split(separator).list.get(1).cast(pl.Int64).alias("pos"))
        pl_df = pl_df.sort("pos").drop("pos")
        dfs.append(pl_df)

    cis_nominal_df = pl.concat(dfs) if len(dfs) > 1 else dfs[0]

    return cis_nominal_df


def process_dataframe(df):
    # Copy the dataframe to avoid modifying the original one
    df_processed = df.copy()

    for column in df_processed.columns:
        # Identify binary columns (including boolean)
        if len(df_processed[column].unique()) == 2:
            # Convert to numeric if not already
            if df_processed[column].dtype == "bool":
                df_processed[column] = df_processed[column].astype(int)
            else:
                # Map non-numeric binary columns to 0 and 1
                unique_values = df_processed[column].unique()
                df_processed[column] = df_processed[column].map({unique_values[0]: 0, unique_values[1]: 1})

        # Identify and one-hot encode categorical columns
        elif df_processed[column].dtype == "object" or df_processed[column].dtype.name == "category":
            df_processed = pd.get_dummies(df_processed, columns=[column], dtype=np.int8)

    return df_processed


def apply_cct(group_df):
    pvals = group_df.get_column("pval_nominal").to_numpy()
    if np.all(np.isnan(pvals)):
        pval_acat = np.nan
        min_pval = np.nan
        min_beta = np.nan
        min_variant = None
    else:
        pval_acat = cct(pvals)
        min_idx = group_df.get_column("pval_nominal").arg_min()
        min_pval = group_df.get_column("pval_nominal")[min_idx]
        min_beta = group_df.get_column("slope")[min_idx]
        min_variant = group_df.get_column("variant_id")[min_idx]

    result = pl.DataFrame(
        {
            "phenotype_id": [group_df.get_column("phenotype_id")[0]],
            "pval_acat": [pval_acat],
            "pval_nominal": [min_pval],
            "slope": [min_beta],
            "variant_id": [min_variant],
        }
    )
    return result


def write_z_file(gene_df, prefix):
    gene = gene_df.get_column("phenotype_id")[0]

    # Determine separator from the first variant
    first_variant = gene_df["variant_id"][0]
    separator = ":" if ":" in first_variant else "_"

    # Split the variant ID and extract components
    gene_df = gene_df.with_columns(pl.col("variant_id").str.split(separator).alias("parts"))

    # Extract each component from the parts list
    gene_df = gene_df.with_columns(
        pl.col("parts").list.get(0).alias("chromosome"),
        pl.col("parts").list.get(1).cast(pl.Int64).alias("position"),
        pl.col("parts").list.get(2).alias("allele1"),
        pl.col("parts").list.get(3).alias("allele2"),
        (0.5 - (0.5 - pl.col("af")).abs()).alias("maf"),
    )

    gene_df = gene_df.rename({"variant_id": "rsid", "slope": "beta", "slope_se": "se", "pval_nominal": "p"})
    gene_df = gene_df.select(
        ["rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se", "p"]
    ).drop_nulls()
    gene_df.write_csv(f"{prefix}.{gene}.z", separator=" ", null_value="NA", include_header=True)
    return pl.DataFrame({"gene": [gene], "status": ["success"]})


def main(args):
    # load phenotypes and covariates
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(args.expression_bed)
    covariates_df = pd.read_csv(args.covar, sep="\t", index_col=0)
    covariates_df = process_dataframe(covariates_df)

    # filter to chrom
    phenotype_df = phenotype_df.loc[phenotype_pos_df.chr.isin(args.chromosomes), :]
    phenotype_pos_df = phenotype_pos_df.loc[phenotype_pos_df.chr.isin(args.chromosomes), :]

    # keep only samples that exist in bed and covar
    pheno_samples = phenotype_df.columns.tolist()
    covar_samples = covariates_df.index.tolist()
    if len(pheno_samples) != len(covar_samples):
        keep = np.intersect1d(pheno_samples, covar_samples).tolist()
        n_samples = len(keep)
        print(
            f"# samples don't match, keeping {n_samples} intersected samples -- pheno:{len(pheno_samples)}, covar:{len(covar_samples)}"
        )
        phenotype_df = phenotype_df.loc[:, keep]
        covariates_df = covariates_df.loc[keep, :]
    else:
        keep = None
        n_samples = len(pheno_samples)

    # compute a new MAF threshold
    maf_threshold = np.max([args.mac_threshold / (2 * n_samples), args.maf_threshold])

    # PLINK reader for genotypes
    genotype_df, variant_df = genotypeio.load_genotypes(
        args.pfile, select_samples=keep, dosages=(not args.dominant and not args.recessive)
    )

    if (args.dominant or args.recessive) and (maf_threshold > 0):
        maf = 0.5 - np.abs(0.5 - genotype_df.sum(axis=1) / (2 * genotype_df.shape[1]))
        genotype_df = genotype_df.loc[maf > maf_threshold, :]
        variant_df = variant_df.loc[maf > maf_threshold, :]
        maf_threshold = 0

    if args.dominant:
        print("Dominant model is specified. Genotypes are now coded 0, 1, 1.")
        genotype_df = (genotype_df > 0).astype(int)
    if args.recessive:
        print("Recessive model is specified. Genotypes are now coded 0, 0, 1.")
        genotype_df = (genotype_df == 2).astype(int)

    # always use chr prefix
    if not variant_df.chrom.iloc[0].startswith("chr"):
        variant_df["chrom"] = "chr" + variant_df.chrom.astype(str)
    # map PAR1/2 to X
    variant_df.loc[variant_df.chrom.str.endswith(("PAR1", "PAR2")), "chrom"] = "chrX"

    # cis-QTL: nominal p-values for all variant-phenotype pairs
    if args.mode == "cis_nominal":
        cis.map_nominal(
            genotype_df,
            variant_df,
            phenotype_df,
            phenotype_pos_df,
            args.prefix,
            covariates_df=covariates_df,
            maf_threshold=maf_threshold,
            window=args.window,
        )

        if args.write_tsv:
            # returns polars DF
            cis_nominal_df = read_cis_qtl_pairs(args.prefix, variant_df)
            write_bgz(cis_nominal_df, f"{args.prefix}.cis_qtl_pairs.tsv.gz", sep="\t", na_rep="NA")

            if args.acat:
                acat_df = cis_nominal_df.group_by("phenotype_id").map_groups(apply_cct)
                write_bgz(acat_df, f"{args.prefix}.cis_qtl_acat.tsv.gz", sep="\t", na_rep="NA")

            if args.write_z_files:
                sig_genes = cis_nominal_df.filter(pl.col("pval_nominal") < 5e-8).select("phenotype_id").unique()

                if len(sig_genes) > 0:
                    filtered_df = cis_nominal_df.filter(
                        pl.col("phenotype_id").is_in(sig_genes.get_column("phenotype_id").to_list())
                    )
                    results = filtered_df.group_by("phenotype_id").map_groups(lambda df: write_z_file(df, args.prefix))
                    print(f"Wrote z-files for {len(results)} significant genes.")

    elif args.mode == "cis_interaction":
        interaction_df = pd.read_csv(args.interaction, sep="\t", index_col=0)
        interaction_df = process_dataframe(interaction_df).dropna()

        pheno_samples = phenotype_df.columns.tolist()
        interaction_samples = interaction_df.index.tolist()
        if len(pheno_samples) != len(interaction_samples):
            keep = np.intersect1d(pheno_samples, interaction_samples).tolist()
            print(
                f"# samples don't match, keeping {len(keep)} intersected samples -- interaction:{len(interaction_samples)}"
            )
            phenotype_df = phenotype_df.loc[:, keep]
            covariates_df = covariates_df.loc[keep, :]
            interaction_df = interaction_df.loc[keep, :]

        cis.map_nominal(
            genotype_df,
            variant_df,
            phenotype_df,
            phenotype_pos_df,
            args.prefix,
            covariates_df=covariates_df,
            interaction_df=interaction_df,
            run_eigenmt=True,
            maf_threshold=maf_threshold,
            window=args.window,
            write_top=True,
            write_stats=args.write_tsv,
        )

        if args.write_tsv:
            cis_nominal_df = read_cis_qtl_pairs(args.prefix, variant_df)
            write_bgz(cis_nominal_df, f"{args.prefix}.cis_qtl_pairs.tsv.gz", index=False, sep="\t", na_rep="NA")

    elif args.mode == "cis":
        cis_df = cis.map_cis(
            genotype_df,
            variant_df,
            phenotype_df,
            phenotype_pos_df,
            covariates_df=covariates_df,
            maf_threshold=maf_threshold,
            window=args.window,
            seed=args.seed,
        )

        if args.chromosomes == CHROMOSOMES:
            tensorqtl.calculate_qvalues(cis_df)
            if "pval_nominal_threshold" not in cis_df.columns:
                cis_df["pval_nominal_threshold"] = np.nan

        write_bgz(cis_df, f"{args.prefix}.cis_qtl_egenes.tsv.gz", index_label="phenotype_id", sep="\t", na_rep="NA")
    elif args.mode == "trans":
        trans_df = trans.map_trans(
            genotype_df,
            phenotype_df,
            covariates_df=covariates_df,
            return_sparse=True,
            pval_threshold=args.pval_threshold,
            maf_threshold=maf_threshold,
            batch_size=args.batch_size,
        )

        # remove cis-variants within 5 Mbp
        trans_df = trans.filter_cis(trans_df, phenotype_pos_df, variant_df, window=5000000)
        trans_df.to_parquet(f"{args.prefix}.trans_qtl_pairs.parquet")

        if args.write_tsv:
            write_bgz(trans_df, f"{args.prefix}.trans_qtl_pairs.tsv.gz", index=False, sep="\t", na_rep="NA")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", type=str, required=True)
    parser.add_argument("--type", type=str, choices=["gex", "atac"], default="gex")
    parser.add_argument(
        "--mode", type=str, choices=["cis_nominal", "cis_interaction", "cis", "trans"], default="cis_nominal"
    )
    parser.add_argument("--expression-bed", type=str, required=True)
    parser.add_argument("--covar", type=str, required=True)
    parser.add_argument("--interaction", type=str)
    parser.add_argument("--pfile", type=str, required=True)
    parser.add_argument("--chromosomes", type=str, choices=CHROMOSOMES, default=CHROMOSOMES, nargs="+")
    parser.add_argument("--window", type=int, default=1000000)
    parser.add_argument("--maf-threshold", type=float, default=0)
    parser.add_argument("--mac-threshold", type=int, default=0)
    parser.add_argument("--pval-threshold", type=float, default=1e-5)
    parser.add_argument("--seed", type=int, default=123456)
    parser.add_argument("--batch-size", type=int, default=20000)
    parser.add_argument("--dominant", action="store_true")
    parser.add_argument("--recessive", action="store_true")
    parser.add_argument("--write-tsv", action="store_true")
    parser.add_argument("--write-z-files", action="store_true")
    parser.add_argument("--acat", action="store_true")
    args = parser.parse_args()

    if args.chromosomes == "all":
        args.chromosomes = CHROMOSOMES

    if args.mode == "cis_interaction" and args.interaction is None:
        parser.error("--interaction is required for mode 'cis_interaction'")

    if args.dominant and args.recessive:
        parser.error("Cannot specify both --dominant and --recessive")

    main(args)
