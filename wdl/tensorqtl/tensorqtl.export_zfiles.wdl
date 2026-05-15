version 1.0

workflow tensorqtl_export_zfiles {
    input {
        Array[String] chromosomes
        String output_directory_pattern
        File cell_type_list
        String sumstats_pattern
        String acat_pattern
        Float acat_q_threshold = 0.05
        Boolean debug = false
        String docker_suite
        String zones
    }

    Array[String] cell_types = if debug then [read_lines(cell_type_list)[0]] else read_lines(cell_type_list)

    scatter (cell_type in cell_types) {
        String output_directory_stripped = sub(sub(output_directory_pattern, "\\{CELL_TYPE\\}", cell_type), "[/\\s]+$", "")
        File sumstats = sub(sumstats_pattern, "\\{CELL_TYPE\\}", cell_type)
        File acat = sub(acat_pattern, "\\{CELL_TYPE\\}", cell_type)
        String prefix = "~{basename(sumstats, '.cis_qtl_pairs.tsv.gz')}"

        scatter (chrom in chromosomes) {
            String prefix_chrom = "~{prefix}.~{chrom}"
            call export_zfiles {
                input:
                    chrom = chrom,
                    sumstats = sumstats,
                    acat = acat,
                    acat_q_threshold = acat_q_threshold,
                    prefix_chrom = prefix_chrom,
                    docker = docker_suite,
                    zones = zones
            }

            if (export_zfiles.n_zfiles > 0) {
                call copy {
                    input:
                        output_dir = output_directory_stripped,
                        tar = export_zfiles.out_tar,
                        docker = docker_suite,
                        zones = zones
                }
            }
        }
    }
}

task export_zfiles {
    input {
        String chrom
        File sumstats
        File sumstats_tbi = "~{sumstats}.tbi"
        File acat
        Float acat_q_threshold
        String prefix_chrom
        String docker
        String zones
    }

    command <<<
        set -e

        ulimit -s 65536

        tabix -h ~{sumstats} ~{chrom} |  awk '
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
            print
        }
        NR > 1 && $col["slope_se"] != "NA" {
            print
        }
        ' > ~{chrom}.txt

        echo "Finished processing ~{chrom}.txt"
        wc -l ~{chrom}.txt

        cat << "__EOF__" > export_saige_qtl_zfiles.py
        #!/usr/bin/env python3
        import argparse
        import polars as pl


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
            acat_df = pl.read_csv(args.acat, separator="\t", schema_overrides={"qval": pl.Float64})
            cis_nominal_df = pl.read_csv(args.cis_nominal, separator="\t")
            sig_genes = acat_df.filter(pl.col("qval") < ~{acat_q_threshold}).select("phenotype_id").unique()


            filtered_df = cis_nominal_df.filter(
                pl.col("phenotype_id").is_in(sig_genes.get_column("phenotype_id").to_list())
            )
            if filtered_df.height > 0:
                results = filtered_df.group_by("phenotype_id").map_groups(lambda df: write_z_file(df, args.prefix))
                print(f"Wrote z-files for {len(results)} significant genes.")
            else:
                print("No significant genes found.")


        if __name__ == "__main__":
            parser = argparse.ArgumentParser()
            parser.add_argument("--prefix", type=str, required=True)
            parser.add_argument("--cis-nominal", type=str, required=True)
            parser.add_argument("--acat", type=str, required=True)
            args = parser.parse_args()

            main(args)
        __EOF__

        python3 export_saige_qtl_zfiles.py  \
        --prefix ~{prefix_chrom} \
        --acat ~{acat} \
        --cis-nominal ~{chrom}.txt && \
        ls ~{prefix_chrom}.*.z | wc -l > n_zfiles.txt && \
        tar czvf ~{prefix_chrom}.tar.gz ~{prefix_chrom}.*.z && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        Int n_zfiles = read_int("n_zfiles.txt")
        File out_tar = "~{prefix_chrom}.tar.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "32 GB"
        disks: "local-disk 30 HDD"
        zones: zones
        preemptible: 2
    }
}

task copy {
    input {
        String output_dir
        File tar
        String docker
        String zones
    }

    command <<<
        set -e

        ulimit -s 65536

        tar xzvf ~{tar} && \
        ls *.z > zfiles.txt && \
        cat zfiles.txt | gcloud storage cp -I ~{output_dir}/ && \
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
        disks: "local-disk 30 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}