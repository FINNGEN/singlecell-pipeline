version development

task preprocess_pseudobulk {
    input {
        File gex_features
        File expression_bed
        File nonzero_prop
        File covar
        String sample_id
        Int num_genes_per_chunk
        Boolean filter_genes
        Float nonzero_prop_threshold = 0
        File? subset_genes
        String prefix
        String chunk_list
        String docker
        String zones
    }

    command <<<
        set -e

        python3 << "__EOF__"
        import csv
        import math
        import pandas as pd
        df = pd.read_csv("~{expression_bed}", sep="\t")
        df_exp = df.iloc[:, 4:].T
        df_exp.columns = df.gene_id
        df_nonzero = pd.read_csv("~{nonzero_prop}", sep="\t")[["gene_id", "nonzero_prop"]]
        df_covar = pd.read_csv("~{covar}", sep="\t")

        out = df_covar.merge(df_exp, how="inner", left_on="~{sample_id}", right_index=True)
        out.to_csv("~{prefix}.SAIGE.pheno.txt.gz", sep="\t", index=False)

        features = pd.read_csv(
            "~{gex_features}",
            names=["gene_id", "symbol", "type", "chrom", "start", "end", "gene_type"],
            header=None,
            sep="\t",
        )
        df = df.merge(features, on="gene_id").merge(df_nonzero, on="gene_id")

        nonzero_prop_threshold = ~{nonzero_prop_threshold}
        if nonzero_prop_threshold > 0:
            # filter to genes with nonzero prop > threshold
            df = df.loc[df.nonzero_prop > nonzero_prop_threshold, :]
        if ~{true='True' false='False' filter_genes}:
            # filter to protein-coding genes
            df = df.loc[(df.gene_type == "protein_coding"), :]
        if ~{true='True' false='False' defined(subset_genes)}:
            # subset to genes in the provided list
            subset_genes = pd.read_csv("~{subset_genes}", header=None).iloc[:, 0]
            df = df.loc[df.gene_id.isin(subset_genes), :]

        # for peaks, we need to sanitize "-" in gene_id to "_" for SAIGE-QTL
        df["gene_id"] = df["gene_id"].str.replace("-", "_", regex=False)

        # chunk genes
        num_genes_per_chunk = ~{num_genes_per_chunk}
        list_of_lists = []
        for chromosome, group in df.groupby("chr"):
            gene_ids = group["gene_id"].tolist()
            total_genes = len(gene_ids)
            min_chunks_needed = math.ceil(total_genes / num_genes_per_chunk)
            genes_per_chunk_adjusted = math.ceil(total_genes / min_chunks_needed)
            genes_per_chunk_adjusted = min(genes_per_chunk_adjusted, total_genes)
            genes_per_chunk_adjusted = min(genes_per_chunk_adjusted, num_genes_per_chunk)

            for i in range(0, len(gene_ids), genes_per_chunk_adjusted):
                chunk = gene_ids[i:i+genes_per_chunk_adjusted]
                list_of_lists.append([chromosome] + chunk)

        # first column is chromosome - then gene_ids
        padded_list_of_lists = [row + ["_NULL_"] * (num_genes_per_chunk + 1 - len(row)) for row in list_of_lists]
        with open("~{prefix}.chunks.txt", "w", newline="") as f:
            csv.writer(f, delimiter="\t").writerows(padded_list_of_lists)
        __EOF__

        if [[ ~{defined(subset_genes)} == "false" ]]
        then
            gcloud storage cp ~{prefix}.chunks.txt ~{chunk_list}
        fi
    >>>

    output {
        File out = prefix + ".SAIGE.pheno.txt.gz"
        File out_chunks = prefix + ".chunks.txt"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "16 GB"
        disks: "local-disk 5 HDD"
        zones: zones
        preemptible: 2
    }
}

task preprocess_sc {
    input {
        String h5ad_pattern
        Array[String] chunk
        File covar
        Array[String] covariates
        File? cell_covar
        File obs
        File barcode
        String obs_sample_id
        String sample_id
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e

        if [[ ~{defined(cell_covar)} == "true" ]]
        then
            if [[ "~{cell_covar}" == *.gz ]]
            then
                mv ~{cell_covar} cell_covar.txt.gz
            else
                gzip -c ~{cell_covar} > cell_covar.txt.gz
            fi
        fi

        python3 << "__EOF__"
        import numpy as np
        import pandas as pd
        from gcs_anndata import GCSAnnData
        from google.cloud import storage

        # sample-level covar
        df_covar = pd.read_csv("~{covar}", sep="\t")

        # obs covariates
        covariates = "~{sep=',' covariates}".split(",")
        obs_cols = ["barcode", "~{obs_sample_id}"]
        if "log_total_counts" in covariates:
            total_counts_col = "total_counts"
        elif "log_total_counts_wo_highly_expressed" in covariates:
            total_counts_col = "total_counts_wo_highly_expressed"
        else:
            total_counts_col = None
        obs_cols += [total_counts_col] if total_counts_col else []
        obs_cols += list(set(covariates) - set(obs_cols) - set(df_covar.columns) - set([f"log_{total_counts_col}"]))
        df_obs = pd.read_csv("~{obs}", sep="\t", usecols=obs_cols)
        if total_counts_col is not None:
            df_obs[[f"log_{total_counts_col}"]] = np.log(df_obs[[total_counts_col]])

        chunk = "~{sep=',' chunk}".split(",")
        chromosome = chunk[0]
        gene_ids = [x for x in chunk[1:] if x != "_NULL_"]
        adata = GCSAnnData("~{h5ad_pattern}".replace("{CHROM}", chromosome))
        df_exp = pd.DataFrame(adata.get_columns(gene_ids).toarray(), index=df_obs.barcode, columns=gene_ids).reset_index()
        barcodes = pd.read_csv("~{barcode}", header=None).iloc[:, 0]

        if ~{if defined(cell_covar) then "True" else "False"}:
            # cell-level covar
            df_cell_covar = pd.read_csv("cell_covar.txt.gz", sep="\t")
            df_obs = df_obs.merge(df_cell_covar, how="inner", on="barcode")
            overlapping_covars = df_covar.columns.intersection(df_cell_covar.columns)
            df_covar = df_covar.drop(columns=overlapping_covars)

        # export a phenofile for SAIGE-QTL
        df_exp = df_exp.merge(df_obs, how="inner", on="barcode")
        out = df_covar.merge(df_exp, how="inner", left_on="~{sample_id}", right_on="~{obs_sample_id}")
        out = out.loc[out.barcode.isin(barcodes), :]
        # for peaks, we need to sanitize "-" in gene_id to "_" for SAIGE-QTL
        out.columns = out.columns.str.replace("-", "_", regex=False)
        out.to_csv("~{prefix}.SAIGE.pheno.txt.gz", sep="\t", index=False)

        # export PEER columns
        out.columns[out.columns.str.startswith("PEER")].to_frame().to_csv("peer_covariates.txt", index=False, header=False)

        # export number of cells
        with open("n_cells.txt", "w") as f:
            f.write(str(len(out.index)))

        # list gene_ids to process
        gene_ids = df_exp.columns.to_series()
        gene_ids = gene_ids[gene_ids.str.startswith('ENSG') | gene_ids.str.match(r'^chr[\dX]+[-_]\d+[-_]\d+$')]
        gene_ids = gene_ids.str.replace("-", "_", regex=False)
        gene_ids.to_csv("~{prefix}.gene_ids.txt", index=False, header=False)
        __EOF__
    >>>

    output {
        File out = prefix + ".SAIGE.pheno.txt.gz"
        File out_gene_ids = prefix + ".gene_ids.txt"
        Int n_cells = read_int("n_cells.txt")
        Array[String] peer_covariates = read_lines("peer_covariates.txt")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{ceil(200 * size(obs, 'GB')) + 10} GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task check_step1_null {
    input {
        Array[String] gene_ids
        String out_dir
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e

        python3 << "__EOF__"
        import pandas as pd
        from google.cloud import storage

        # list gene_ids to process
        gene_ids = "~{sep=',' gene_ids}".split(",")
        rda_files = [f"~{out_dir}/step1/~{prefix}.{x}.rda" for x in gene_ids]
        vr_files = [f"~{out_dir}/step1/~{prefix}.{x}.varianceRatio.txt" for x in gene_ids]

        storage_client = storage.Client()
        bucket = storage_client.get_bucket("~{out_dir}".split("/")[2])
        rda_exists = [bucket.blob(x.split("/", 3)[-1]).exists() for x in rda_files]
        vr_exists = [bucket.blob(x.split("/", 3)[-1]).exists() for x in vr_files]

        df = pd.DataFrame(
            {
                "gene_ids": gene_ids,
                "rda_exists": rda_exists,
                "vr_exists": vr_exists
            }
        )
        df = df[~df.rda_exists | ~df.vr_exists]
        df.gene_ids.to_csv("~{prefix}.gene_ids.txt", index=False, header=False)
        __EOF__
    >>>

    output {
        File out_gene_ids = prefix + ".gene_ids.txt"
    }

    meta {
        volatile: true
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 5 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task step1_null {
    input {
        String trait_type
        String bfile
        File bed = bfile + ".bed"
        File bim = bfile + ".bim"
        File fam = bfile + ".fam"
        File pheno
        String sample_id
        String gene_id
        Array[String] covariates
        Array[String] sample_covariates
        Array[String] categorical_covariates = []
        Array[String] dynamic_covariates = []
        String offset = ""
        Int n_cells
        String out_dir
        String prefix
        Boolean isCovariateOffset = false
        Boolean isStoreSigma = true
        Boolean useGRMtoFitNULL = false
        Float tol = 0.00001
        Int maxiter = 20
        Int maxiterPCG = 50000
        String tauInit = "0,0"
        String docker
        String zones
    }

    Int task_size = ceil(1.5 * n_cells / 100000 + 2 * length(dynamic_covariates))
    Int memory = if (task_size < 8) then 8 else task_size
    String offset_arg = if offset != "" then "--offsetCol=~{offset}" else ""

    command <<<
        set -e

        n_cpu=$(grep -c ^processor /proc/cpuinfo)

        step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL=FALSE \
        --useGRMtoFitNULL=~{true='TRUE' false='FALSE' useGRMtoFitNULL} \
        --phenoFile=~{pheno} \
        --phenoCol=~{gene_id} \
        --covarColList=~{sep="," covariates} \
        --sampleCovarColList=~{sep="," sample_covariates} \
        ~{true='--qCovarColList=' false='' length(categorical_covariates) > 0}~{sep="," categorical_covariates} \
        ~{true='--dynamicCovarColList=' false='' length(dynamic_covariates) > 0}~{sep="," dynamic_covariates} \
       --tauInit=~{tauInit} \
        ~{offset_arg} \
        --sampleIDColinphenoFile=~{sample_id} \
        --traitType=~{trait_type} \
        --outputPrefix=~{prefix} \
        --skipVarianceRatioEstimation=FALSE \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=~{true='TRUE' false='FALSE' isCovariateOffset} \
        --isCovariateTransform=TRUE \
        --skipModelFitting=FALSE \
        --tol=~{tol} \
        --maxiter=~{maxiter} \
        --maxiterPCG=~{maxiterPCG} \
        --nThreads=$n_cpu \
        --bedFile=~{bed} \
        --bimFile=~{bim} \
        --famFile=~{fam} \
        --IsOverwriteVarianceRatioFile=TRUE \
        --isStoreSigma=~{true='TRUE' false='FALSE' isStoreSigma} \
        --isShrinkModelOutput=TRUE \
        --cellIDColinphenoFile="barcode" \
        --isExportResiduals=TRUE && \
        bgzip ~{prefix}.residuals.txt && \
        bgzip ~{prefix}.sample.residuals.txt && \
        gcloud storage cp \
        ~{prefix}.rda \
        ~{prefix}.residuals.txt.gz \
        ~{prefix}.sample.residuals.txt.gz \
        ~{prefix}.varianceRatio.txt \
        ~{out_dir}/step1/ && \
        touch _SUCCESS
    >>>

    output {
        String out_gene_id = gene_id
        File out_modelfile = "~{prefix}.rda"
        File out_varianceratio = "~{prefix}.varianceRatio.txt"
        File out_residuals = "~{prefix}.residuals.txt.gz"
        File out_sample_residuals = "~{prefix}.sample.residuals.txt.gz"
        File success = "_SUCCESS"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{memory} GB"
        disks: "local-disk 100 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task rechunk_step2_list {
    input {
        File chunk_list
        Int step2_chunk_size
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e

        python3 << "__EOF__"
        import pandas as pd

        chunk_size = ~{step2_chunk_size}
        df = pd.read_csv("~{chunk_list}", sep="\t", header=None)

        rows = []
        for _, row in df.iterrows():
            chrom = row.iloc[0]
            gene_ids = [x for x in row.iloc[1:] if pd.notna(x) and x != "_NULL_"]
            for i in range(0, len(gene_ids), chunk_size):
                rows.append([chrom] + gene_ids[i:i+chunk_size])

        pd.DataFrame(rows).fillna("_NULL_").to_csv("~{prefix}.rechunked.txt", sep="\t", index=False, header=False)
        __EOF__
    >>>

    output {
        File out_chunk_list = "~{prefix}.rechunked.txt"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 1 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task generate_multi_gene_file {
    input {
        String step1_rda_pattern
        String step1_vr_pattern
        Array[String] chunk
        Array[String] bgenfiles
        File expression_bed
        Int cis_window = 1000000
        Boolean trans = false
        String prefix
        String tmp_bucket
        Boolean overwrite
        Boolean skip_non_existing
        String docker
        String zones
    }

    command <<<
        set -e

        python3 << "__EOF__"
        import os
        import re
        import sqlite3
        import pandas as pd
        from google.cloud import storage

        chunk = "~{sep=',' chunk}".split(",")
        chromosome = chunk[0]
        gene_ids = [x for x in chunk[1:] if x != "_NULL_"]

        # check existence of step1 and step2 outputs
        rda_files = ["~{step1_rda_pattern}".replace("{GENE_ID}", x) for x in gene_ids]
        vr_files = ["~{step1_vr_pattern}".replace("{GENE_ID}", x) for x in gene_ids]
        step2_gene_id_file = "~{tmp_bucket}/~{prefix}.gene_ids.txt"
        step2_ht_path = step2_gene_id_file.replace(".gene_ids.txt", ".ht")

        storage_client = storage.Client()
        step1_bucket = storage_client.get_bucket("~{step1_rda_pattern}".split("/")[2])
        step2_bucket = storage_client.get_bucket("~{tmp_bucket}".split("/")[2])
        rda_exists = [step1_bucket.blob(x.split("/", 3)[-1]).exists() for x in rda_files]
        vr_exists = [step1_bucket.blob(x.split("/", 3)[-1]).exists() for x in vr_files]
        step2_ht_exits = all([step2_bucket.blob(x.split("/", 3)[-1]).exists() for x in [step2_gene_id_file, f"{step2_ht_path}/_SUCCESS"]])

        df = pd.DataFrame(
            {
                "gene_ids": gene_ids,
                "modelfiles": rda_files,
                "varianceratiofiles": vr_files,
                "rda_exists": rda_exists,
                "vr_exists": vr_exists
            }
        )

        print(df[~df.rda_exists | ~df.vr_exists])
        df[~df.rda_exists | ~df.vr_exists].gene_ids.to_csv("~{prefix}.skip_gene_ids.txt", index=False, header=False)
        if (not ~{true='True' false='False' skip_non_existing}) and (not all(rda_exists) or not all(vr_exists)):
            raise ValueError("Some step1 files do not exist")
        if (not any(rda_exists) or not any(vr_exists)):
            raise ValueError("All step1 files do not exist")

        step2_ht_out = ""
        if (not ~{true='True' false='False' overwrite}) and step2_ht_exits:
            ht_gene_ids = set(pd.read_csv(step2_gene_id_file, header=None).iloc[:, 0])
            if set(df.gene_ids) == ht_gene_ids:
                print(f"Step2 ht already exists: {step2_ht_path}")
                step2_ht_out = step2_ht_path

        with open("step2_ht.txt", "w") as f:
            f.write(step2_ht_out)

        df = df[df.rda_exists & df.vr_exists]
        df.gene_ids.to_csv("~{prefix}.gene_ids.txt", index=False, header=False)
        df.modelfiles.to_csv("~{prefix}.modelfiles.txt", index=False, header=False)
        df.varianceratiofiles.to_csv("~{prefix}.varianceratiofiles.txt", index=False, header=False)

        # compute union cis region from expression_bed
        cis_window = ~{cis_window} + 1000
        bed = pd.read_csv("~{expression_bed}", sep="\t", usecols=["chr", "start", "end", "gene_id"])
        bed = bed.rename(columns={"#chr": "chr"})
        # for peaks, we need to sanitize "-" in gene_id to "_" for SAIGE-QTL
        bed["gene_id"] = bed["gene_id"].str.replace("-", "_", regex=False)

        bed_chunk = bed[bed.gene_id.isin(df.gene_ids.tolist())]
        cis_start = max(1, int(bed_chunk.start.min()) - cis_window)
        cis_end = int(bed_chunk.end.max()) + cis_window

        with open("~{prefix}.cis_region.txt", "w") as f:
            f.write(f"{chromosome}\t{cis_start}\t{cis_end}\n")

        # identify bgen files overlapping the cis region
        all_bgenfiles = "~{sep=',' bgenfiles}".split(",")
        trans = ~{true='True' false='False' trans}

        def extract_bgen_chrom(bgen_path):
            """Extract chromosome from bgen filename.
            e.g. FG_EA5_batch1_5.chr1.00.bgen -> chr1
            """
            basename = os.path.basename(bgen_path).replace(".bgen", "")
            basename = re.sub(r'^.*?\.(?=chr)', '', basename)
            chrom = re.sub(r'\.[0-9]+$', '', basename)
            return chrom

        def extract_bgen_chunk_index(bgen_path):
            """Extract chunk index from bgen filename.
            e.g. FG_EA5_batch1_5.chr1.00.bgen -> 0
            """
            basename = os.path.basename(bgen_path).replace(".bgen", "")
            m = re.search(r'\.(\d+)$', basename)
            return int(m.group(1)) if m else 0

        if trans:
            overlapping_bgens = all_bgenfiles
        else:
            # filter to bgens matching the chunk's chromosome, sorted by chunk index
            chrom_bgens = sorted(
                [b for b in all_bgenfiles if extract_bgen_chrom(b) == chromosome],
                key=extract_bgen_chunk_index
            )

            # download .bgi files from GCS and query for cis-region overlap
            # bgens are sorted by position, so overlapping bgens are contiguous:
            # skip bgens entirely before cis region, collect overlapping ones,
            # stop at the first bgen entirely after cis region
            os.makedirs("bgi_tmp", exist_ok=True)
            overlapping_bgens = []
            found_overlap = False
            for bgen_path in chrom_bgens:
                bgi_gcs_path = f"{bgen_path}.bgi"
                local_bgi = os.path.join("bgi_tmp", os.path.basename(bgi_gcs_path))

                bucket_name = bgi_gcs_path.split("/")[2]
                blob_path = bgi_gcs_path.split("/", 3)[-1]
                bucket = storage_client.get_bucket(bucket_name)
                bucket.blob(blob_path).download_to_filename(local_bgi)

                conn = sqlite3.connect(local_bgi)
                cursor = conn.cursor()
                cursor.execute(
                    "SELECT COUNT(*) FROM Variant WHERE position >= ? AND position <= ?",
                    (cis_start, cis_end)
                )
                count = cursor.fetchone()[0]
                conn.close()

                if count > 0:
                    found_overlap = True
                    overlapping_bgens.append(bgen_path)
                    print(f"  {os.path.basename(bgen_path)}: {count} variants in cis region")
                else:
                    print(f"  {os.path.basename(bgen_path)}: no variants in cis region, skipping")
                    if found_overlap:
                        # past the cis region, no need to check remaining bgens
                        break

            print(f"Cis region: {chromosome}:{cis_start}-{cis_end}")
            print(f"Chromosome bgens: {len(chrom_bgens)}, overlapping cis region: {len(overlapping_bgens)}")

        with open("~{prefix}.bgenfiles.txt", "w") as f:
            for b in overlapping_bgens:
                f.write(b + "\n")
        __EOF__
    >>>

    output {
        File out_gene_ids = prefix + ".gene_ids.txt"
        File out_modelfiles = prefix + ".modelfiles.txt"
        File out_varianceratiofiles = prefix + ".varianceratiofiles.txt"
        File out_skip_gene_ids = prefix + ".skip_gene_ids.txt"
        File out_cis_region = prefix + ".cis_region.txt"
        Array[String] out_bgenfiles = read_lines(prefix + ".bgenfiles.txt")
        String out_ht = read_string("step2_ht.txt")
    }

    meta {
        volatile: true
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task step2_test {
    input {
        Array[String] gene_ids
        Array[File] modelfiles
        Array[File] varianceratiofiles
        File bgen
        File bgi = bgen + ".bgi"
        File samplefile = sub(bgen, "\\.bgen$", ".sample")
        File cis_region
        Float min_maf = 0
        Int min_mac = 20
        Int markers_per_chunk = 10000
        Boolean is_gxe = false
        Float pval_cutoff_for_gxe = 1
        Boolean is_permute_e = false
        Boolean is_permute_ginge = false
        String prefix
        String docker
        String zones
    }

    String pval_cutoff_for_gxe_arg = if is_gxe then "--pval_cutoff_for_gxe=~{pval_cutoff_for_gxe}" else ""
    String permute_e_arg = if is_permute_e then "--is_permute_e=TRUE" else ""
    String permute_ginge_arg = if is_permute_ginge then "--is_permute_ginge=TRUE" else ""
    Int memory = if is_gxe then ceil(2 * size(modelfiles, 'GB')) + 10 else ceil(6 * size(modelfiles, 'GB')) + 10

    command <<<
        set -e

        python3 << "__EOF__"
        import pandas as pd
        pd.DataFrame(
            {
                "gene_ids": "~{sep=',' gene_ids}".split(","),
                "modelfiles": "~{sep=',' modelfiles}".split(","),
                "varianceratiofiles": "~{sep=',' varianceratiofiles}".split(",")
            }
        ).to_csv("multi_gene_file.txt", sep="\t", index=False, header=False)
        __EOF__

        step2_tests_qtl.R \
        --bgenFile=~{bgen} \
        --bgenFileIndex=~{bgi} \
        --sampleFile=~{samplefile} \
        --AlleleOrder=ref-first \
        --rangestoIncludeFile=~{cis_region} \
        --minMAF=~{min_maf} \
        --minMAC=~{min_mac} \
        --SAIGEOutputFile=~{prefix}.SAIGE.step2 \
        --LOCO=FALSE \
        --SPAcutoff=2 \
        --markers_per_chunk=~{markers_per_chunk} \
        --is_fastTest=TRUE \
        --pval_cutoff_for_fastTest=0.05 \
        --GMMATmodel_varianceRatio_multiTraits_File=multi_gene_file.txt \
        ~{pval_cutoff_for_gxe_arg} \
        ~{permute_e_arg} \
        ~{permute_ginge_arg}

        # SAIGE-QTL omits the _{phenotype} suffix when only one trait is tested,
        # which breaks the output glob and downstream split("step2_") parsing.
        if [ ~{length(gene_ids)} -eq 1 ]; then
            mv "~{prefix}.SAIGE.step2" "~{prefix}.SAIGE.step2_~{gene_ids[0]}"
        fi
    >>>

    output {
        Array[File] out_results = glob(prefix + ".SAIGE.step2_*")
        File out_index = prefix + ".SAIGE.step2.index"
    }

    runtime {
        docker: docker
        cpu: 1
        # 45 for mk, 80 for mk2, 200 for mk2 & peer45, 8 for w/ total counts
        memory: "~{memory} GB"
        disks: "local-disk ~{ceil(size(modelfiles, 'GB')) + 10} HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task write_step2_ht {
    input {
        Array[String] gene_ids
        Array[String] result_paths
        File expression_bed
        String cell_type
        Array[String] dynamic_covariates = []
        Boolean is_gxe = length(dynamic_covariates) > 0
        String prefix
        String tmp_bucket
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
        from hail.utils import new_temp_file

        hl.default_reference("GRCh38")
        cis_window = 1000000
        ref = hl.default_reference()
        contig_lengths = hl.dict(ref.lengths)
        # gxe
        is_gxe = ~{true='True' false='False' is_gxe}
        dynamic_cols = ["Beta_ge", "seBeta_ge", "pval_ge"]
        dynamic_covariates = "~{sep=',' dynamic_covariates}".split(",")

        ht_interval = hl.import_table("~{expression_bed}", impute=True, force=True)
        # for peaks
        ht_interval = ht_interval.annotate(gene_id=ht_interval.gene_id.replace("-", "_"))
        ht_interval = ht_interval.key_by(
            interval=hl.locus_interval(
                ht_interval.chr,
                hl.max(1, ht_interval.start - cis_window),
                hl.min(ht_interval.end + cis_window, contig_lengths[ht_interval.chr]),
            )
        ).select("gene_id").cache()

        paths = "~{sep=',' result_paths}".split(",")
        type_dict = {
            "CHR": hl.tstr,
            "POS": hl.tint32,
            "MarkerID": hl.tstr,
            "Allele1": hl.tstr,
            "Allele2": hl.tstr,
            "AC_Allele2": hl.tfloat64,
            "AF_Allele2": hl.tfloat64,
            "MissingRate": hl.tfloat64,
            "BETA": hl.tfloat64,
            "SE": hl.tfloat64,
            "Tstat": hl.tfloat64,
            "var": hl.tfloat64,
            "p.value": hl.tfloat64,
            "p.value.NA": hl.tfloat64,
            "Is.SPA": hl.tbool,
            "N": hl.tint32
        }

        if is_gxe:
            type_dict.update({
                "Beta_ge": hl.tstr,
                "seBeta_ge": hl.tstr,
                "pval_ge": hl.tstr,
                "pval_noSPA_ge": hl.tstr,
                "pval_ge_SKAT": hl.tstr,
            })

        ht = hl.import_table(paths, source_file_field="source", types=type_dict)

        if is_gxe:
            annot = {"pval_ge_SKAT": hl.float64(ht.pval_ge_SKAT)}
            for col in dynamic_cols:
                split_arr = ht[col].split(",")
                annot.update({f"{col}_{covar}": hl.float64(split_arr[i]) for i, covar in enumerate(dynamic_covariates)})
            ht = ht.annotate(**annot).drop(*dynamic_cols)

        ht = ht.annotate(**hl.parse_variant(ht.MarkerID.replace("_", ":"), reference_genome="GRCh38"))
        ht = ht.transmute(phenotype_id=ht.source.split("step2_")[-1], cell_type="~{cell_type}")
        ht = ht.key_by("locus", "alleles", "phenotype_id", "cell_type")
        ht = ht.checkpoint(new_temp_file())
        ht = ht.annotate(cis=ht_interval.index(ht.locus, all_matches=True).any(lambda x: x.gene_id == ht.phenotype_id))
        ht = ht.checkpoint(new_temp_file())
        ht2 = ht.drop("CHR", "POS", "MarkerID", "Allele1", "Allele2", "MissingRate").checkpoint("~{tmp_bucket}/~{prefix}.ht", overwrite=True)
        ht2.export("~{tmp_bucket}/~{prefix}.txt.bgz", parallel="header_per_shard")

        # export cis
        ht = ht.rename({'CHR' : '#CHR'}).key_by()
        cis_cols = ["#CHR", "POS", "cell_type", "phenotype_id", "MarkerID", "AF_Allele2", "BETA", "SE", "p.value"]
        if is_gxe:
            cis_cols += [f"{col}_{covar}" for covar in dynamic_covariates for col in dynamic_cols] + ["pval_ge_SKAT"]
        ht.filter(ht.cis).select(*cis_cols).export("~{tmp_bucket}/~{prefix}.SAIGE.cis_qtl_pairs.txt.bgz")

        # save gene_ids
        df = pd.DataFrame({"gene_ids": "~{sep=',' gene_ids}".split(",")})
        hl.Table.from_pandas(df).export("~{tmp_bucket}/~{prefix}.gene_ids.txt", header=False)
        __EOF__

        PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=$mem pyspark-shell" \
        python3 script.py && \
        touch _SUCCESS
    >>>

    output {
        String out_ht = "~{tmp_bucket}/~{prefix}.ht"
        String out_cis = "~{tmp_bucket}/~{prefix}.SAIGE.cis_qtl_pairs.txt.bgz"
        File success = "_SUCCESS"
    }

    runtime {
        docker: docker
        cpu: 32
        memory: "~{ceil(25 * size(result_paths, 'GB')) + 10} GB"
        disks: "local-disk 50 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task combine_cis {
    input {
        String out_dir
        Array[String] step2_ht_paths
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e

        mem=$(grep MemTotal /proc/meminfo | awk '{printf("%dg\n", $2 / 1024 / 1024)}')
        echo "$mem"

        cat << "__EOF__" > script.py
        import hail as hl
        import pandas as pd
        from hail.utils import new_temp_file

        # export step2_ht_paths
        step2_ht_paths = [x for x in "~{sep=',' step2_ht_paths}".split(",") if x != ""]
        df = pd.DataFrame({"step2_ht_paths": step2_ht_paths})
        df.to_csv("step2_ht_paths.txt", sep="\t", index=False, header=False)

        # combine cis files
        cis_files = [x.replace(".ht", ".SAIGE.cis_qtl_pairs.txt.bgz") for x in step2_ht_paths]
        ht = hl.import_table(cis_files, impute=True)
        ht = ht.annotate(CHR_NUM=hl.int(ht["#CHR"].replace("chr", "").replace("X", "23")))
        ht = ht.checkpoint(new_temp_file())
        ht = ht.order_by(ht.CHR_NUM, ht.POS, ht.phenotype_id)
        ht.drop("CHR_NUM").export("~{out_dir}/step2_cis/~{prefix}.SAIGE.cis_qtl_pairs.txt.bgz")
        __EOF__

        PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=$mem pyspark-shell" \
        python3 script.py && \
        touch _SUCCESS
    >>>

    output {
        File out_step2_ht_list = "step2_ht_paths.txt"
        String out_cis = out_dir + "/step2_cis/" + prefix + ".SAIGE.cis_qtl_pairs.txt.bgz"
        File success = "_SUCCESS"
    }

    runtime {
        docker: docker
        cpu: 32
        memory: "256 GB"
        disks: "local-disk 100 HDD"
        zones: zones
        preemptible: 2
    }
}

task step3_gene_pvalue {
    input {
        File cis
        Array[String] group_by = ["phenotype_id", "cell_type"]
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e

        mv ~{cis} ~{prefix}.gz

        cat << "__EOF__" > script.R
        library(dplyr)

        vars_to_group_by = unlist(stringr::str_split("~{sep=',' group_by}", ","))

        # run ACAT
        df = data.table::fread("~{prefix}.gz", data.table=FALSE)
        df.out =
            dplyr::group_by(df, dplyr::across(tidyselect::all_of(vars_to_group_by))) %>%
            dplyr::summarize(
                ACAT_p = SAIGEQTL::get_CCT_pvalue(p.value),
                top_pval = min(p.value),
                top_MarkerID = MarkerID[which.min(p.value)]
            ) %>%
            dplyr::ungroup()

        data.table::fwrite(df.out, "~{prefix}.SAIGE.acat_p.txt", quote = F, row.names = F, sep = "\t")
        __EOF__

        Rscript script.R && \
        touch _SUCCESS
    >>>

    output {
        File out_acat_p = prefix + ".SAIGE.acat_p.txt"
        File success = "_SUCCESS"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{ceil(16 * size(cis, 'GB')) + 10} GB"
        disks: "local-disk 50 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task step3_gxe_pvalue {
    input {
        File cis
        Array[String] group_by = ["phenotype_id", "cell_type"]
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e

        mv ~{cis} ~{prefix}.gz

        cat << "__EOF__" > script.R
        library(dplyr)
        library(tidyr)
        library(SAIGEQTL)

        df <- data.table::fread("~{prefix}.gz", data.table = FALSE)

        vars_to_group_by <- unlist(strsplit("~{sep=',' group_by}", ","))

        # Auto-detect p-value columns: p.value + pval_ge_* + pval_ge_SKAT
        pval_cols <- c("p.value", grep("^pval_ge_", colnames(df), value = TRUE))
        pval_cols <- unique(pval_cols)
        cat("Processing p-value columns:", paste(pval_cols, collapse = ", "), "\n")

        # Pivot to long format: one row per gene x pval_column x variant
        df_out <-
            tidyr::pivot_longer(
                df,
                cols = tidyselect::all_of(pval_cols),
                names_to = "pval_column",
                values_to = "pval"
            ) %>%
            dplyr::mutate(pval = as.numeric(pval)) %>%
            dplyr::group_by(dplyr::across(tidyselect::all_of(c(vars_to_group_by, "pval_column")))) %>%
            dplyr::filter(!is.na(pval)) %>%
            dplyr::summarize(
                ACAT_p = tryCatch(SAIGEQTL::get_CCT_pvalue(pval), error = function(e) NA_real_),
                top_MarkerID = MarkerID[which.min(pval)],
                top_pval = min(pval),
                .groups = "drop"
            ) %>%
            dplyr::group_by(pval_column) %>%
            dplyr::mutate(ACAT_q = qvalue::qvalue(ACAT_p)$qvalue) %>%
            dplyr::ungroup()

        data.table::fwrite(df_out, "~{prefix}.SAIGE.acat_p.txt", quote = FALSE, sep = "\t")
        __EOF__

        Rscript script.R && \
        touch _SUCCESS
    >>>

    output {
        File out_acat_p = prefix + ".SAIGE.acat_p.txt"
        File success = "_SUCCESS"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{ceil(16 * size(cis, 'GB')) + 10} GB"
        disks: "local-disk 50 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task postprocess {
    input {
        String out_dir
        File cis
        File acat_p
        File step2_ht_list
        Array[File] skip_gene_ids
        String prefix
        Boolean is_gxe = false
        String docker
        String zones
    }

    command <<<
        set -e

        n_cpu=$(grep -c ^processor /proc/cpuinfo)

        cat ~{sep=' ' skip_gene_ids} | sort > ~{prefix}.skip_gene_ids.txt
        tabix -@${n_cpu} -s 1 -b 2 -e 2 ~{cis}

        if ~{true='true' false='false' is_gxe}
        then
            cat ~{acat_p} | bgzip -c > ~{prefix}.SAIGE.acat.txt.gz
        else
            cat << "__EOF__" > qvalue.R
            library(dplyr)

            df = data.table::fread("~{acat_p}", data.table=FALSE) %>%
                dplyr::mutate(ACAT_q = qvalue::qvalue(ACAT_p)$qvalue)

            write.table(df, "~{prefix}.SAIGE.acat.txt", quote = F, row.names = F, sep = "\t")
            __EOF__

            Rscript qvalue.R && \
            bgzip ~{prefix}.SAIGE.acat.txt
        fi

        gcloud storage cp ~{prefix}.skip_gene_ids.txt ~{out_dir}/step2/~{prefix}.skip_gene_ids.txt && \
        gcloud storage cp ~{step2_ht_list} ~{out_dir}/step2/~{prefix}.step2_ht_paths.txt && \
        gcloud storage cp ~{cis}.tbi ~{out_dir}/step2_cis/ && \
        gcloud storage cp ~{prefix}.SAIGE.acat.txt.gz ~{out_dir}/step3/ && \
        touch _SUCCESS
    >>>

    output {
        Array[String] out_ht_paths = read_lines(step2_ht_list)
        File out_acat_q = "~{prefix}.SAIGE.acat.txt.gz"
        File success = "_SUCCESS"
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
