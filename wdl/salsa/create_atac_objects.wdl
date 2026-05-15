version 1.0

workflow create_atac_objects {
    input {
        File obs
        File pool_list
        File annotation_rds
        File peaks
        String bam_pattern
        String metadata_pattern
        String fragments_pattern
        String h5ad_pattern
        String merged_object_prefix
        Boolean debug = false
        Boolean overwrite = false
        Int n_count_peaks_lower_bound = 2000
        Int n_count_peaks_upper_bound = 100000
        Float blacklist_ratio_threshold = 0.01
        Float nucleosome_signal_threshold = 2
        Float tss_enrichment_threshold = 4
        Float pct_reads_in_peaks_threshold = 40
        Array[String] chromosomes
        String docker_signac
        String docker_suite
        String zones
    }

    call prepare_params {
        input:
            obs = obs,
            pool_list = pool_list,
            bam_pattern = bam_pattern,
            metadata_pattern = metadata_pattern,
            fragments_pattern = fragments_pattern,
            h5ad_pattern = h5ad_pattern,
            chromosomes = chromosomes,
            overwrite = overwrite,
            docker = docker_suite,
            zones = zones
    }

    if (prepare_params.flag_save) {
        Array[File] param_files = if debug then [prepare_params.out_param_files[0]] else prepare_params.out_param_files

        scatter (param_file in param_files) {
            Array[Array[String]] params = if debug then [read_tsv(param_file)[29]] else read_tsv(param_file)

            scatter (param in params) {
                String fgid2 = param[0]
                File bam = param[1]
                File metadata = param[2]
                File fragments = param[3]
                File fragments_tbi = "~{fragments}.tbi"
                File h5ad = param[4]

                call save_h5ad {
                    input:
                        fgid = fgid2,
                        obs = obs,
                        annotation_rds = annotation_rds,
                        metadata = metadata,
                        fragments = fragments,
                        fragments_tbi = fragments_tbi,
                        peaks = peaks,
                        chromosomes = chromosomes,
                        n_count_peaks_lower_bound = n_count_peaks_lower_bound,
                        n_count_peaks_upper_bound = n_count_peaks_upper_bound,
                        blacklist_ratio_threshold = blacklist_ratio_threshold,
                        nucleosome_signal_threshold = nucleosome_signal_threshold,
                        tss_enrichment_threshold = tss_enrichment_threshold,
                        pct_reads_in_peaks_threshold = pct_reads_in_peaks_threshold,
                        docker = docker_signac,
                        zones = zones
                }

                call copy {
                    input:
                        output_dir = sub(h5ad_pattern, "/[^/]+$", ""),
                        h5ad = flatten([[save_h5ad.out_h5ad], save_h5ad.out_chrom_h5ad]),
                        docker = docker_suite,
                        zones = zones
                }
            }
        }
    }

    if (prepare_params.flag_merge) {
        scatter (chrom in chromosomes) {
            call prepare_h5ad_list {
                input:
                    prefix = basename(merged_object_prefix),
                    chrom = chrom,
                    h5ad = if prepare_params.flag_save then flatten(select_all([prepare_params.out_chrom_h5ad, save_h5ad.out_chrom_h5ad])) else prepare_params.out_chrom_h5ad,
                    docker = docker_suite,
                    zones = zones
            }

            scatter (h5ad_list in prepare_h5ad_list.out_chunk_lists) {
                call merge_h5ad as merge_h5ad_chunk {
                    input:
                        prefix = basename(merged_object_prefix),
                        chrom = chrom,
                        h5ad = read_lines(h5ad_list),
                        min_cells = 0,
                        docker = docker_suite,
                        zones = zones
                }
            }

            call merge_h5ad {
                input:
                    prefix = basename(merged_object_prefix),
                    chrom = chrom,
                    h5ad = merge_h5ad_chunk.out_h5ad,
                    docker = docker_suite,
                    zones = zones
            }
        }

        call copy as copy_merged {
            input:
                output_dir = sub(merged_object_prefix, "/[^/]+$", ""),
                h5ad = merge_h5ad.out_h5ad,
                docker = docker_suite,
                zones = zones
        }
    }
}

task prepare_params {
    input {
        File obs
        File pool_list
        String bam_pattern
        String metadata_pattern
        String fragments_pattern
        String h5ad_pattern
        Array[String] chromosomes
        Boolean overwrite
        Int chunk_size = 200
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        # for each pool_ffid
        python3 << "__EOF__"
        import math
        import numpy as np
        import pandas as pd
        import re
        from google.cloud import storage

        df = pd.read_csv("~{obs}", sep="\t")
        df["FINNGENID2"] = df["barcode"].str.split("-", n=2).str[0]
        df = df[["FINNGENID", "FINNGENID2", "ffid_2", "pool_name"]].drop_duplicates()

        # filter to QCed pools
        pools = pd.read_csv("~{pool_list}", header=None).iloc[:, 0]
        df = df.loc[df["pool_name"].isin(pools), :]
        # format paths
        df["bam"] = df["ffid_2"].apply(lambda x: re.sub(r"\{FFID\}", x, "~{bam_pattern}"))
        df["metadata"] = df["pool_name"].apply(lambda x: re.sub(r"\{POOL\}", x, "~{metadata_pattern}"))
        df["fragments"] = df["FINNGENID2"].apply(lambda x: re.sub(r"\{FINNGENID\}", x, "~{fragments_pattern}"))
        df["h5ad"] = df["FINNGENID2"].apply(lambda x: re.sub(r"\{FINNGENID\}", x, "~{h5ad_pattern}".replace("{CHROM}", "chr22")))
        df = df[['FINNGENID2', 'bam', "metadata", "fragments", "h5ad"]]

        # filter out non-existent bams
        storage_client = storage.Client()
        bucket = storage_client.get_bucket("~{bam_pattern}".split("/")[2])
        bam_exists = np.array([bucket.blob(x.split("/", 3)[-1]).exists() for x in df.bam])
        print(f"{np.sum(bam_exists)}/{len(bam_exists)} bams exist")

        # filter out finished samples
        bucket = storage_client.get_bucket("~{fragments_pattern}".split("/")[2])
        fragments_exists = np.array([bucket.blob(x.split("/", 3)[-1]).exists() for x in df.fragments])
        print(f"{np.sum(fragments_exists)}/{len(fragments_exists)} fragments exist")

        bucket = storage_client.get_bucket("~{h5ad_pattern}".split("/")[2])
        h5ad_exists = np.array([bucket.blob(x.split("/", 3)[-1]).exists() for x in df.h5ad])
        print(f"{np.sum(h5ad_exists)}/{len(h5ad_exists)} h5ad exist")

        overwrite = ~{true='True' false='False' overwrite}
        df_to_merge = df.loc[(not overwrite) & h5ad_exists, :]
        if len(df_to_merge.index) > 0:
            chromosomes = "~{sep=',' chromosomes}".split(",")
            pd.concat([df_to_merge.h5ad.str.replace(".chr22.", f".{chrom}.", regex=False) for chrom in chromosomes], axis=0).to_csv(
                "h5ad_to_merge.txt", sep="\t", index=False, header=False
            )
            print(f"Found {len(df_to_merge.index)} h5ad to merge")
        else:
            pd.DataFrame().to_csv("h5ad_to_merge.txt")

        df = df.loc[bam_exists & fragments_exists & (overwrite | ~h5ad_exists), :]
        df.to_csv("params.tsv", sep="\t", index=False, header=False)
        print(f"Processing {len(df.index)} samples for h5ad")

        with open("flag_merge.txt", "w") as f:
            f.write("true" if np.sum(bam_exists) == np.sum(fragments_exists) else "false")
        with open("flag_save.txt", "w") as f:
            f.write("true" if len(df.index) > 0 else "false")

        if len(df.index) > 0:
            # chunk params
            chunk_size= ~{chunk_size}
            total_lines = len(df.index)
            n_chunks = math.ceil(total_lines / chunk_size)
            lines_per_chunk = math.ceil(total_lines / n_chunks)

            for i in range(n_chunks):
                start_line = i * lines_per_chunk
                end_line = min(start_line + lines_per_chunk, total_lines)
                df.iloc[start_line:end_line].to_csv(f"params.{i}.tsv", sep="\t", index=False, header=False)
        else:
            pd.DataFrame().to_csv("params.0.tsv")
        __EOF__
    >>>

    output {
        File out_params = "params.tsv"
        Array[File] out_param_files = glob("params.*.tsv")
        Array[String] out_chrom_h5ad = read_lines("h5ad_to_merge.txt")
        Boolean flag_merge = read_boolean("flag_merge.txt")
        Boolean flag_save = read_boolean("flag_save.txt")
    }

    meta {
        volatile: true
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

task save_h5ad {
    input {
        String fgid
        File obs
        File annotation_rds
        File metadata
        File fragments
        File fragments_tbi
        File peaks
        Array[String] chromosomes
        Int n_count_peaks_lower_bound
        Int n_count_peaks_upper_bound
        Float blacklist_ratio_threshold
        Float nucleosome_signal_threshold
        Float tss_enrichment_threshold
        Float pct_reads_in_peaks_threshold
        Int min_cells = 0
        Int min_features = 200
        String docker
        String zones
    }

    # Int disk_space = ceil(extra_disk_space + 12 * size(fragments, "GB"))

    command <<<
        set -e

        cat << "__EOF__" > save_h5ad.R
        library(dplyr)
        library(GenomicRanges)

        metadata_file_path <- "~{metadata}"
        fragments_file_path <- "~{fragments}"
        peak_file_path <- "peaks.bed"

        obs <- data.table::fread("~{obs}", data.table = FALSE) %>%
            dplyr::mutate(
                ffid_barcode = dplyr::if_else(
                    stringr::str_detect(pool_name, "^pool"),
                    stringr::str_c(pool_name, "_", ffid, "-", stringr::str_remove(stringr::str_split_fixed(barcode, "-", 2)[, 2], "-1$")),
                    stringr::str_c(ffid_2, "-", stringr::str_remove(stringr::str_split_fixed(barcode, "-", 2)[, 2], "-1$"))
                )
            ) %>%
            dplyr::select(ffid_barcode, ffid_2, pool_name, predicted.celltype.l1, predicted.celltype.l2)
        metadata <- data.table::fread(metadata_file_path, data.table = FALSE) %>%
            dplyr::filter(
                nCount_peaks > ~{n_count_peaks_lower_bound} &
                nCount_peaks < ~{n_count_peaks_upper_bound} &
                blacklist_ratio < ~{blacklist_ratio_threshold} &
                nucleosome_signal < ~{nucleosome_signal_threshold} &
                TSS.enrichment > ~{tss_enrichment_threshold}
            ) %>%
            dplyr::inner_join(obs, by = "ffid_barcode") %>%
            dplyr::rename(ffid = ffid_2) %>%
            tibble::column_to_rownames("barcode")

        peaks <- rtracklayer::import(peak_file_path, format = 'BED')
        frags <- Signac::CreateFragmentObject(path = fragments_file_path)
        filtered_counts <- Signac::FeatureMatrix(fragments = frags, features = peaks, cells = rownames(metadata))

        # Print the dimensions of the filtered counts matrix
        cat("Dimensions of filtered counts:", dim(filtered_counts), "\n")

        chrom_assay <- Signac::CreateChromatinAssay(
            counts = filtered_counts,
            sep = c(":", "-"),
            fragments = fragments_file_path,
            min.cells = ~{min_cells},
            min.features = ~{min_features}
        )

        pbmc <- Seurat::CreateSeuratObject(
            counts = chrom_assay,
            assay = "peaks",
            meta.data = metadata
        )

        # Search for the Ensembl 98 EnsDb for Homo sapiens on AnnotationHub
        # ah <- AnnotationHub::AnnotationHub()
        # query(ah, "EnsDb.Hsapiens.v98")
        # ensdb_v98 <- ah[["AH75011"]]

        # extract gene annotations from EnsDb
        # annotations <- Signac::GetGRangesFromEnsDb(ensdb = ensdb_v98)

        # load annotations from RDS
        annotations <- readRDS("~{annotation_rds}")

        # change to UCSC style since the data was mapped to hg38
        seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
        genome(annotations) <- "hg38"

        # add the gene information to the object
        Signac::Annotation(pbmc) <- annotations

        # compute nucleosome signal score per cell
        pbmc <- Signac::NucleosomeSignal(object = pbmc)

        # compute TSS enrichment score per cell
        pbmc <- Signac::TSSEnrichment(object = pbmc)

        # add fraction of reads in peaks
        pbmc$pct_reads_in_peaks <- pbmc$atac_peak_region_fragments / pbmc$atac_fragments * 100

        # add blacklist ratio
        pbmc$blacklist_ratio <- Signac::FractionCountsInRegion(
            object = pbmc,
            assay = "peaks",
            regions = Signac::blacklist_hg38_unified
        )

        # apply QC again
        pbmc <- subset(
            x = pbmc,
            subset = nCount_peaks > ~{n_count_peaks_lower_bound} &
                nCount_peaks < ~{n_count_peaks_upper_bound} &
                blacklist_ratio < ~{blacklist_ratio_threshold} &
                nucleosome_signal < ~{nucleosome_signal_threshold} &
                TSS.enrichment > ~{tss_enrichment_threshold} &
                pct_reads_in_peaks > ~{pct_reads_in_peaks_threshold}
        )
        SeuratDisk::SaveH5Seurat(pbmc, filename = "~{fgid}.post_qc.h5seurat")
        SeuratDisk::Convert("~{fgid}.post_qc.h5seurat", dest = "h5ad")
        __EOF__

        cat << "__EOF__" > chunk.py
        import anndata as ad
        import numpy as np
        import scanpy as sc
        import scipy.sparse as sp

        adata = ad.read_h5ad("~{fgid}.post_qc.tmp.h5ad")
        adata.obs = adata.obs.rename(columns={"ffid_barcode": "barcode"}).set_index('barcode')
        adata.obs = adata.obs.drop(['orig.ident'], axis=1)
        adata.raw = None
        adata.X = sp.csc_matrix(adata.X)
        adata.X = adata.X.astype(np.uint32)
        adata.write("~{fgid}.post_qc.h5ad")

        chromosomes = "~{sep=',' chromosomes}".split(",")

        for chrom in chromosomes:
            idx = adata.var_names.str.startswith(f"{chrom}-")
            adata[:, idx].write(f"~{fgid}.post_qc.{chrom}.h5ad")
        __EOF__

        touch ~{fragments_tbi} && \
        zcat ~{peaks} | cut -f1-3 > peaks.bed && \
        Rscript save_h5ad.R && \
        mv ~{fgid}.post_qc.h5ad ~{fgid}.post_qc.tmp.h5ad && \
        python3 chunk.py && \
        rm ~{fgid}.post_qc.tmp.h5ad && \
        touch _SUCCESS

    >>>

    output {
        File out_success = "_SUCCESS"
        File out_h5ad = "~{fgid}.post_qc.h5ad"
        Array[File] out_chrom_h5ad = glob("~{fgid}.post_qc.*.h5ad")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{ceil(12 * size(fragments, 'GB')) + 10} GB"
        disks: "local-disk ~{ceil(12 * size(fragments, 'GB')) + 10} HDD"
        zones: zones
        preemptible: 2
    }
}

task prepare_h5ad_list {
    input {
        String prefix
        String chrom
        Array[String] h5ad
        Int chunk_size = 200
        String docker
        String zones
    }

    command <<<
        set -e

        echo "~{sep=' ' h5ad}" | tr " " "\n" | grep ".post_qc.~{chrom}.h5ad" > ~{prefix}.~{chrom}.h5ad_list.txt

        python3 << "__EOF__"
        import math
        import pandas as pd

        df = pd.read_csv("~{prefix}.~{chrom}.h5ad_list.txt", header=None)

        chunk_size = ~{chunk_size}
        total_lines = len(df.index)
        n_chunks = math.ceil(total_lines / chunk_size)
        lines_per_chunk = math.ceil(total_lines / n_chunks)

        for i in range(n_chunks):
            start_line = i * lines_per_chunk
            end_line = min(start_line + lines_per_chunk, total_lines)
            df.iloc[start_line:end_line].to_csv(f"h5ad_list.{i}.txt", sep="\t", index=False, header=False)
        __EOF__

    >>>

    output {
        File out_list = "~{prefix}.~{chrom}.h5ad_list.txt"
        Array[File] out_chunk_lists = glob("h5ad_list.*.txt")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}

task merge_h5ad {
    input {
        String prefix
        String chrom
        Array[File] h5ad
        Int min_cells = 10
        String docker
        String zones
    }

    command <<<
        set -e

        cat << "__EOF__" > script.py
        import anndata as ad
        import os.path
        import scanpy as sc

        h5ad_files = "~{sep=',' h5ad}".split(",")

        res = []
        for h5ad in h5ad_files:
            adata = ad.read_h5ad(h5ad)
            res.append(adata)
            print(f"{len(res)} / {len(h5ad_files)}")

        combined_atac = ad.concat(res, axis = 0, join="outer")
        combined_atac.raw = None

        sc.pp.filter_genes(combined_atac, min_cells=~{min_cells})
        combined_atac.write("~{prefix}.~{chrom}.h5ad")
        __EOF__

        python3 script.py && \
        touch _SUCCESS

    >>>

    output {
        File out_success = "_SUCCESS"
        File out_h5ad = "~{prefix}.~{chrom}.h5ad"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{ceil(10 * size(h5ad, 'GB')) + 10} GB"
        disks: "local-disk ~{ceil(3 * size(h5ad, 'GB')) + 10} HDD"
        zones: zones
        preemptible: 2
    }
}

task copy {
    input {
        String output_dir
        Array[String] h5ad
        String docker
        String zones
    }

    command {
        set -e
        gcloud storage cp ~{sep=" " h5ad} ~{output_dir} && \
        touch _SUCCESS
    }

    output {
        File out_success = "_SUCCESS"
        String result = output_dir
    }

    meta {
        volatile: true
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}
