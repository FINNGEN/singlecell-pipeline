version 1.0

import "preprocess_h5ad.tasks.wdl" as tasks

workflow preprocess_h5ad {
    input {
        String type
        Array[String] chromosomes
        String h5ad_pattern
        File gex_features
        File tss_bed
        Array[File]? obs_annot
        Array[String]? obs_key
        Array[String]? obs_annot_key
        String obs_sample_id
        String? obs_pool_id
        Boolean update_barcodes
        String? resequenced_suffix
        Float? downsample
        Boolean calculate_qc_metrics
        Array[String] vars_to_regress
        Array[String] pseudobulk_cell_type_annots
        String pseudobulk_all_cell_annot
        Array[String] pseudobulk_agg_methods
        Int? pseudobulk_min_cells
        Float? pseudobulk_min_sample_prop
        Float? pseudobulk_min_cell_prop
        Boolean chunk_sctransform
        Boolean chunk_by_genes
        Int split_num_genes
        String docker_archr
        String docker_suite
        String docker_seurat
        String zones
    }

    scatter (chrom in chromosomes) {
        call tasks.update_h5ad {
            input:
                h5ad = sub(h5ad_pattern, "\\{CHROM\\}", chrom),
                type = type,
                features = gex_features,
                obs_annot = obs_annot,
                obs_key = obs_key,
                obs_annot_key = obs_annot_key,
                obs_sample_id = obs_sample_id,
                obs_pool_id = obs_pool_id,
                update_barcodes = update_barcodes,
                resequenced_suffix = resequenced_suffix,
                vars_to_regress = vars_to_regress,
                downsample = downsample,
                calculate_qc_metrics = calculate_qc_metrics,
                docker = docker_suite,
                zones = zones
        }

        call tasks.generate_peak_bed {
            input:
                var = update_h5ad.out_var,
                docker = docker_suite,
                zones = zones
        }

        scatter (cell_type_annot in pseudobulk_cell_type_annots) {
            scatter (agg_method in pseudobulk_agg_methods) {
                call tasks.pseudobulk {
                    input:
                        h5ad = update_h5ad.out_h5ad,
                        tss_bed = generate_peak_bed.out_bed,
                        sample_id = obs_sample_id,
                        agg_method = agg_method,
                        cell_type_annot = cell_type_annot,
                        all_cell_annot = pseudobulk_all_cell_annot,
                        min_cells = pseudobulk_min_cells,
                        min_sample_prop = pseudobulk_min_sample_prop,
                        min_cell_prop = pseudobulk_min_cell_prop,
                        normalize_by_peak_length = true,
                        docker = docker_suite,
                        zones = zones
                }
            }
        }

        if (chunk_by_genes) {
            call tasks.split_ensgene {
                input:
                    gex_features = gex_features,
                    chrom = chrom,
                    num_genes = split_num_genes,
                    docker = docker_suite,
                    zones = zones
            }

            scatter (ensgene_list in split_ensgene.out_lists) {
                Array[String] gene_ids = read_lines(ensgene_list)

                call tasks.chunk_gene_id {
                    input:
                        h5ad = update_h5ad.out_h5ad,
                        type = type,
                        gene_ids = gene_ids,
                        tss_bed = tss_bed,
                        docker = docker_suite,
                        zones = zones
                }
            }
        }
    }

    call tasks.merge_atac_obs {
        input:
            obs = update_h5ad.out_obs,
            prefix = "obs",
            docker = docker_suite,
            zones = zones
    }

    scatter (h5ad in update_h5ad.out_h5ad) {
        call tasks.update_atac_h5ad {
            input:
                h5ad = h5ad,
                obs = merge_atac_obs.out_obs,
                docker = docker_suite,
                zones = zones
        }
    }

    # original: n_chroms x n_cell_types x n_agg_methods x n_annots
    # transposed: n_annots x n_chrom x n_agg_methods x n_cell_types
    Array[Array[Array[Array[File]]]] outer_inv_bed = transpose(pseudobulk.out_inv_bed)
    Array[Array[Array[Array[File]]]] outer_raw_bed = transpose(pseudobulk.out_raw_bed)
    Array[Array[Array[Array[File]]]] outer_nonzero_prop = transpose(pseudobulk.out_nonzero_prop)
    Array[Array[Array[Array[File]]]] outer_barcode = transpose(pseudobulk.out_barcode)

    scatter (i in range(length(pseudobulk_cell_type_annots))) {
        # original: n_chrom x n_agg_methods x n_cell_types
        # transposed: n_agg_methods x n_chrom x n_cell_types
        Array[Array[Array[File]]] inner_inv_bed = transpose(outer_inv_bed[i])
        Array[Array[Array[File]]] inner_raw_bed = transpose(outer_raw_bed[i])
        Array[Array[Array[File]]] inner_nonzero_prop = transpose(outer_nonzero_prop[i])
        Array[Array[Array[File]]] inner_barcode = transpose(outer_barcode[i])

        scatter (j in range(length(pseudobulk_agg_methods))) {
            # original: n_chrom x n_cell_types
            # transposed: n_cell_types x n_chrom
            Array[Array[File]] inv_bed = transpose(inner_inv_bed[j])
            Array[Array[File]] raw_bed = transpose(inner_raw_bed[j])
            Array[Array[File]] nonzero_prop = transpose(inner_nonzero_prop[j])
            Array[Array[File]] barcode = transpose(inner_barcode[j])

            scatter (k in range(length(inv_bed))) {
                call tasks.combine_pseudobulk {
                    input:
                        inv_bed = inv_bed[k],
                        raw_bed = raw_bed[k],
                        nonzero_prop = nonzero_prop[k],
                        barcode = [barcode[k][0]],
                        agg_method = pseudobulk_agg_methods[j],
                        docker = docker_suite,
                        zones = zones
                }

                call calculate_gc_content {
                    input:
                        inv_bed = combine_pseudobulk.out_inv_bed,
                        docker = docker_archr,
                        zones = zones
                }
            }
        }
    }
}

task calculate_gc_content {
    input {
        File inv_bed
        String prefix = basename(inv_bed, ".bed.gz")
        String docker
        String zones
    }

    Int disk_space = ceil(size(inv_bed, 'GB')) + 10

    command <<<
        set -e

        cat << "__EOF__" > script.R
        library(ArchR)
        library(BSgenome.Hsapiens.UCSC.hg38)
        library(dplyr)

        inv <- data.table::fread("~{inv_bed}", data.table = F)

        # compute GC content
        ArchR:::.requirePackage("Biostrings",source="bioc")
        ArchR::addArchRGenome("hg38")
        genomeAnnotation <- ArchR::getGenomeAnnotation()
        BSgenome <- eval(parse(text = genomeAnnotation$genome))
        BSgenome <- ArchR::validBSgenome(BSgenome)

        nucFreq <- BSgenome::alphabetFrequency(BSgenome::getSeq(BSgenome, inv$chr, start = inv$start, end = inv$end))
        inv <- dplyr::mutate(inv, percentage_gc_content = rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq))

        write.table(inv, "~{prefix}.gc_content.bed", sep = "\t", quote = F, row.names = F)
        __EOF__

        Rscript script.R && \
        bgzip ~{prefix}.gc_content.bed && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_bed = "~{prefix}.gc_content.bed.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "16 GB"
        disks: "local-disk ~{disk_space} HDD"
        zones: zones
        preemptible: 2
    }
}