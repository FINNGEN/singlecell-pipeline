version 1.0

import "preprocess_h5ad.tasks.wdl" as tasks

workflow downsample_h5ad {
    input {
        Array[String] chromosomes
        File tss_bed
        File gex_h5ad
        String atac_h5ad_pattern
        File barcode_file_list
        Boolean debug = false
        String docker_suite
        String zones
    }

    Array[File] barcode_files = if debug then [read_lines(barcode_file_list)[0]] else read_lines(barcode_file_list)

    scatter (chrom in chromosomes) {
        call tasks.subset_h5ad as subset_atac_h5ad {
            input:
                h5ad = sub(atac_h5ad_pattern, "\\{CHR\\}", chrom),
                barcode_files = barcode_files,
                recompute_metrics = false,
                docker = docker_suite,
                zones = zones
        }
    }

    scatter (barcode_file in barcode_files) {
        call gather_atac_h5ad {
            input:
                h5ad = flatten(subset_atac_h5ad.out_h5ad),
                obs = flatten(subset_atac_h5ad.out_obs),
                condition = basename(barcode_file, ".txt.gz"),
                docker = docker_suite,
                zones = zones
        }

        call tasks.merge_atac_obs {
            input:
                obs = gather_atac_h5ad.out_obs,
                prefix = "obs",
                docker = docker_suite,
                zones = zones
        }

        scatter (h5ad in gather_atac_h5ad.out_h5ad) {
            call tasks.update_atac_h5ad {
                input:
                    h5ad = h5ad,
                    obs = merge_atac_obs.out_obs,
                    docker = docker_suite,
                    zones = zones
            }

            call tasks.generate_peak_bed {
                input:
                    var = update_atac_h5ad.out_var,
                    docker = docker_suite,
                    zones = zones
            }

            call tasks.pseudobulk as pseudobulk_atac {
                input:
                    h5ad = update_atac_h5ad.out_h5ad,
                    tss_bed = generate_peak_bed.out_bed,
                    sample_id = "FINNGENID",
                    agg_method = "sum",
                    cell_type_annot = "predicted.celltype.l1",
                    all_cell_annot = "CD4_CD8_T",
                    min_cells = 10,
                    min_sample_prop = 0.5,
                    min_cell_prop = 0,
                    normalize_by_peak_length = true,
                    docker = docker_suite,
                    zones = zones
            }
        }

        call filter_CD4_CD8_T {
            input:
                inv_bed = flatten(pseudobulk_atac.out_inv_bed),
                raw_bed = flatten(pseudobulk_atac.out_raw_bed),
                nonzero_prop = flatten(pseudobulk_atac.out_nonzero_prop),
                barcode = flatten(pseudobulk_atac.out_barcode),
                docker = docker_suite,
                zones = zones
        }

        call tasks.combine_pseudobulk as combine_pseudobulk_atac {
            input:
                inv_bed = filter_CD4_CD8_T.out_inv_bed,
                raw_bed = filter_CD4_CD8_T.out_raw_bed,
                nonzero_prop = filter_CD4_CD8_T.out_nonzero_prop,
                barcode = [filter_CD4_CD8_T.out_barcode[0]],
                agg_method = "sum",
                docker = docker_suite,
                zones = zones
        }
    }
}

task gather_atac_h5ad {
    input {
        Array[String] h5ad
        Array[String] obs
        String condition
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        echo "~{sep=',' h5ad}" | tr ',' '\n' | grep ~{condition} > h5ad_list.txt
        echo "~{sep=',' obs}" | tr ',' '\n' | grep ~{condition} > obs_list.txt
    >>>

    output {
        Array[File] out_h5ad = read_lines("h5ad_list.txt")
        Array[File] out_obs = read_lines("obs_list.txt")
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

task filter_CD4_CD8_T {
    input {
        Array[String] inv_bed
        Array[String] raw_bed
        Array[String] nonzero_prop
        Array[String] barcode
        String pattern = "predicted.celltype.l1.CD4_CD8_T"
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        echo "~{sep=',' inv_bed}" | tr ',' '\n' | grep ~{pattern} > inv_bed.txt
        echo "~{sep=',' raw_bed}" | tr ',' '\n' | grep ~{pattern} > raw_bed.txt
        echo "~{sep=',' nonzero_prop}" | tr ',' '\n' | grep ~{pattern} > nonzero_prop.txt
        echo "~{sep=',' barcode}" | tr ',' '\n' | grep ~{pattern} > barcode.txt
    >>>

    output {
        Array[File] out_inv_bed = read_lines("inv_bed.txt")
        Array[File] out_raw_bed = read_lines("raw_bed.txt")
        Array[File] out_nonzero_prop = read_lines("nonzero_prop.txt")
        Array[File] out_barcode = read_lines("barcode.txt")
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
