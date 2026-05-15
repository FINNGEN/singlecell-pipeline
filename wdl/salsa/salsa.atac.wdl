version 1.0

import "salsa.tasks.wdl" as tasks

workflow salsa {
    input {
        File obs
        File reference
        File chrom_info
        String vcf_pattern
        String bam_pattern
        String fragments_pattern
        String counts_pattern
        String reference_dir = "refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
        Boolean debug = false
        Boolean overwrite = false
        String docker_salsa
        String docker_snapatac2
        String docker_suite
        String zones
    }

    call tasks.prepare_params {
        input:
            obs = obs,
            vcf_pattern = vcf_pattern,
            bam_pattern = bam_pattern,
            fragments_pattern = fragments_pattern,
            counts_pattern = counts_pattern,
            overwrite = overwrite,
            docker = docker_suite,
            zones = zones
    }

    Array[File] param_files = if debug then [prepare_params.out_param_files[0]] else prepare_params.out_param_files

    scatter (param_file in param_files) {
        Array[Array[String]] params = if debug then [read_tsv(param_file)[29]] else read_tsv(param_file)

        scatter (param in params) {
            String fgid2 = param[0]
            File vcf = param[1]
            File bam = param[2]

            call tasks.salsa_wasp {
                input:
                    fgid = fgid2,
                    vcf = vcf,
                    bam = bam,
                    reference = reference,
                    chrom_info = chrom_info,
                    modality = "atac",
                    reference_dir = reference_dir,
                    docker = docker_salsa,
                    zones = zones
            }

            call tasks.make_fragment_file {
                input:
                    bam = salsa_wasp.out_bam,
                    output_directory = sub(fragments_pattern, "/[^/]+$", ""),
                    docker = docker_snapatac2,
                    zones = zones
            }

            scatter (chunked_barcodes in salsa_wasp.out_chunked_barcodes) {
                call tasks.chunk_wasp_bam {
                    input:
                        bam = salsa_wasp.out_bam,
                        barcodes = chunked_barcodes,
                        docker = docker_suite,
                        zones = zones
                }

                call tasks.salsa_count {
                    input:
                        fgid = fgid2,
                        vcf = vcf,
                        bam = chunk_wasp_bam.out_bam,
                        bai = chunk_wasp_bam.out_bai,
                        reference = reference,
                        modality = "atac",
                        reference_dir = reference_dir,
                        docker = docker_salsa,
                        zones = zones
                }
            }

            call tasks.merge_counts {
                input:
                    fgid = fgid2,
                    single_cell_tables = salsa_count.out_single_cell_table,
                    single_cell_phased_tables = salsa_count.out_single_cell_phased_table,
                    output_directory = sub(counts_pattern, "/[^/]+$", ""),
                    docker = docker_suite,
                    zones = zones
            }
        }
    }
}
