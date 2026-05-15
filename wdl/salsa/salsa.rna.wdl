version 1.0

import "salsa.tasks.wdl" as tasks

workflow salsa {
    input {
        File obs
        File reference
        File chrom_info
        String vcf_pattern
        String bam_pattern
        String counts_pattern
        String reference_dir = "refdata-gex-GRCh38-2020-A"
        Boolean debug = false
        Boolean overwrite = false
        String docker_salsa
        String docker_suite
        String zones
        # STAR needs the genome resident in memory (~30 GB for GRCh38).
        Int salsa_wasp_memory_gb_floor = 40
    }

    call tasks.prepare_params {
        input:
            obs = obs,
            vcf_pattern = vcf_pattern,
            bam_pattern = bam_pattern,
            counts_pattern = counts_pattern,
            overwrite = overwrite,
            docker = docker_suite,
            zones = zones
    }

    Array[File] param_files = if debug then [prepare_params.out_param_files[0]] else prepare_params.out_param_files

    scatter (param_file in param_files) {
        Array[Array[String]] params = if debug then [read_tsv(param_file)[0]] else read_tsv(param_file)

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
                    modality = "rna",
                    reference_dir = reference_dir,
                    memory_gb_floor = salsa_wasp_memory_gb_floor,
                    docker = docker_salsa,
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
                        modality = "rna",
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
