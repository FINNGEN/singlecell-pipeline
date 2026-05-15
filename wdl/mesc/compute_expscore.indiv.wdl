version 1.0

import "mesc.tasks.wdl" as tasks

workflow compute_expscore {
    input {
        Array[String] chromosomes
        File cell_type_list
        String bed_pattern
        String covar_pattern
        String exp_bfile_pattern
        String geno_bfile_pattern
        Boolean debug = false
        String docker_suite
        String docker_ldsc
        String zones
    }

    Array[String] cell_types = if debug then [read_lines(cell_type_list)[0]] else read_lines(cell_type_list)

    scatter (cell_type in cell_types) {
        File bed = sub(bed_pattern, "\\{CELL_TYPE\\}", cell_type)
        File covar = sub(covar_pattern, "\\{CELL_TYPE\\}", cell_type)
        String prefix = "~{basename(bed, '.mean.inv.bed.gz')}.mesc"

        call tasks.preprocess_indiv_expr {
            input:
                bed = bed,
                covar = covar,
                prefix = prefix,
                docker = docker_suite,
                zones = zones
        }

        scatter (chrom in chromosomes) {
            String prefix_chrom = "~{prefix}.~{sub(chrom, '^chr', '')}"

            String exp_bfile = sub(exp_bfile_pattern, "\\{CHR\\}", chrom)
            File exp_bed = "~{exp_bfile}.bed"
            File exp_bim = "~{exp_bfile}.bim"
            File exp_fam = "~{exp_bfile}.fam"
            String geno_bfile = sub(geno_bfile_pattern, "\\{CHR\\}", chrom)
            File geno_bed = "~{geno_bfile}.bed"
            File geno_bim = "~{geno_bfile}.bim"
            File geno_fam = "~{geno_bfile}.fam"

            call tasks.preprocess_chunks {
                input:
                    chrom = chrom,
                    expression_matrix = preprocess_indiv_expr.out_expression_matrix,
                    docker = docker_suite,
                    zones = zones
            }

            scatter (chunk in preprocess_chunks.out_chunks) {
                String prefix_chunk = "~{prefix}.~{basename(chunk, '.txt')}"
                String prefix_chunk_chrom = "~{prefix_chunk}.~{sub(chrom, '^chr', '')}"

                call tasks.compute_expscore_indiv_chunk {
                    input:
                        chunk = chunk,
                        chrom = chrom,
                        exp_bed = exp_bed,
                        exp_bim = exp_bim,
                        exp_fam = exp_fam,
                        geno_bed = geno_bed,
                        geno_bim = geno_bim,
                        geno_fam = geno_fam,
                        expression_matrix = preprocess_indiv_expr.out_expression_matrix,
                        covar = preprocess_indiv_expr.out_covar,
                        prefix = prefix_chunk,
                        prefix_chrom = prefix_chunk_chrom,
                        docker = docker_ldsc,
                        zones = zones
                }
            }

            call tasks.compute_expscore_from_lasso {
                input:
                    chrom = chrom,
                    geno_bed = geno_bed,
                    geno_bim = geno_bim,
                    geno_fam = geno_fam,
                    lasso = compute_expscore_indiv_chunk.out_lasso,
                    hsq = compute_expscore_indiv_chunk.out_hsq,
                    prefix = prefix,
                    prefix_chrom = prefix_chrom,
                    docker = docker_ldsc,
                    zones = zones
            }
        }

        call tasks.tar_expscore {
            input:
                prefix = prefix,
                expscore_files = flatten([
                    compute_expscore_from_lasso.out_G,
                    compute_expscore_from_lasso.out_ave,
                    compute_expscore_from_lasso.out_expscore,
                    compute_expscore_from_lasso.out_gannot,
                    compute_expscore_from_lasso.out_hsq,
                    compute_expscore_from_lasso.out_lasso
                ]),
                hsq_files = compute_expscore_from_lasso.out_hsq,
                docker = docker_suite,
                zones = zones
        }
    }
}
