version 1.0

task preprocess_indiv_expr {
    input {
        File bed
        File covar
        String prefix
        String docker
        String zones
    }

    command <<<
        set -e

        cat << "__EOF__" > script.R
        library(dplyr)

        covar <- data.table::fread("~{covar}", data.table = FALSE) %>%
            dplyr::mutate(FID = 0) %>%
            dplyr::select(FID, IID, tidyselect::everything())
        expr_mat <- data.table::fread("~{bed}", data.table = FALSE)

        samples <- intersect(covar$IID, colnames(expr_mat)[5:ncol(expr_mat)])
        
        covar <- dplyr::filter(covar, IID %in% samples)
        expr_mat <-
            dplyr::filter(expr_mat, chr != "chrX") %>%
            dplyr::rename(GENE = gene_id, CHR = chr, GENE_COORD = end) %>%
            dplyr::mutate(CHR = stringr::str_replace(CHR, "^chr", "")) %>%
            dplyr::select(GENE, CHR, GENE_COORD, tidyselect::all_of(covar$IID))

        write.table(covar, file = "~{prefix}.covar.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
        write.table(expr_mat, file = "~{prefix}.expr_mat.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
        __EOF__

        Rscript script.R && \
        bgzip ~{prefix}.covar.tsv && \
        bgzip ~{prefix}.expr_mat.tsv && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_covar = "~{prefix}.covar.tsv.gz"
        File out_expression_matrix = "~{prefix}.expr_mat.tsv.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "8 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}


task preprocess_chunks {
    input {
        String chrom
        File expression_matrix
        Int chunk_size = 200
        Boolean debug = false
        String docker
        String zones
    }

    command <<<
        set -e

        python3 << "__EOF__"
        import math
        import pandas as pd
        import re

        chrom = "~{chrom}".replace("chr", "")
        df = pd.read_csv("~{expression_matrix}", sep="\t")
        df = df.loc[df.CHR.astype(str) == chrom, :]

        # chunk params
        chunk_size= ~{chunk_size}
        total_lines = len(df.index)
        n_chunks = math.ceil(total_lines / chunk_size)
        lines_per_chunk = math.ceil(total_lines / n_chunks)

        if ~{true='True' false='False' debug}:
            n_chunks = 1

        for i in range(n_chunks):
            start_line = i * lines_per_chunk
            end_line = min(start_line + lines_per_chunk, total_lines)
            df.GENE.iloc[start_line:end_line].to_csv(f"chunk.{i}.tsv", sep="\t", index=False, header=False)
        __EOF__
    >>>

    output {
        Array[File] out_chunks = glob("chunk.*.tsv")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "8 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}

task compute_expscore {
    input {
        File sumstats
        String prefix
        String prefix_chrom
        String docker
        String zones
    }

    command <<<
        set -e

        export PYENV_ROOT="$HOME/.pyenv"
        [[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"
        eval "$(pyenv init - bash)"

        python ~/mesc/run_mesc.py \
        --compute-expscore-sumstat \
        --eqtl-sumstat ~{sumstats} \
        --out ~{prefix}
    >>>

    output {
        File out_G = "~{prefix_chrom}.G"
        File out_ave = "~{prefix_chrom}.ave_h2cis"
        File out_expscore = "~{prefix_chrom}.expscore.gz"
        File out_gannot = "~{prefix_chrom}.gannot.gz"
        File out_hsq = "~{prefix_chrom}.hsq"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "8 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task compute_expscore_indiv {
    input {
        String chrom
        File exp_bed
        File exp_bim
        File exp_fam
        File geno_bed
        File geno_bim
        File geno_fam
        File expression_matrix
        File covar
        String prefix
        String prefix_chrom
        String docker
        String zones
    }

    command <<<
        set -e

        export PYENV_ROOT="$HOME/.pyenv"
        [[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"
        eval "$(pyenv init - bash)"

        mkdir -p exp_bfile geno_bfile
        ln -s ~{exp_bed} ~{exp_bim} ~{exp_fam} exp_bfile/
        ln -s ~{geno_bed} ~{geno_bim} ~{geno_fam} geno_bfile/
        exp_bfile=exp_bfile/$(basename ~{exp_bed} .bed)
        geno_bfile=geno_bfile/$(basename ~{geno_bed} .bed)

        gunzip -c ~{covar} > ~{basename(covar, ".gz")}

        python ~/mesc/run_mesc.py \
        --compute-expscore-indiv \
        --plink-path $(which plink) \
        --expression-matrix ~{expression_matrix} \
        --exp-bfile $exp_bfile \
        --geno-bfile $geno_bfile \
        --chr ~{sub(chrom, "^chr", "")} \
        --covariates ~{basename(covar, ".gz")} \
        --out ~{prefix}
    >>>

    output {
        File out_G = "~{prefix_chrom}.G"
        File out_ave = "~{prefix_chrom}.ave_h2cis"
        File out_expscore = "~{prefix_chrom}.expscore.gz"
        File out_gannot = "~{prefix_chrom}.gannot.gz"
        File out_hsq = "~{prefix_chrom}.hsq"
        File out_lasso = "~{prefix_chrom}.lasso"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "32 GB"
        disks: "local-disk 20 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task compute_expscore_indiv_chunk {
    input {
        File chunk
        String chrom
        File exp_bed
        File exp_bim
        File exp_fam
        File geno_bed
        File geno_bim
        File geno_fam
        File expression_matrix
        File covar
        String prefix
        String prefix_chrom
        String docker
        String zones
    }

    command <<<
        set -e

        export PYENV_ROOT="$HOME/.pyenv"
        [[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"
        eval "$(pyenv init - bash)"

        mkdir -p exp_bfile geno_bfile
        ln -s ~{exp_bed} ~{exp_bim} ~{exp_fam} exp_bfile/
        ln -s ~{geno_bed} ~{geno_bim} ~{geno_fam} geno_bfile/
        exp_bfile=exp_bfile/$(basename ~{exp_bed} .bed)
        geno_bfile=geno_bfile/$(basename ~{geno_bed} .bed)

        gunzip -c ~{covar} > ~{basename(covar, ".gz")}

        python ~/mesc/run_mesc.py \
        --compute-expscore-indiv \
        --plink-path $(which plink) \
        --expression-matrix ~{expression_matrix} \
        --exp-bfile $exp_bfile \
        --geno-bfile $geno_bfile \
        --chr ~{sub(chrom, "^chr", "")} \
        --covariates ~{basename(covar, ".gz")} \
        --out ~{prefix} \
        --gene-list ~{chunk} \
        --est-lasso-only
    >>>

    output {
        File out_hsq = "~{prefix_chrom}.hsq"
        File out_lasso = "~{prefix_chrom}.lasso"
    }

    runtime {
        docker: docker
        cpu: 4
        memory: "8 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task compute_expscore_from_lasso {
    input {
        String chrom
        File geno_bed
        File geno_bim
        File geno_fam
        Array[File] lasso
        Array[File] hsq
        String prefix
        String prefix_chrom
        String docker
        String zones
    }

    command <<<
        set -e

        export PYENV_ROOT="$HOME/.pyenv"
        [[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"
        eval "$(pyenv init - bash)"

        mkdir -p geno_bfile
        ln -s ~{geno_bed} ~{geno_bim} ~{geno_fam} geno_bfile/
        geno_bfile=geno_bfile/$(basename ~{geno_bed} .bed)

        python ~/mesc/run_mesc.py \
        --compute-expscore-from-lasso \
        --geno-bfile $geno_bfile \
        --chr ~{sub(chrom, "^chr", "")} \
        --lasso-files ~{sep=" " lasso} \
        --hsq-files ~{sep=" " hsq} \
        --out ~{prefix}
    >>>

    output {
        File out_G = "~{prefix_chrom}.G"
        File out_ave = "~{prefix_chrom}.ave_h2cis"
        File out_expscore = "~{prefix_chrom}.expscore.gz"
        File out_gannot = "~{prefix_chrom}.gannot.gz"
        File out_hsq = "~{prefix_chrom}.hsq"
        File out_lasso = "~{prefix_chrom}.lasso"
    }

    runtime {
        docker: docker
        cpu: 4
        memory: "8 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task tar_expscore {
    input {
        String prefix
        Array[File] expscore_files
        Array[File] hsq_files
        String docker
        String zones
    }

    command <<<
        set -e

        mkdir -p ~{prefix}
        cp ~{sep=' ' expscore_files} ~{prefix}/
        tar czvf ~{prefix}.tar.gz ~{prefix}

        awk 'NR == 1 || FNR > 1' ~{sep=' ' hsq_files} | bgzip -c > ~{prefix}.hsq.gz
    >>>

    output {
        File out_tar = "~{prefix}.tar.gz"
        File out_hsq = "~{prefix}.hsq.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task preprocess_phenotypes {
    input {
        File phenotype_summary
        Float h2_z_threshold = 2
        Int num_gw_significant_threshold = 10
        Int chunk_size = 25
        Boolean debug = false
        String docker
        String zones
    }

    command <<<
        set -e

        python3 << "__EOF__"
        import math
        import pandas as pd
        import re

        df = pd.read_csv("~{phenotype_summary}", sep="\t")
        df = df.loc[(df.H2_Z > ~{h2_z_threshold}) & (df.num_gw_significant > ~{num_gw_significant_threshold}), :]

        # chunk params
        chunk_size= ~{chunk_size}
        total_lines = len(df.index)
        n_chunks = math.ceil(total_lines / chunk_size)
        lines_per_chunk = math.ceil(total_lines / n_chunks)

        if ~{true='True' false='False' debug}:
            n_chunks = 1

        for i in range(n_chunks):
            start_line = i * lines_per_chunk
            end_line = min(start_line + lines_per_chunk, total_lines)
            df[["phenocode", "ldsc_munged"]].iloc[start_line:end_line].to_csv(f"pheno.{i}.tsv", sep="\t", index=False, header=False)
            df[["ldsc_munged"]].iloc[start_line:end_line].to_csv(f"sumstats.{i}.tsv", sep="\t", index=False, header=False)
        __EOF__
    >>>

    output {
        Array[File] out_phenos = glob("pheno.*.tsv")
        Array[File] out_sumstats = glob("sumstats.*.tsv")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}

task run_mesc {
    input {
        File phenotype_list
        String cell_type
        Array[File] sumstats
        File expscore
        String prefix
        String expscore_prefix
        String docker
        String zones
    }

    command <<<
        set -e

        export PYENV_ROOT="$HOME/.pyenv"
        [[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"
        eval "$(pyenv init - bash)"

        tar xzvf ~{expscore} --strip-components=1

        mkdir -p sumstats out
        ln -s ~{sep=' ' sumstats} sumstats/

        cat << "__EOF__" > script.awk
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            print "Phenotype", "Cell_Type", $0
        }
        NR > 1 {
            print phenotype, cell_type, $0
        }
        __EOF__

        while IFS=$'\t' read -r pheno sumstats_path
        do
            echo "Processing phenotype: $pheno"
            prefix=$pheno.~{expscore_prefix}
            fname=$(basename $sumstats_path)

            python ~/mesc/run_mesc.py \
            --h2med sumstats/$fname \
            --exp-chr ~{expscore_prefix} \
            --out $prefix

            awk -f script.awk -v phenotype=$pheno -v cell_type=~{cell_type} $prefix.all.h2med > out/$prefix.all.h2med
            awk -f script.awk -v phenotype=$pheno -v cell_type=~{cell_type} $prefix.categories.h2med > out/$prefix.categories.h2med
        done < ~{phenotype_list}

        awk 'NR == 1 || FNR > 1' out/*.all.h2med > ~{prefix}.all.h2med
        awk 'NR == 1 || FNR > 1' out/*.categories.h2med > ~{prefix}.categories.h2med
    >>>

    output {
        File out_all_h2med = "~{prefix}.all.h2med"
        File out_categories_h2med = "~{prefix}.categories.h2med"
        Array[File] out_log = glob("*.log")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "8 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task combine_results {
    input {
        String prefix
        Array[File] all_h2med
        Array[File] categories_h2med
        Boolean bgzip = false
        String docker
        String zones
    }

    command <<<
        set -e

        awk 'NR == 1 || FNR > 1' ~{sep=' ' all_h2med} > ~{prefix}.all.h2med
        awk 'NR == 1 || FNR > 1' ~{sep=' ' categories_h2med} > ~{prefix}.categories.h2med

        if [[ ~{bgzip} == "true" ]]
        then
            bgzip ~{prefix}.all.h2med
            bgzip ~{prefix}.categories.h2med
        fi
    >>>

    output {
        File out_all_h2med = if bgzip then "~{prefix}.all.h2med.gz" else "~{prefix}.all.h2med"
        File out_categories_h2med = if bgzip then "~{prefix}.categories.h2med.gz" else "~{prefix}.categories.h2med"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
    }
}
