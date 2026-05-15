version 1.0

workflow scanpy_umap {
    input {
        File h5ad
        String docker_suite
        String zones
    }

    call run_umap {
        input:
            h5ad = h5ad,
            docker = docker_suite,
            zones = zones
    }
}

task run_umap {
    input {
        File h5ad
        String prefix = basename(h5ad, ".h5ad")
        String docker
        String zones
    }

    command <<<
        set -e

        export n_cpu=$(grep -c ^processor /proc/cpuinfo)

        cat << "__EOF__" > script.py
        import numpy as np
        import pandas as pd
        import anndata as ad
        import scanpy as sc


        adata = ad.read_h5ad("~{h5ad}")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata._inplace_subset_var(adata.var.highly_variable)
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pp.neighbors(adata)
        sc.tl.umap(adata, random_state=42)

        adata.X = None
        adata.write_h5ad("~{prefix}.umap.h5ad")
        __EOF__

        python3 script.py && \
        touch _SUCCESS

    >>>

    output {
        File out_success = "_SUCCESS"
        File out_h5ad = "~{prefix}.umap.h5ad"
    }

    runtime {
        docker: docker
        cpu: 96
        memory: "624 GB"
        disks: "local-disk 300 HDD"
        zones: zones
        preemptible: 2
    }
}
