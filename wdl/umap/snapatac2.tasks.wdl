version 1.0

task prepare_h5ad_list {
    input {
        String obs
        String h5ad_pattern
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        python3 << "__EOF__"
        import pandas as pd
        import re

        df = pd.read_csv("~{obs}", sep="\t")
        df["FINNGENID2"] = df["barcode"].str.split("-", n=2).str[0]
        df["h5ad"] = df["FINNGENID2"].apply(lambda x: re.sub(r"\{FINNGENID\}", x, "~{h5ad_pattern}"))
        df = df[["FINNGENID", "h5ad"]].drop_duplicates().reset_index(drop=True)
        df = df.groupby("FINNGENID")["h5ad"].apply(lambda x: '\t'.join(x)).reset_index(drop=True)
        df.to_csv("h5ad.tsv", index=False, header=False)
        __EOF__
    >>>

    output {
        Array[Array[File]] out_h5ad = read_tsv("h5ad.tsv")
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

task convert_to_csr {
    input {
        Array[File] h5ad
        String prefix = sub(basename(h5ad[0], '.h5ad'), "_[0-9]+\\.", ".")
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        python3 << "__EOF__"
        import anndata as ad

        h5ad_files = "~{sep=',' h5ad}".split(",")

        res = []
        for h5ad in h5ad_files:
            adata = ad.read_h5ad(h5ad)
            res.append(adata)
            print(f"{len(res)} / {len(h5ad_files)}")

        adata = ad.concat(res, axis = 0, join="outer")
        adata.X = adata.X.tocsr()
        adata.write_h5ad("~{prefix}.csr.h5ad")
        __EOF__
    >>>

    output {
        File out_h5ad = "~{prefix}.csr.h5ad"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{ceil(3 * size(h5ad, 'GB')) + 10} GB"
        disks: "local-disk ~{ceil(2 * size(h5ad, 'GB')) + 10} HDD"
        zones: zones
        preemptible: 2
    }
}
