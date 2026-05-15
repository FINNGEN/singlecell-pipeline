version 1.0

# Shared SALSA tasks used by both salsa.atac.wdl and salsa.rna.wdl.
# Modality-specific behavior is controlled by the `modality` input to salsa_wasp
# and salsa_count ("atac" or "rna"). The reference tarball's top-level directory
# name is passed via `reference_dir`:
#   ATAC: refdata-cellranger-arc-GRCh38-2020-A-2.0.0 (BWA index + fasta)
#   RNA:  refdata-gex-GRCh38-2020-A                  (STAR index + fasta + .dict)
# ATAC-only: make_fragment_file lives here for convenience; it's called only
# from salsa.atac.wdl.

task prepare_params {
    input {
        File obs
        String vcf_pattern
        String bam_pattern
        String counts_pattern
        String? fragments_pattern
        Boolean overwrite
        Int chunk_size = 200
        String docker
        String zones
    }

    String fragments_pattern_nullable = select_first([fragments_pattern, ""])

    command <<<
        set -e -o pipefail

        # for each pool_ffid
        python3 << "__EOF__"
        import math
        import numpy as np
        import pandas as pd
        import re
        from google.cloud import storage

        FRAGMENTS_PATTERN = "~{fragments_pattern_nullable}"
        has_fragments = bool(FRAGMENTS_PATTERN)

        df = pd.read_csv("~{obs}", sep="\t")
        df["FINNGENID2"] = df["barcode"].str.split("-", n=2).str[0]
        df = df[["FINNGENID", "FINNGENID2", "ffid_2"]].drop_duplicates()
        df["vcf"] = df["FINNGENID"].apply(lambda x: re.sub(r"\{FINNGENID\}", x, "~{vcf_pattern}"))
        df["bam"] = df["ffid_2"].apply(lambda x: re.sub(r"\{FFID\}", x, "~{bam_pattern}"))
        if has_fragments:
            df["fragments"] = df["FINNGENID2"].apply(lambda x: re.sub(r"\{FINNGENID\}", x, FRAGMENTS_PATTERN))
        df["counts"] = df["FINNGENID2"].apply(lambda x: re.sub(r"\{FINNGENID\}", x, "~{counts_pattern}"))
        cols = ["FINNGENID2", "vcf", "bam"] + (["fragments"] if has_fragments else []) + ["counts"]
        df = df[cols]

        storage_client = storage.Client()

        # filter out non-existent bams
        bucket = storage_client.get_bucket("~{bam_pattern}".split("/")[2])
        bam_exists = np.array([bucket.blob(x.split("/", 3)[-1]).exists() for x in df.bam])
        print(f"{np.sum(bam_exists)}/{len(bam_exists)} bams exist")

        # filter out finished samples
        bucket = storage_client.get_bucket("~{counts_pattern}".split("/")[2])
        counts_exists = np.array([bucket.blob(x.split("/", 3)[-1]).exists() for x in df.counts])
        print(f"{np.sum(counts_exists)}/{len(counts_exists)} count tables exist")

        if has_fragments:
            bucket = storage_client.get_bucket(FRAGMENTS_PATTERN.split("/")[2])
            fragments_exists = np.array([bucket.blob(x.split("/", 3)[-1]).exists() for x in df.fragments])
            print(f"{np.sum(fragments_exists)}/{len(fragments_exists)} fragments exist")

        overwrite = ~{true='True' false='False' overwrite}
        if has_fragments:
            keep = bam_exists & (overwrite | ~fragments_exists | ~counts_exists)
        else:
            keep = bam_exists & (overwrite | ~counts_exists)
        df = df.loc[keep, :]
        df.to_csv("params.tsv", sep="\t", index=False, header=False)
        print(f"Processing {len(df)} samples for SALSA")

        # chunk params
        chunk_size = ~{chunk_size}
        total_lines = len(df.index)
        n_chunks = max(1, math.ceil(total_lines / chunk_size))
        lines_per_chunk = math.ceil(total_lines / n_chunks) if total_lines > 0 else 1

        for i in range(n_chunks):
            start_line = i * lines_per_chunk
            end_line = min(start_line + lines_per_chunk, total_lines)
            df.iloc[start_line:end_line].to_csv(f"params.{i}.tsv", sep="\t", index=False, header=False)
        __EOF__
    >>>

    output {
        File out_params = "params.tsv"
        Array[File] out_param_files = glob("params.*.tsv")
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

task salsa_wasp {
    input {
        String fgid
        File vcf
        File bam
        File reference
        File chrom_info
        String modality       # "atac" or "rna"
        String reference_dir  # top-level dir name inside reference tar
        Int chunk_size = 1000
        Int cpu = 8
        # Memory floor: STAR (RNA) needs ~30 GB for the genome load;
        # BWA (ATAC) is proportional to bam size. The larger of the two wins.
        Int memory_gb_floor = 16
        String docker
        String zones
    }

    Int memory_gb = if ceil(size(bam, 'GB')) + 10 > memory_gb_floor then ceil(size(bam, 'GB')) + 10 else memory_gb_floor

    command <<<
        set -e -o pipefail
        n_cpu=$(grep -c ^processor /proc/cpuinfo)
        export TMPDIR=$(pwd) SCRATCH1=$(pwd)

        tar xzvf ~{reference}
        mkdir -p reference && mv ~{chrom_info} reference/

        if [ "~{modality}" = "atac" ]; then
            ref_flag="--atacref ~{reference_dir}"
        else
            ref_flag="--stargenome ~{reference_dir}/star"
        fi

        bash /app/SALSA/step6_wasp.sh \
        --inputvcf ~{vcf} \
        --inputbam ~{bam} \
        --outputdir ~{fgid} \
        --outputbam ~{fgid}.wasp.bam \
        $ref_flag \
        --genotype ~{modality} \
        --library_id ~{fgid} \
        --modality ~{modality} \
        --isphased \
        --threads $n_cpu

        # chunk bam barcodes
        samtools view ~{fgid}/~{fgid}.wasp.bam \
        | cut -f 12- \
        | tr "\t" "\n" \
        | grep  "^CB:Z:" \
        | cut -d ':' -f3 \
        | sort | uniq > ~{fgid}.wasp.barcodes.txt

        python3 << "__EOF__"
        import math

        with open("~{fgid}.wasp.barcodes.txt", 'r') as file:
            lines = file.readlines()

        chunk_size = ~{chunk_size}
        total_lines = len(lines)
        n_chunks = max(1, math.ceil(total_lines / chunk_size))
        lines_per_chunk = math.ceil(total_lines / n_chunks) if total_lines > 0 else 1

        for i in range(n_chunks):
            with open(f"~{fgid}.wasp.{i:02d}.barcodes.txt", 'w') as file:
                start_line = i * lines_per_chunk
                end_line = min(start_line + lines_per_chunk, total_lines)
                file.writelines(lines[start_line:end_line])
        __EOF__
    >>>

    output {
        File out_bam = "~{fgid}/~{fgid}.wasp.bam"
        File out_bai = "~{fgid}/~{fgid}.wasp.bam.bai"
        File out_barcodes = "~{fgid}.wasp.barcodes.txt"
        Array[File] out_chunked_barcodes = glob("~{fgid}.wasp.*.barcodes.txt")
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{ceil(2 * size(vcf, 'GB') + 15 * size(bam, 'GB') + 3 * size(reference, 'GB')) + 20} HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task chunk_wasp_bam {
    input {
        String bam
        File barcodes
        String prefix = basename(barcodes, '.barcodes.txt')
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        n_cpu=$(grep -c ^processor /proc/cpuinfo)
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view \
        -D CB:~{barcodes} \
        --write-index \
        -@ $n_cpu \
        -o ~{prefix}.bam##idx##~{prefix}.bam.bai \
        ~{bam} && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_bam = "~{prefix}.bam"
        File out_bai = "~{prefix}.bam.bai"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "8 GB"
        disks: "local-disk ~{ceil(2 * size(bam, 'GB')) + 10} HDD"
        zones: zones
        preemptible: 2
        maxRetries: 1
        noAddress: true
    }
}

task salsa_count {
    input {
        String fgid
        File vcf
        File vcf_tbi = "~{vcf}.tbi"
        File bam
        File bai
        File reference
        String modality
        String reference_dir
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail
        n_cpu=$(grep -c ^processor /proc/cpuinfo)
        export TMPDIR=$(pwd) SCRATCH1=$(pwd)

        tar xzvf ~{reference}

        bash /app/SALSA/step7_gatk_alleleCount.sh \
        --inputvcf ~{vcf} \
        --inputbam ~{bam} \
        --outputdir ~{fgid} \
        --genotype ~{modality} \
        --library_id ~{fgid} \
        --modality ~{modality} \
        --reference ~{reference_dir} \
        --single_cell_counts \
        --isphased \
        --threads $n_cpu

        if [[ -s ~{fgid}/~{fgid}.counts.single_cell.table ]]
        then
            touch _SUCCESS
        fi
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_single_cell_table = "~{fgid}/~{fgid}.counts.single_cell.table"
        File out_single_cell_phased_table = "~{fgid}/~{fgid}.counts.single_cell.phased.table"
    }

    runtime {
        docker: docker
        cpu: 8
        memory: "16 GB"
        disks: "local-disk ~{ceil(2 * size(vcf, 'GB') + 30 * size(bam, 'GB') + 2 * size(reference, 'GB')) + 10} HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

task merge_counts {
    input {
        String fgid
        Array[File] single_cell_tables
        Array[File] single_cell_phased_tables
        String output_directory
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        # When `cat`-piped, awk sees a single stream so FNR==NR; using FNR>1 would
        # leak headers from files 2..N into the data. Match on the literal header
        # column name instead.
        cat ~{sep=" " single_cell_tables} | awk '
        NR == 1 {print; next}
        $1 == "contig" {next}
        {print | "sort -k1,1V -k2,2n"}
        ' | bgzip -c > ~{fgid}.counts.single_cell.table.gz

        cat ~{sep=" " single_cell_phased_tables} | awk '
        NR == 1 {print; next}
        $1 == "variant_id" {next}
        {print | "sort -k1,1V -k16,16"}
        ' | bgzip -c > ~{fgid}.counts.single_cell.phased.table.gz

        gcloud storage cp ~{fgid}.counts.single_cell.table.gz ~{fgid}.counts.single_cell.phased.table.gz ~{output_directory}/ && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_single_cell_table = "~{fgid}.counts.single_cell.table.gz"
        File out_single_cell_phased_table = "~{fgid}.counts.single_cell.phased.table.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "16 GB"
        disks: "local-disk ~{2 * ceil(size(single_cell_tables, 'GB') + size(single_cell_phased_tables, 'GB')) + 10} HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}

# ATAC-only: reconstruct a fragments BED from the paired-end WASP BAM.
# Lives here for convenience; only called from salsa.atac.wdl.
task make_fragment_file {
    input {
        File bam
        String output_directory
        String prefix = basename(bam, '.bam')
        String docker
        String zones
    }

    command <<<
        set -e -o pipefail

        cat << "__EOF__" > script.py
        import snapatac2
        snapatac2.pp.make_fragment_file(
            bam_file="~{bam}",
            output_file="~{prefix}.fragments.bed",
            is_paired=True,
            barcode_tag="CB",
            source="10x"
        )
        __EOF__

        python3 script.py && \
        sort -V -k1,1 -k2,2 -k3,3 ~{prefix}.fragments.bed | grep '^chr' | bgzip -c > ~{prefix}.fragments.bed.gz && \
        tabix -p bed ~{prefix}.fragments.bed.gz && \
        gcloud storage cp ~{prefix}.fragments.bed.gz* ~{output_directory}/ && \
        touch _SUCCESS
    >>>

    output {
        File out_success = "_SUCCESS"
        File out_fragments = "~{prefix}.fragments.bed.gz"
        File out_fragments_tbi = "~{prefix}.fragments.bed.gz.tbi"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "16 GB"
        disks: "local-disk ~{2 * ceil(size(bam, 'GB')) + 10} HDD"
        zones: zones
        preemptible: 2
        noAddress: true
    }
}
