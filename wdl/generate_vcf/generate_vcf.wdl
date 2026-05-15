version 1.0

workflow generate_vcf {
    input {
        Array[String] chromosomes
        File mapping
        File? subset_samples
        String outname
        File bgen_list
        Array[String] bgen_files = read_lines(bgen_list)
        File genotype_vcf_list
        Array[String] genotype_vcfs = read_lines(genotype_vcf_list)
        File info
        Boolean update_ids
        String docker
        String docker_suite
        String zones
    }

    scatter (bgen in bgen_files) {
        call extract_samples {
            input:
                bgen = bgen,
                mapping = mapping,
                info = info,
                subset_samples = subset_samples,
                update_ids = update_ids,
                docker = docker,
                zones = zones
        }
    }

    call concat_vcf {
        input:
            vcfs = extract_samples.out_vcf,
            outname = outname,
            docker = docker,
            zones = zones
    }

    call make_grm {
        input:
            pgen = concat_vcf.out_pgen,
            pvar = concat_vcf.out_pvar,
            psam = concat_vcf.out_psam,
            docker = docker,
            zones = zones
    }

    scatter (chrom in chromosomes) {
        call chunk_chrom {
            input:
                pgen = concat_vcf.out_pgen,
                pvar = concat_vcf.out_pvar,
                psam = concat_vcf.out_psam,
                chrom = chrom,
                docker = docker,
                zones = zones
        }

        scatter (snplist in chunk_chrom.out_snplist) {
            call chunk_bgen {
                input:
                    pgen = chunk_chrom.out_pgen,
                    pvar = chunk_chrom.out_pvar,
                    psam = chunk_chrom.out_psam,
                    snplist = snplist,
                    docker = docker,
                    zones = zones
            }
        }
    }

    # for pca_kinship and rasqual
    scatter (genotype_vcf in genotype_vcfs) {
        call subset_genotype_vcf {
            input:
                genotype_vcf = genotype_vcf,
                mapping = mapping,
                # finngen_R12_chr1.vcf.gz
                outname = "~{outname}.~{sub(sub(basename(genotype_vcf, '.vcf.gz'), 'finngen_R12_', ''), 'chr23', 'chrX')}",
                docker = docker,
                zones = zones
        }

        # RASQUAL
        call filter_het_snps {
            input:
                genotype_vcf = subset_genotype_vcf.out_vcf,
                docker = docker_suite,
                zones = zones
        }
    }

    # RASQUAL
    call combine_site_vcfs {
        input:
            site_vcfs = filter_het_snps.out_vcf_sites,
            outname = outname,
            docker = docker_suite,
            zones = zones
    }

    call concat_phased_vcfs {
        input:
            vcfs = subset_genotype_vcf.out_vcf,
            outname = outname,
            docker = docker,
            zones = zones
    }

    scatter (fgid in concat_phased_vcfs.out_samples) {
        call subset_phased_vcfs {
            input:
                vcf = concat_phased_vcfs.out_vcf,
                fgid = fgid,
                docker = docker_suite,
                zones = zones
        }
    }

    output {
        File out_vcf = concat_vcf.out_vcf
        File out_vcf_tbi = concat_vcf.out_vcf_tbi
        File out_lexsorted_vcf = concat_vcf.out_lexsorted_vcf
        File out_lexsorted_vcf_tbi = concat_vcf.out_lexsorted_vcf_tbi
        File out_pgen = concat_vcf.out_pgen
        File out_pvar = concat_vcf.out_pvar
        File out_psam = concat_vcf.out_psam
        File out_grm_bed = make_grm.out_bed
        File out_grm_bim = make_grm.out_bim
        File out_grm_fam = make_grm.out_fam
        Array[File] out_per_chrom_pgen = chunk_chrom.out_pgen
        Array[File] out_per_chrom_pvar = chunk_chrom.out_pvar
        Array[File] out_per_chrom_psam = chunk_chrom.out_psam
        Array[File] out_chunk_bgen = flatten(select_all(chunk_bgen.out_bgen))
        Array[File] out_chunk_bgi = flatten(select_all(chunk_bgen.out_bgi))
        Array[File] out_chunk_sample = flatten(select_all(chunk_bgen.out_sample))
    }
}

task extract_samples {
    input {
        File bgen
        File bgen_bgi = sub(bgen, ".bgen$", ".bgen.bgi")
        File bgen_sample = sub(bgen, ".bgen$", ".bgen.sample")
        File mapping
        File info
        File? subset_samples
        String prefix = sub(basename(bgen, ".bgen"), "_23.", "_X.")
        Boolean update_ids

        String docker
        String zones
    }

    command <<<
        set -e

        tail -n+2 ~{mapping} > mapping.txt
        awk '{print $1, $1, $2, $2}' mapping.txt > update_ids.txt

        if ~{true='true' false='false' defined(subset_samples)}
        then
            awk '
            BEGIN {
                OFS = "\t"
            }
            FNR == NR {
                s[$1] = 1
                next
            }
            FNR < NR && $1 in s {
                print $1, $1 > "keep.txt"
            }
            ' ~{subset_samples} mapping.txt
            cat keep.txt
        else
            awk '{print $1, $1}' mapping.txt > keep.txt
        fi

        zcat ~{info} | grep -v '^ID' | sed -e 's/:/_/g' -e 's/^23/X/' -e 's/^/chr/' > info.extract

        plink2 \
        --bgen ~{bgen} ref-first \
        --sample ~{bgen_sample} \
        --keep keep.txt \
        --extract info.extract \
        ~{true='--update-ids update_ids.txt' false='' update_ids} \
        --mac 1 \
        --export vcf bgz id-paste='iid' vcf-dosage='GP' \
        --output-chr chrMT \
        --out from_bgen

        bcftools +fill-tags -Oz from_bgen.vcf.gz -- -t AN,AC,AF > ~{prefix}.vcf.gz
        tabix ~{prefix}.vcf.gz
    >>>

    output {
        File out_vcf = prefix + ".vcf.gz"
        File out_vcf_tbi = prefix + ".vcf.gz.tbi"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "16 GB"
        disks: "local-disk 20 HDD"
        zones: zones
        preemptible: 2
    }
}

task concat_vcf {
    input {
        Array[File] vcfs
        String outname

        String docker
        String zones
    }

    command <<<
        set -e

        export TMPDIR=$(pwd)

        python << "__EOF__"  > files.txt
        import os.path
        import re

        files = "~{sep=' ' vcfs}".split(" ")
        # Process filenames like finngen_R12_1.0.vcf.gz or finngen_R12_FIXED_PAR_X_both_PAR.vcf.gz
        def parse_filename(filepath):
            filename = os.path.basename(filepath)
            # Extract the part after the prefix
            parsed = re.sub(r"^finngen_R12_(?:FIXED_PAR_)?", "", filename)
            if parsed.startswith("X_"):
                # Handle X chromosome special case
                return [filepath, ["X", "100"]]
            else:
                # Regular case like 1.0.vcf.gz
                parts = parsed.split(".")
                return [filepath, parts[:2]]

        s = [parse_filename(x) for x in files]
        # sort by lex order for the autosomal chroms, then numerical order for the batches
        auto_lexsorted = sorted(filter(lambda x: x[1][0] != "X", s), key = lambda x: [x[1][0], int(x[1][1])])
        auto_numsorted = sorted(filter(lambda x: x[1][0] != "X", s), key = lambda x: [int(x[1][0].replace("X", "23")), int(x[1][1])])
        chrX = filter(lambda x: x[1][0] == "X", s)

        with open("files_auto_lexsorted.txt", "w") as f:
            f.write('\n'.join([x for x, _ in auto_lexsorted]))
        with open("files_auto_numsorted.txt", "w") as f:
            f.write('\n'.join([x for x, _ in auto_numsorted]))
        with open("files_chrX.txt", "w") as f:
            f.write('\n'.join([x for x, _ in chrX]))
        __EOF__

        echo "files_auto_lexsorted.txt"
        cat files_auto_lexsorted.txt
        echo "files_auto_numsorted.txt"
        cat files_auto_numsorted.txt
        echo "files_chrX.txt"
        cat files_chrX.txt

        bcftools concat -f files_auto_numsorted.txt | bgzip -c > ~{outname}.autosomes.vcf.gz
        bcftools concat -f files_chrX.txt | bcftools sort --temp-dir $TMPDIR/bcftools - | bgzip -c > ~{outname}.chrX.vcf.gz

        bcftools concat -f files_auto_lexsorted.txt | bgzip -c > ~{outname}.lexsorted.vcf.gz
        bcftools concat ~{outname}.autosomes.vcf.gz ~{outname}.chrX.vcf.gz | bgzip -c > ~{outname}.vcf.gz

        tabix ~{outname}.lexsorted.vcf.gz
        tabix ~{outname}.vcf.gz

        plink2 \
        --vcf ~{outname}.vcf.gz 'dosage=GP' \
        --output-chr chrMT \
        --mac 1 \
        --make-pfile \
        --out ~{outname}

        plink2 \
        --pfile ~{outname} \
        --make-bed \
        --out ~{outname}
    >>>

    output {
        File out_vcf = "~{outname}.vcf.gz"
        File out_vcf_tbi = "~{outname}.vcf.gz.tbi"
        File out_lexsorted_vcf = "~{outname}.lexsorted.vcf.gz"
        File out_lexsorted_vcf_tbi = "~{outname}.lexsorted.vcf.gz.tbi"
        File out_pgen = "~{outname}.pgen"
        File out_pvar = "~{outname}.pvar"
        File out_psam = "~{outname}.psam"
        File out_bed = "~{outname}.bed"
        File out_bim = "~{outname}.bim"
        File out_fam = "~{outname}.fam"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "14 GB"
        disks: "local-disk 100 HDD"
        zones: zones
        preemptible: 2
    }
}

task make_grm {
    input {
        File pgen
        File pvar
        File psam
        File include_variants
        Float geno_missing
        Float maf
        String ld
        String outname = basename(pgen, ".pgen") + "_GRM_LD_0.2"

        String docker
        String zones
    }

    command <<<
        set -e

        plink2 \
        --pgen ~{pgen} \
        --pvar ~{pvar} \
        --psam ~{psam} \
        --extract ~{include_variants} \
        --geno ~{geno_missing} \
        --maf ~{maf} \
        --indep-pairwise ~{ld} \
        --chr 1-22 \
        --out ~{outname}

        plink2 \
        --pgen ~{pgen} \
        --pvar ~{pvar} \
        --psam ~{psam} \
        --extract ~{outname}.prune.in \
        --make-bed \
        --out ~{outname}

    >>>

    output {
        File out_bed = outname + ".bed"
        File out_bim = outname + ".bim"
        File out_fam = outname + ".fam"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "14 GB"
        disks: "local-disk 100 HDD"
        zones: zones
        preemptible: 2
    }
}

task chunk_chrom {
    input {
        File pgen
        File pvar
        File psam
        String chrom
        String outname = basename(pgen, ".pgen") + "." + chrom
        Int chunk_size

        String docker
        String zones
    }

    command <<<
        set -e

        # chunk pfile
        plink2 \
        --pgen ~{pgen} \
        --pvar ~{pvar} \
        --psam ~{psam} \
        --chr ~{chrom} \
        --make-pfile \
        --out ~{outname}

        # export bgen
        plink2 \
        --pfile ~{outname} \
        --output-chr chrMT \
        --export bgen-1.2 ref-first bits=8 id-paste=iid sample-v2 \
        --out ~{outname}

        awk '
        NR == 1 {
            print $1"_1", $1"_2", $2, $3
        }
        NR > 1 {
            print $1, $1, $2, $3
        }' ~{outname}.sample > updated.sample
        mv updated.sample ~{outname}.sample

        bgenix -index -g ~{outname}.bgen

        # export chunked snplist
        grep -v "#" ~{outname}.pvar | cut -f3 > ~{chrom}.txt
        # split -l ~{chunk_size} -d --additional-suffix .txt ~{chrom}.txt ~{outname}.
        python3 << "__EOF__"
        import math

        with open("~{chrom}.txt", 'r') as file:
            lines = file.readlines()

        chunk_size= ~{chunk_size}
        total_lines = len(lines)
        num_chunks = math.ceil(total_lines / chunk_size)
        lines_per_chunk = math.ceil(total_lines / num_chunks)

        for i in range(num_chunks):
            with open(f"~{outname}.{i:02d}.txt", 'w') as file:
                start_line = i * lines_per_chunk
                end_line = min(start_line + lines_per_chunk, total_lines)
                file.writelines(lines[start_line:end_line])
        __EOF__
    >>>

    output {
        File out_pgen = outname + ".pgen"
        File out_pvar = outname + ".pvar"
        File out_psam = outname + ".psam"
        File out_bgen = outname + ".bgen"
        File out_bgi = outname + ".bgen.bgi"
        File out_sample = outname + ".sample"
        Array[File] out_snplist = glob(outname + ".*.txt")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "14 GB"
        disks: "local-disk 100 HDD"
        zones: zones
        preemptible: 2
    }
}

task chunk_bgen {
    input {
        File pgen
        File pvar
        File psam
        File snplist
        String outname = basename(snplist, ".txt")

        String docker
        String zones
    }

    command <<<
        set -e

        plink2 \
        --pgen ~{pgen} \
        --pvar ~{pvar} \
        --psam ~{psam} \
        --output-chr chrMT \
        --mac 1 \
        --extract ~{snplist} \
        --export bgen-1.2 ref-first bits=8 id-paste=iid sample-v2 \
        --out ~{outname}

        awk '
        NR == 1 {
            print $1"_1", $1"_2", $2, $3
        }
        NR > 1 {
            print $1, $1, $2, $3
        }' ~{outname}.sample > updated.sample
        mv updated.sample ~{outname}.sample

        bgenix -index -g ~{outname}.bgen

    >>>

    output {
        File out_bgen = outname + ".bgen"
        File out_bgi = outname + ".bgen.bgi"
        File out_sample = outname + ".sample"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "14 GB"
        disks: "local-disk 20 HDD"
        zones: zones
        preemptible: 2
    }
}

task subset_genotype_vcf {
    input {
        File genotype_vcf
        File mapping
        String outname
        Int extra_disk_space = 10
        String docker
        String zones
    }
    Int disk_space = ceil(extra_disk_space + 1.2 * size(genotype_vcf, "GB"))

    command <<<
        set -e -o pipefail

        tail -n+2 ~{mapping} | cut -f1 | sort | uniq > keep.txt

        bcftools view -Oz -i 'INFO/INFO > 0.6' -S keep.txt --force-samples ~{genotype_vcf} | \
        bcftools annotate -Oz -x INFO > ~{outname}.vcf.gz
        tabix ~{outname}.vcf.gz
    >>>

    output {
        File out_vcf = outname + ".vcf.gz"
        File out_vcf_tbi = outname + ".vcf.gz.tbi"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "64 GB"
        disks: "local-disk ~{disk_space} HDD"
        zones: zones
        preemptible: 0
    }
}

task filter_het_snps {
    input {
        String genotype_vcf
        Int min_mac = 20
        String outname = "~{basename(genotype_vcf, '.vcf.gz')}.snps.het.mac~{min_mac}"
        Int extra_disk_space = 10
        String docker
        String zones
    }
    Int disk_space = ceil(extra_disk_space + 1.2 * size(genotype_vcf, "GB"))

    command <<<
        set -e -o pipefail
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        bcftools view -Oz -v snps -g het --min-ac ~{min_mac}:minor ~{genotype_vcf} > ~{outname}.vcf.gz
        tabix ~{outname}.vcf.gz

        bcftools view -Oz -HG ~{outname}.vcf.gz > ~{outname}.sites.vcf.gz
    >>>

    output {
        File out_vcf = "~{outname}.vcf.gz"
        File out_vcf_tbi = "~{outname}.vcf.gz.tbi"
        File out_vcf_sites = "~{outname}.sites.vcf.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "16 GB"
        disks: "local-disk ~{disk_space} HDD"
        zones: zones
        preemptible: 0
    }
}

task concat_phased_vcfs {
    input {
        Array[File] vcfs
        String outname
        Int extra_disk_space = 10
        String docker
        String zones
    }
    Int disk_space = ceil(extra_disk_space + 2 * size(vcfs, "GB"))

    command <<<
        set -e -o pipefail
        echo ~{sep=" " vcfs} | tr " " "\n" | sort -V > list.txt

        bcftools concat -Oz -f list.txt > ~{outname}.vcf.gz
        tabix ~{outname}.vcf.gz

        bcftools query -l ~{outname}.vcf.gz | sort > ~{outname}.samples
    >>>

    output {
        File out_vcf = "~{outname}.vcf.gz"
        File out_vcf_tbi = "~{outname}.vcf.gz.tbi"
        Array[String] out_samples = read_lines("~{outname}.samples")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "16 GB"
        disks: "local-disk ~{disk_space} HDD"
        zones: zones
        preemptible: 0
    }
}

task subset_phased_vcfs {
    input {
        File vcf
        String fgid
        Int extra_disk_space = 10
        String docker
        String zones
    }
    Int disk_space = ceil(extra_disk_space + size(vcf, "GB"))

    command <<<
        set -e -o pipefail

        bcftools view -Oz -v snps -g het --min-ac 1:minor -s ~{fgid} ~{vcf} > ~{fgid}.vcf.gz
        tabix ~{fgid}.vcf.gz
    >>>

    output {
        File out_vcf = "~{fgid}.vcf.gz"
        File out_vcf_tbi = "~{fgid}.vcf.gz.tbi"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "4 GB"
        disks: "local-disk ~{disk_space} HDD"
        zones: zones
        preemptible: 0
    }
}

task combine_site_vcfs {
    input {
        Array[File] site_vcfs
        String outname
        Int extra_disk_space = 10
        String docker
        String zones
    }
    Int disk_space = ceil(extra_disk_space + 2 * size(site_vcfs, "GB"))

    command <<<
        set -e -o pipefail
        zcat ~{sep=" " site_vcfs} | bgzip -c > ~{outname}.all.snps.het.mac20.sites.vcf.gz
    >>>

    output {
        File out_vcf = "~{outname}.all.snps.het.mac20.sites.vcf.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "16 GB"
        disks: "local-disk ~{disk_space} HDD"
        zones: zones
        preemptible: 0
    }
}