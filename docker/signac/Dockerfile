FROM gcr.io/google.com/cloudsdktool/google-cloud-cli:slim AS builder

ENV HTSLIB_VERSION=1.21
RUN apt-get update --fix-missing && \
    apt-get install -y \
    libbz2-dev \
    libcurl4-openssl-dev \
    libdeflate-dev \
    liblzma-dev \
    unzip \
    wget

RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 -P ~ && \
    tar xjvf ~/htslib-${HTSLIB_VERSION}.tar.bz2 -C ~ && \
    cd ~/htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr --enable-libcurl --enable-gcs --with-libdeflate && \
    make && make install

RUN wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
    chmod +x bedGraphToBigWig && \
    mv bedGraphToBigWig /usr/bin

FROM europe-docker.pkg.dev/finngen-refinery-dev/eu.gcr.io/seurat:5.1.0-mk
LABEL maintainer="mkanai@broadinstitute.org"

WORKDIR /app
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# tools and htslib
RUN apt-get update --fix-missing && \
    apt-get install -y \
    bedtools \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    parallel \
    python3-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV MACS2_VERSION=2.2.9.1
RUN pip3 install macs2==${MACS2_VERSION}

ENV ABC_VERSION=698e361c7f184d09f43bdb310360c28ef3a252f4
RUN mkdir -p /app/ABC && \
    wget -P /app/ABC \
    https://raw.githubusercontent.com/mkanai/ABC-Enhancer-Gene-Prediction/${ABC_VERSION}/reference/hg38/GRCh38_unified_blacklist.bed \
    https://raw.githubusercontent.com/mkanai/ABC-Enhancer-Gene-Prediction/${ABC_VERSION}/workflow/scripts/makeCandidateRegions.py \
    https://raw.githubusercontent.com/mkanai/ABC-Enhancer-Gene-Prediction/${ABC_VERSION}/workflow/scripts/neighborhoods.py \
    https://raw.githubusercontent.com/mkanai/ABC-Enhancer-Gene-Prediction/${ABC_VERSION}/workflow/scripts/peaks.py \
    https://raw.githubusercontent.com/mkanai/ABC-Enhancer-Gene-Prediction/${ABC_VERSION}/workflow/scripts/tools.py && \
    pip3 install \
    pyBigWig \
    pyranges \
    pysam && \
    python3 /app/ABC/makeCandidateRegions.py --help

ENV SIGNAC_VERSION=1.14.0
RUN Rscript -e "options(repos=c(CRAN='https://cloud.r-project.org')); setRepositories(ind=1:3); install.packages('Signac')" && \
    Rscript -e "BiocManager::install(c('AnnotationHub', 'biovizBase', 'ensembldb'))"

COPY --from=builder /usr/bin/bedGraphToBigWig /usr/bin/bedGraphToBigWig
COPY --from=builder /usr/bin/bgzip /usr/bin/bgzip
COPY --from=builder /usr/bin/tabix /usr/bin/tabix

CMD [ "/bin/bash" ]
