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

FROM europe-docker.pkg.dev/finngen-refinery-dev/eu.gcr.io/seurat:5.1.0-mk
LABEL maintainer="mkanai@broadinstitute.org"

WORKDIR /app
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# tools and htslib
RUN apt-get update --fix-missing && \
    apt-get install -y \
    bedtools \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    parallel \
    python3-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV MACS2_VERSION=2.2.9.1
RUN pip3 install macs2==${MACS2_VERSION}

ENV SIGNAC_VERSION=1.0.2
RUN Rscript -e "remotes::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories())" && \
    Rscript -e "remotes::install_github('PhanstielLab/bedtoolsr')" && \
    Rscript -e "BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38'))" && \
    Rscript -e "ArchR::installExtraPackages()"

COPY --from=builder /usr/bin/bgzip /usr/bin/bgzip
COPY --from=builder /usr/bin/tabix /usr/bin/tabix

COPY R/iterative_overlap_peak_merge.R /usr/local/bin/

RUN chmod a+x /usr/local/bin/*.R

CMD [ "/bin/bash" ]
