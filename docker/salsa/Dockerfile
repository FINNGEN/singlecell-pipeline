# Modified from https://github.com/p4rkerw/SALSA/blob/main/docker/Dockerfile
FROM broadinstitute/gatk:4.6.0.0 AS builder

RUN apt-get update --fix-missing && \
    apt-get install -y \
    libbz2-dev \
    libcurl4-openssl-dev \
    libdeflate-dev \
    libgsl-dev \
    liblzma-dev \
    libperl-dev \
    pv \
    zlib1g-dev

ENV HTSLIB_VERSION=1.21
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 -P ~ && \
    tar xjvf ~/htslib-${HTSLIB_VERSION}.tar.bz2 -C ~ && \
    cd ~/htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr --enable-libcurl --enable-gcs --with-libdeflate && \
    make && make install

ENV SAMTOOLS_VERSION=1.21
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 -P ~ && \
    tar xjvf ~/samtools-${SAMTOOLS_VERSION}.tar.bz2 -C ~ && \
    cd ~/samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr --without-curses && \
    make && make install

ENV BCFTOOLS_VERSION=1.21
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 -P ~ && \
    tar xjvf ~/bcftools-${BCFTOOLS_VERSION}.tar.bz2 -C ~ && \
    cd ~/bcftools-${BCFTOOLS_VERSION} && \
    ./configure --prefix=/usr --enable-libgsl --enable-perl-filters && \
    make && make install && \
    rm -r ~/bcftools-${BCFTOOLS_VERSION}*

ENV BWA_VERSION=0.7.18
RUN wget https://github.com/lh3/bwa/archive/refs/tags/v${BWA_VERSION}.tar.gz -P ~ && \
    tar xzvf ~/v${BWA_VERSION}.tar.gz -C ~ && \
    cd ~/bwa-${BWA_VERSION} && \
    make && \
    mv bwa /usr/bin && \
    rm -r ~/bwa-${BWA_VERSION}*

ENV STAR_VERSION=2.7.11b
RUN wget https://github.com/alexdobin/STAR/releases/download/${STAR_VERSION}/STAR_${STAR_VERSION}.zip && \
    unzip STAR_${STAR_VERSION}.zip && \
    mv STAR_${STAR_VERSION}/Linux_x86_64/STAR /usr/bin/STAR

ENV SUBSET_BAM_VERSION=1.1.0
RUN wget https://github.com/10XGenomics/subset-bam/releases/download/v${SUBSET_BAM_VERSION}/subset-bam_linux && \
  chmod +x subset-bam_linux && \
  mv subset-bam_linux /usr/bin/subset-bam

RUN wget https://raw.githubusercontent.com/aertslab/single_cell_toolkit/refs/heads/master/subset_bam_per_cb.sh && \
    chmod +x subset_bam_per_cb.sh && \
    mv subset_bam_per_cb.sh /usr/bin/subset_bam_per_cb.sh

FROM broadinstitute/gatk:4.6.0.0
LABEL maintainer="mkanai@broadinstitute.org"

WORKDIR /app
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update --fix-missing && \
    apt-get install -y \
    libbz2-dev \
    libcurl4-openssl-dev \
    libdeflate-dev \
    libgsl-dev \
    liblzma-dev \
    libperl-dev \
    pv \
    zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder \
    /usr/bin/bcftools \
    /usr/bin/bgzip \
    /usr/bin/bwa \
    /usr/bin/samtools \
    /usr/bin/STAR \
    /usr/bin/subset-bam \
    /usr/bin/subset_bam_per_cb.sh \
    /usr/bin/tabix \
    /usr/bin/
COPY --from=builder /usr/lib/libhts.* /usr/lib/
COPY --from=builder /usr/lib/pkgconfig/htslib.pc /usr/lib/pkgconfig/
COPY --from=builder /usr/include/htslib /usr/include/

ENV BCFTOOLS_PLUGINS=/usr/libexec/bcftools/
COPY --from=builder /usr/libexec/bcftools/ /usr/libexec/bcftools/

ENV LD_LIBRARY_PATH="/lib/x86_64-linux-gnu:/opt/miniconda/lib"

# install r bioconda pysam bwa STAR hdf5 for WASP dependencies
# pysam hdf5 required for WASP snp2h5 install and compile
RUN conda config --add channels r && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install hdf5 && \
    pip install pysam && \
    pip install tables

# clone the wasp v0.3.4 repo and compile snp2h5 in /opt
RUN cd /opt && \
  git clone -b v0.3.4 https://github.com/bmvdgeijn/WASP.git && \
  cat WASP/snp2h5/Makefile | sed 's/$(HDF_INSTALL)/\/opt\/miniconda/g' > WASP/snp2h5/Makefile2 && \
  mv WASP/snp2h5/Makefile2 WASP/snp2h5/Makefile && \
  export PATH=/opt/miniconda/bin/lib:$PATH && \
  cd WASP/snp2h5 && \
  make && \
  ln -s /opt/WASP/snp2h5/snp2h5 /usr/bin/snp2h5

ENV SALSA_VERSION=c60d4cafb129d11d2ea0a32f71d0081297cfefea
RUN git clone https://github.com/mkanai/SALSA && \
    cd SALSA && \
    git checkout ${SALSA_VERSION} && \
    rm -r .git Tutorials