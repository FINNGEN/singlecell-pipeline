FROM gcr.io/google.com/cloudsdktool/google-cloud-cli:slim AS builder

RUN apt-get update --fix-missing && \
    apt-get install -y \
    libbz2-dev \
    libcurl4-openssl-dev \
    libdeflate-dev \
    libgsl-dev \
    liblzma-dev \
    libperl-dev \
    unzip \
    wget

ENV HTSLIB_VERSION=1.23
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 -P ~ && \
    tar xjvf ~/htslib-${HTSLIB_VERSION}.tar.bz2 -C ~ && \
    cd ~/htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr --enable-libcurl --enable-gcs --with-libdeflate && \
    make && make install

ENV SAMTOOLS_VERSION=${HTSLIB_VERSION}
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 -P ~ && \
    tar xjvf ~/samtools-${SAMTOOLS_VERSION}.tar.bz2 -C ~ && \
    cd ~/samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr --without-curses && \
    make && make install && \
    rm -r ~/samtools-${SAMTOOLS_VERSION}*

ENV BCFTOOLS_VERSION=${HTSLIB_VERSION}
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 -P ~ && \
    tar xjvf ~/bcftools-${BCFTOOLS_VERSION}.tar.bz2 -C ~ && \
    cd ~/bcftools-${BCFTOOLS_VERSION} && \
    ./configure --prefix=/usr --enable-libgsl --enable-perl-filters && \
    make && make install && \
    rm -r ~/bcftools-${BCFTOOLS_VERSION}*


ENV PLINK2_VERSION=20260110
RUN wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_${PLINK2_VERSION}.zip -P ~ && \
    unzip ~/plink2_linux_avx2_${PLINK2_VERSION}.zip -d /usr/bin

ENV SMR_VERSION=1.4.0
RUN wget https://yanglab.westlake.edu.cn/software/smr/download/smr-${SMR_VERSION}-linux-x86_64.zip -P ~ && \
    unzip ~/smr-${SMR_VERSION}-linux-x86_64.zip -d ~ && \
    mv ~/smr-${SMR_VERSION}-linux-x86_64/smr /usr/bin

FROM gcr.io/google.com/cloudsdktool/google-cloud-cli:slim
LABEL maintainer="mkanai@broadinstitute.org"

WORKDIR /app
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install packages - use latest R for bookworm to avoid dependency issues
RUN echo "deb https://cloud.r-project.org/bin/linux/debian bookworm-cran40/" >> /etc/apt/sources.list && \
    gpg --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' && \
    gpg --armor --export '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' | tee /etc/apt/trusted.gpg.d/cran_debian_key.asc && \
    apt-get update --fix-missing && \
    apt-get install -y \
    gawk \
    libcurl4-openssl-dev \
    libgsl-dev \
    liblapack3 \
    libopenblas-dev \
    libssl-dev \
    parallel \
    pigz \
    python3-pip \
    r-base-core \
    wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 install --break-system-packages --no-cache-dir \
    "zarr<3.0.0" \
    anndata \
    cnmf \
    cython \
    gcsfs \
    numpy \
    pandas \
    pegasusio \
    polars \
    qtl \
    scanpy[skmisc] \
    starcatpy && \
    pip3 install --break-system-packages "bgzip>=0.5.0" && \
    git clone https://github.com/mkanai/gcs_anndata && \
    cd gcs_anndata && \
    git checkout 02756044c8ea1eae1798af95c2e8617270ae9fbc && \
    pip3 install --break-system-packages . && \
    cd .. && \
    rm -rf gcs_anndata

ARG GITHUB_PAT=""
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(Rscript -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site && \
    Rscript -e "install.packages('pak')" && \
    Rscript -e "install.packages(c( \
    'argparse', \
    'data.table', \
    'dplyr', \
    'furrr', \
    'purrr', \
    'R.utils', \
    'tidyr'))" && \
    Rscript -e "pak::pkg_install(c( \
    'immunogenomics/lisi', \
    'jdstorey/qvalue', \
    'mkanai/fasthurdle@c2720845d10e72fbbb522a4880fc7acf49774f98')); pak::cache_clean()" && \
    Rscript -e "pak::pkg_install(c('ChangSuBiostats/scMultiMap', 'jaredhuling/fastglm', 'mkanai/GCSanndataR', 'scverse/anndataR', 'progress')); \
    pak::pkg_install(c('mkanai/multilinkeR@73eb51ec9ea230fcf86d075b14b6c9494c546bfe'), dependencies=FALSE); \
    pak::cache_clean()" && \
    echo "Sys.setenv(RETICULATE_PYTHON = '$(which python3)')" > ~/.Rprofile

COPY --from=builder \
    /usr/bin/bcftools \
    /usr/bin/bgzip \
    /usr/bin/plink2 \
    /usr/bin/samtools \
    /usr/bin/smr \
    /usr/bin/tabix \
    /usr/bin/

ENV BCFTOOLS_PLUGINS=/usr/libexec/bcftools/
COPY --from=builder /usr/libexec/bcftools/ /usr/libexec/bcftools/

COPY python/chunk_h5ad.py \
    python/pseudobulk.py \
    python/update_h5ad.py \
    R/compute_loco_expression_pcs.R \
    /usr/local/bin/

RUN chmod a+x \
    /usr/local/bin/*.py \
    /usr/local/bin/*.R