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

ENV HTSLIB_VERSION=1.22
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 -P ~ && \
    tar xjvf ~/htslib-${HTSLIB_VERSION}.tar.bz2 -C ~ && \
    cd ~/htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr --enable-libcurl --enable-gcs --with-libdeflate && \
    make && make install

ENV LDSTORE_VERSION=v2.0
RUN wget --quiet http://www.christianbenner.com/ldstore_${LDSTORE_VERSION}_x86_64.tgz -P ~ && \
    tar xf ~/ldstore_${LDSTORE_VERSION}_x86_64.tgz -C ~ && \
    cp ~/ldstore_${LDSTORE_VERSION}_x86_64/ldstore_${LDSTORE_VERSION}_x86_64 /usr/bin/ldstore


FROM gcr.io/google.com/cloudsdktool/google-cloud-cli:slim
LABEL maintainer="mkanai@broadinstitute.org"

WORKDIR /app
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install packages
RUN apt-get update --fix-missing && \
    apt-get install -y \
    cmake \
    lbzip2 \
    libcurl4-openssl-dev \
    libgsl-dev \
    liblapack3 \
    libopenblas-dev \
    libsqlite3-dev \
    parallel \
    pigz \
    python3-pip \
    r-base-core \
    wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


RUN pip3 install --break-system-packages  \
    "bgzip>=0.5.0" \
    gcsfs \
    numpy \
    pandas

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
    Rscript -e "pak::pkg_install('stephenslab/susieR@v0.14.2'); pak::cache_clean()" && \
    echo "Sys.setenv(RETICULATE_PYTHON = '$(which python3)')" > ~/.Rprofile

RUN Rscript -e "pak::pkg_install('mkanai/rbcor@d2be97234f633def9dcfe304fc5430e3a8caf160')" && \
    pip3 install --break-system-packages \
    git+https://github.com/mkanai/ldcov@v0.5.0

COPY ./R/run_susieR.R /usr/local/bin/
RUN chmod +x /usr/local/bin/run_susieR.R

COPY --from=builder \
    /usr/bin/bgzip \
    /usr/bin/ldstore \
    /usr/bin/tabix \
    /usr/bin/
