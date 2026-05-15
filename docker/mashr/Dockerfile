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
    libssl-dev \
    parallel \
    pigz \
    python3-pip \
    r-base-core \
    wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


ARG GITHUB_PAT=""
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
    'remotes', \
    'R.utils', \
    'tidyr'))" && \
    Rscript -e "remotes::install_version('Polychrome', version='1.5.1')" && \
    Rscript -e "remotes::install_github('willwerscheid/flashier')" && \
    echo "Sys.setenv(RETICULATE_PYTHON = '$(which python3)')" > ~/.Rprofile

RUN Rscript -e "pak::pkg_install(c( \
    'barkasn/fastSave', \
    'mkanai/mashr@86289a578b22f8faf886715fb0806e70ee8313f8'));  pak::cache_clean()"

COPY --from=builder \
    /usr/bin/bgzip \
    /usr/bin/tabix \
    /usr/bin/
