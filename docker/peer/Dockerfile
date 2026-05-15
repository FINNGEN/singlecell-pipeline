# original source: https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/Dockerfile
# Dockerfile for GTEx QTL pipeline
FROM ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive
LABEL maintainer="mkanai@broadinstitute.org"

WORKDIR /app
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        build-essential \
        curl \
        lbzip2 \
        libboost-all-dev \
        libcurl3-dev \
        libgsl-dev \
        libhdf5-serial-dev \
        openjdk-8-jdk \
        python3 \
        python3-pip \
        r-base-core \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/


# workaround for PEER, see https://github.com/mz2/peer/issues/4
RUN apt-get update \
    && apt-get install -y \
        gcc-5 \
        g++-5 \
        gfortran-5 \
        cmake \
        swig \
    && rm -rf /var/lib/apt/lists/* \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 50 --slave /usr/bin/g++ g++ /usr/bin/g++-5

# R
# https://github.com/downloads/PMBio/peer/R_peer_source_1.3.tgz is 404 now
RUN wget https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/master/qtl/R_peer_source_1.3.tgz && \
    R CMD INSTALL R_peer_source_1.3.tgz && \
    rm R_peer_source_1.3.tgz && \
    echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/rlang/rlang_1.0.6.tar.gz')" && \
    Rscript -e "install.packages(c('data.table', 'R.utils', 'findpython', 'generics', 'jsonlite'))" && \
    Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/R6/R6_2.5.1.tar.gz')" && \
    Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/cli/cli_3.6.2.tar.gz')" && \
    Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/glue/glue_1.6.2.tar.gz')" && \
    Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/lifecycle/lifecycle_1.0.3.tar.gz')" && \
    Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/argparse/argparse_2.1.5.tar.gz')" && \
    Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/vctrs/vctrs_0.5.2.tar.gz')" && \
    Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/withr/withr_2.5.2.tar.gz')" && \
    Rscript -e "install.packages(c('pillar', 'tibble', 'tidyselect'))" && \
    Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_1.0.10.tar.gz')" && \
    Rscript -e "library(argparse); library(data.table); library(dplyr)"

COPY R/compute_peer_factor.R /usr/local/bin/
RUN chmod +x /usr/local/bin/compute_peer_factor.R