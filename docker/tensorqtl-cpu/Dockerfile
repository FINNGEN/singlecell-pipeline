FROM gcr.io/google.com/cloudsdktool/google-cloud-cli:slim
LABEL maintainer="mkanai@broadinstitute.org"

WORKDIR /app
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install packages
RUN apt-get update --fix-missing && \
    apt-get install -y \
    libcurl4-openssl-dev \
    liblapack3 \
    libopenblas-dev \
    python3-pip \
    r-base-core && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --break-system-packages \
    "bgzip>=0.5.0" \
    cython \
    gcsfs \
    numpy \
    pandas \
    Pgenlib \
    polars \
    rpy2 && \
    pip install --break-system-packages git+https://github.com/mkanai/pycct && \
    pip install --break-system-packages \
    torch --extra-index-url https://download.pytorch.org/whl/cpu

# R
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(Rscript -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site && \
    Rscript -e "install.packages('pak')" && \
    Rscript -e "pak::pkg_install(c( \
    'jdstorey/qvalue'))" && \
    Rscript -e "pak::cache_clean()"

# ENV TENSORQTL_VERSION pval_indiv_perm_mk
# RUN git clone https://github.com/mkanai/tensorqtl && \
ENV TENSORQTL_VERSION=0c4db65a0cdc47f3b824ae530b89d270ef5e0096
RUN git clone https://github.com/broadinstitute/tensorqtl && \
    cd tensorqtl && \
    git checkout ${TENSORQTL_VERSION} && \
    pip install --break-system-packages pandas-plink && \
    pip install --break-system-packages . && \
    rm -r /app/tensorqtl

COPY python/run_tensorqtl.py \
    python/calculate_qvalues.py \
    /usr/local/bin/
RUN chmod +x \
    /usr/local/bin/run_tensorqtl.py \
    /usr/local/bin/calculate_qvalues.py