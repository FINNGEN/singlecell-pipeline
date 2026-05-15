# modified from https://github.com/satijalab/seurat-docker/blob/master/latest/Dockerfile
FROM rocker/r-ver:4.4.3
LABEL maintainer="mkanai@broadinstitute.org"

WORKDIR /app
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Set global R options
ENV RETICULATE_MINICONDA_ENABLED=FALSE

# Install Seurat's system dependencies
RUN apt-get update --fix-missing && \
    apt-get install -y \
    curl \
    git \
    libboost-all-dev \
    libcurl4-openssl-dev \
    libfftw3-dev \
    libgeos-dev \
    libgsl-dev \
    libhdf5-dev \
    libopenblas-dev \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
    llvm-14 \
    openjdk-8-jdk \
    pkg-config \
    python3-dev \
    python3-pip \
    wget && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update -y && apt-get install -y google-cloud-cli && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install UMAP
RUN LLVM_CONFIG=/usr/lib/llvm-14/bin/llvm-config pip3 install --break-system-packages llvmlite==0.42.0 && \
    pip3 install --break-system-packages \
    anndata \
    bgzip \
    gcsfs \
    numpy \
    scanpy \
    umap-learn && \
    git clone https://github.com/mkanai/gcs_anndata && \
    cd gcs_anndata && \
    git checkout 0d0c2752ef372aeac0da97fca3e88b306df06d60 && \
    pip3 install --break-system-packages . && \
    cd .. && \
    rm -rf gcs_anndata

# Install FIt-SNE
RUN git clone --branch v1.2.1 https://github.com/KlugerLab/FIt-SNE && \
    cd FIt-SNE && \
    g++ -std=c++11 -O3 src/sptree.cpp src/tsne.cpp src/nbodyfft.cpp -o bin/fast_tsne -pthread -lfftw3 -lm -Wno-address-of-packed-member

# Install dependencies & suggests
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(Rscript -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site && \
    Rscript -e "install.packages('BiocManager')" && \
    Rscript -e "BiocManager::install(c('multtest', 'S4Vectors', 'SummarizedExperiment', 'SingleCellExperiment', 'MAST', 'DESeq2', 'BiocGenerics', 'GenomicRanges', 'IRanges', 'rtracklayer', 'monocle', 'Biobase', 'limma', 'glmGamPoi'))" && \
    Rscript -e "install.packages(c('VGAM', 'R.utils', 'metap', 'Rfast2', 'ape', 'enrichR', 'mixtools'))" && \
    Rscript -e "install.packages(c('spatstat.explore', 'spatstat.geom'))" && \
    Rscript -e "install.packages('hdf5r')" && \
    Rscript -e "install.packages('Matrix')" && \
    Rscript -e "install.packages('rgeos')" && \
    Rscript -e "install.packages('remotes')" && \
    Rscript -e "install.packages('argparse')"

# Install Seurat
ENV SEURAT_VERSION=5.2.1
# Rscript -e "remotes::install_version('Seurat', version='${SEURAT_VERSION}')" && \
RUN Rscript -e "remotes::install_github('mkanai/sctransform@5c860a49fa94703be5072ed0ecc5a3c35bd5bc8e')" && \
    Rscript -e "remotes::install_github('satijalab/seurat@9354a78887e66a3f7d9ba6b726aa44123ad2d4af')" && \
    Rscript -e "remotes::install_github('mojaveazure/seurat-disk')" && \
    Rscript -e "remotes::install_github('bnprks/BPCells/r')" && \
    Rscript -e "remotes::install_github('mkanai/GCSanndataR')" && \
    Rscript -e "library(Seurat)" && \
    Rscript -e "library(SeuratDisk)" && \
    Rscript -e "library(BPCells)"

CMD [ "/bin/bash" ]
