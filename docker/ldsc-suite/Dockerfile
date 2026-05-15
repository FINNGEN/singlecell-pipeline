FROM gcr.io/google.com/cloudsdktool/google-cloud-cli:slim AS builder

RUN apt-get update --fix-missing && \
    apt-get install -y \
    unzip \
    wget

ENV PLINK_VERSION=20241022
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_${PLINK_VERSION}.zip -P ~ && \
    unzip ~/plink_linux_x86_64_${PLINK_VERSION}.zip -d /usr/bin

ENV PLINK2_VERSION=20250420
RUN wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_${PLINK2_VERSION}.zip -P ~ && \
    unzip ~/plink2_linux_avx2_${PLINK2_VERSION}.zip -d /usr/bin

FROM gcr.io/google.com/cloudsdktool/google-cloud-cli:slim AS pyenv_builder
WORKDIR /root

RUN apt-get update --fix-missing && \
    apt-get install -y \
    libbz2-dev \
    libffi-dev \
    liblzma-dev \
    libreadline-dev \
    libsqlite3-dev \
    libssl-dev

ENV PYENV_ROOT=/root/.pyenv
ENV PATH="$PYENV_ROOT/bin:$PATH"
RUN git clone https://github.com/pyenv/pyenv ~/.pyenv && \
    pyenv install 2.7.18 && \
    pyenv global 2.7.18 && \
    (pyenv init >> ~/.bashrc 2>&1 || true) && \
    (pyenv init >> ~/.bash_profile 2>&1 || true) && \
    (pyenv init >> ~/.profile 2>&1 || true)

SHELL ["/bin/bash", "-c"]
RUN eval "$(pyenv init - bash)" && \
    git clone https://github.com/mkanai/ldsc && \
    cd ldsc && \
    git checkout 02426dee04b98f1b90dc6a47391174a3736863c9 && \
    pip install -r requirements.txt

RUN git clone https://github.com/mkanai/mesc && \
    cd mesc && \
    git checkout de79fdcbce857556ef2dd76c4fd9e5f906499285

FROM gcr.io/google.com/cloudsdktool/google-cloud-cli:slim
LABEL maintainer="mkanai@broadinstitute.org"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update --fix-missing && \
    apt-get install -y \
    bedtools \
    libbz2-dev \
    libcurl4-openssl-dev \
    libffi-dev \
    liblapack3 \
    liblzma-dev \
    libopenblas-dev \
    libreadline-dev \
    libsqlite3-dev \
    libssl-dev \
    parallel && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY --from=pyenv_builder /root /root
COPY --from=builder /usr/bin/plink* /usr/bin/
