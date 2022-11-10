from ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq \
  && apt-get install -qq bzip2 gcc g++ make zlib1g-dev wget libncurses5-dev liblzma-dev libbz2-dev pigz libcurl4-openssl-dev \
  && apt-get install -y python3-pip python3-dev python3.8 awscli jq \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3.8 python \
  && pip3 --no-cache-dir install --upgrade pip

RUN python3 -m pip install numpy pandas pysam simplesam
RUN python --version

ENV BWA_VERSION 0.7.17
ENV SAMTOOLS_VERSION 1.16.1
ENV SOUPERCELL_VERSION 2.0

RUN cd /opt/ \
    && wget https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2 \
    && tar -xjf bwa-${BWA_VERSION}.tar.bz2 \
    && rm -f bwa-${BWA_VERSION}.tar.bz2 \
    && cd /opt/bwa-${BWA_VERSION}/ \
    && make

RUN cd /opt/ \
    && wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && rm -rf samtools-${SAMTOOLS_VERSION}.tar.bz2  \
    && cd samtools-${SAMTOOLS_VERSION}/ \
    && make && make install

RUN mkdir -p /opt/subset-bam \
    && cd /opt/subset-bam/ \
    && wget https://github.com/10XGenomics/subset-bam/releases/download/v1.1.0/subset-bam_linux \
    && mv subset-bam_linux subset-bam \
    && chmod +x subset-bam

RUN mkdir -p /opt/conda \
    && cd /opt/conda \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh \
    && bash Miniconda3-py38_4.12.0-Linux-x86_64.sh

RUN mkdir -p /opt/soupercell \
    && cd /opt/soupercell \
    && wget https://github.com/wheaton5/souporcell/archive/refs/tags/${SOUPERCELL_VERSION}.tar.gz \
    && tar xvf soupercell-${SOUPERCELL_VERSION}.tar.gz \
    && conda env create -f /opt/soupercell/soupercell-${SOUPERCELL_VERSION}/souporcell_env.yaml \
    && conda activate souporcell \
    && curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh \
    && cd /path/to/souporcell/souporcell && cargo build --release \
    && cd /path/to/souporcell/troublet && cargo build --release


ENV PATH="/opt/bwa-${BWA_VERSION}/:/opt/samtools-${SAMTOOLS_VERSION}/:/opt/subset-bam/:${PATH}"