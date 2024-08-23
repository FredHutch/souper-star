from ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq \
  && apt-get install -qq bzip2 gcc g++ make zlib1g-dev wget libncurses5-dev liblzma-dev libbz2-dev pigz libcurl4-openssl-dev \
  && apt-get install -y python3-pip python3-dev python3.8 awscli jq tabix\
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3.8 python \
  && pip3 --no-cache-dir install --upgrade pip

RUN python3 -m pip install numpy pandas pysam simplesam scikit-learn matplotlib seaborn
RUN python --version

ENV BWA_VERSION=0.7.18
ENV SAMTOOLS_VERSION=1.16.1
ENV BEDTOOLS_VERSION=2.31.1

RUN cd /opt/ \
    && wget https://github.com/lh3/bwa/archive/refs/tags/v${BWA_VERSION}.tar.gz \
    && tar -xzf v${BWA_VERSION}.tar.gz \
    && rm -f v${BWA_VERSION}.tar.gz \
    && cd /opt/bwa-${BWA_VERSION}/ \
    && make

RUN cd /opt/ \
    && wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && rm -rf samtools-${SAMTOOLS_VERSION}.tar.bz2  \
    && cd samtools-${SAMTOOLS_VERSION}/ \
    && make && make install

RUN cd /opt/ \
    && wget https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz \
    && tar -xzf bedtools-${BEDTOOLS_VERSION}.tar.gz \
    && rm -rf bedtools-${BEDTOOLS_VERSION}.tar.gz \
    && cd bedtools2 \
    && make

ENV PATH=${PATH}:/opt/bedtools2/bin
RUN bedtools --version

RUN mkdir -p /opt/subset-bam \
    && cd /opt/subset-bam/ \
    && wget https://github.com/10XGenomics/subset-bam/releases/download/v1.1.0/subset-bam_linux \
    && mv subset-bam_linux subset-bam \
    && chmod +x subset-bam

ENV PATH="/opt/bwa-${BWA_VERSION}/:/opt/samtools-${SAMTOOLS_VERSION}/:/opt/subset-bam/:${PATH}"