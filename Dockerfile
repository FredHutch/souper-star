FROM r-base

RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site

ENV R_LIBS_USER=/usr/local/lib/R
ENV RETICULATE_MINICONDA_ENABLED=FALSE

RUN apt-get update -qq && \
    apt-get install -y -qq --no-install-recommends \
    gtk-doc-tools \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libgsl-dev \
    libharfbuzz-dev \
    libhdf5-dev \
    libjpeg-dev \
    libmpfr-dev \
    libpng-dev \
    libssl-dev \
    libtiff5-dev \
    libxml2-dev \
    libxt-dev \
    libmagick++-dev \
    libgeos-dev \
    meson \
    bzip2 \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    wget \
    libncurses5-dev \
    liblzma-dev \
    libbz2-dev \
    pigz \
    libcurl4-openssl-dev \
    python3-pip \
    python3-dev \
    python3 \
    awscli \
    jq && \
    rm -rf /var/lib/apt/lists/* \
    && cd /usr/local/bin \
    && ln -s /usr/bin/python3 python \
    && pip3 --no-cache-dir install --upgrade pip

ENV USER=shareseq
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

RUN R --no-echo --no-restore --no-save -e "install.packages(c('devtools','hdf5r','IRkernel','BiocManager','Cairo','magick'))" \
    && R --no-echo --no-restore --no-save -e "BiocManager::install(c('GenomeInfoDbData','GenomicRanges','Rsamtools'), update=F, ask=F)" \
    && R --no-echo --no-restore --no-save -e "devtools::install_github('GreenleafLab/ArchR@v1.0.1', repos = BiocManager::repositories());ArchR::installExtraPackages()" \
    && R --no-echo --no-restore --no-save -e "devtools::install_github('immunogenomics/presto')" \
    && R --no-echo --no-restore --no-save -e "remotes::install_version('Seurat', version = '4.1.1')" \
    && R --no-echo --no-restore --no-save -e "install.packages(c('logr','hexbin'))"

ENV DEBIAN_FRONTEND=noninteractive

RUN python3 -m pip install \
    jupyter \
    papermill \
    numpy \
    pandas \
    pysam \
    simplesam \
    scikit-learn \
    matplotlib \
    seaborn \
    MACS2

ENV BWA_VERSION 0.7.17
ENV SAMTOOLS_VERSION 1.16.1
ENV HTSLIB_VERSION 1.16

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

RUN cd /opt/ \
    && wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    && tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 \
    && rm -rf htslib-${HTSLIB_VERSION}.tar.bz2  \
    && cd htslib-${HTSLIB_VERSION}/ \
    && make && make install

RUN mkdir -p /opt/subset-bam \
    && cd /opt/subset-bam/ \
    && wget https://github.com/10XGenomics/subset-bam/releases/download/v1.1.0/subset-bam_linux \
    && mv subset-bam_linux subset-bam \
    && chmod +x subset-bam

ENV PATH="/opt/bwa-${BWA_VERSION}/:/opt/samtools-${SAMTOOLS_VERSION}/:/opt/htslib-${HTSLIB_VERSION}/:/opt/subset-bam/:${PATH}"