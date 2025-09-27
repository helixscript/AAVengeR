FROM ubuntu:24.04

ENV R_VERSION=4.5.1 \
    CUTADAPT_VERSION=5.1 \
    DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq && apt-get -y install --no-install-recommends \
    ca-certificates \
    build-essential \
    default-jre \
    gfortran \
    libreadline-dev \
    xorg-dev \
    libbz2-dev \
    liblzma-dev \
    curl \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libsodium-dev \
    zlib1g-dev \
    wget \
    git \
    python3-pandas \
    python3-biopython \
    python3-yaml \
    mariadb-server \
    procps \
    nano \
    htop \
    libfribidi-dev \
    libharfbuzz-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgmp-dev \
    libgmp3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* /tmp/* /var/tmp/*

RUN cd / \
    && wget https://github.com/samtools/samtools/releases/download/1.22.1/samtools-1.22.1.tar.bz2 \
    && tar xvf samtools-1.22.1.tar.bz2 \
    && cd samtools-1.22.1 \
    && ./configure --prefix=/usr/bin \
    && make \
    && make install \
    && cd / \
    && rm -rf samtools*

RUN cd / \
    && wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 \
    && tar xvf bcftools-1.22.tar.bz2 \
    && cd bcftools-1.22 \
    && ./configure --prefix=/usr/bin \
    && make \
    && make install \
    && cd / \
    && rm -rf bcftools*

RUN cd / \
    && wget https://github.com/samtools/htslib/releases/download/1.22.1/htslib-1.22.1.tar.bz2 \
    && tar xvf htslib-1.22.1.tar.bz2 \
    && cd htslib-1.22.1 \
    && ./configure --prefix=/usr/bin \
    && make \
    && make install \
    && cd / \
    && rm -rf htslib*

RUN cd / \
    && wget -O o.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-x64-linux.tar.gz \
    && tar xvf o.tar.gz \
    && mv ncbi-blast-2.17.0+/bin/* /usr/bin/ \
    && rm o.tar.gz \
    && rm -rf ncbi-blast*

RUN cd / \
    && wget http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz \
    && tar xvf hmmer-3.4.tar.gz \
    && cd hmmer-3.4 \
    && ./configure \
    && make \
    && make install \
    && cd / \
    && rm -rf hmmer-3.4 \
    && rm hmmer-3.4.tar.gz

RUN cd / \
    && wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz \
    && tar xvf cd-hit-v4.8.1-2019-0228.tar.gz \
    && cd cd-hit-v4.8.1-2019-0228 \
    && make \
    && mv cd-hit cd-hit-est /usr/bin/ \
    && cd / \
    && rm -rf cd-hit-v4.8.1-2019-0228.tar.gz cd-hit-v4.8.1-2019-0228

RUN wget "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/blat/blat" \
    && chmod 755 blat \
    && mv blat /usr/bin/

RUN wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/faToTwoBit \
    && chmod 755 faToTwoBit \
    && mv faToTwoBit /usr/bin/

RUN wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/twoBitToFa \
    && chmod 755 twoBitToFa \
    && mv twoBitToFa /usr/bin/

RUN cd / \ 
    && mkdir tempDir \
    && cd tempDir \
    && wget https://anaconda.org/bioconda/pear/0.9.6/download/linux-64/pear-0.9.6-hb1d24b7_12.tar.bz2 \
    && tar xvf pear-0.9.6-hb1d24b7_12.tar.bz2 \
    && mv bin/pear /usr/bin/pear \
    && cd / \
    && rm -rf tempDir

RUN apt-get update -qq && apt-get -y install --no-install-recommends python3 python3-venv python3-pip \
    && rm -rf /var/lib/apt/lists/* \
    && python3 -m venv /opt/cutadapt \
    && /opt/cutadapt/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/cutadapt/bin/pip install --no-cache-dir "cutadapt==${CUTADAPT_VERSION}" \
    && ln -s /opt/cutadapt/bin/cutadapt /usr/bin/cutadapt

RUN wget -c https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz \
    && tar -xf R-${R_VERSION}.tar.gz \
    && cd R-${R_VERSION} \
    && ./configure --prefix=/opt/R/${R_VERSION} --enable-memory-profiling --enable-R-shlib --with-blas --with-lapack \
    && make -j$(nproc) \
    && make install \
    && cd .. \
    && rm -rf R-${R_VERSION} R-${R_VERSION}.tar.gz

RUN ln -s /opt/R/${R_VERSION}/bin/R /usr/bin/R
RUN ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/bin/Rscript

RUN R -e "install.packages('RcppAlgos', repos='http://cran.case.edu')"
RUN R -e "install.packages('BiocManager', repos='http://cran.case.edu')"
RUN R -e "library('BiocManager'); BiocManager::install('ShortRead')"
RUN R -e "install.packages('devtools', repos='http://cran.case.edu')"
RUN R -e "install.packages('tidyverse', repos='http://cran.case.edu')"
RUN R -e "install.packages('dtplyr', repos='http://cran.case.edu')"
RUN R -e "install.packages('RColorBrewer', repos='http://cran.case.edu')"
RUN R -e "install.packages('optparse', repos='http://cran.case.edu')"
RUN R -e "install.packages('pander', repos='http://cran.case.edu')"
RUN R -e "install.packages('knitr', repos='http://cran.case.edu')"
RUN R -e "install.packages('argparse', repos='http://cran.case.edu')"
RUN R -e "install.packages('RMySQL', repos='http://cran.case.edu')"
RUN R -e "install.packages('vegan', repos='http://cran.case.edu')"
RUN R -e "library('BiocManager'); BiocManager::install('rtracklayer')"
RUN R -e "install.packages('stringdist', repos='http://cran.case.edu')"
RUN R -e "install.packages('igraph', repos='http://cran.case.edu')"
RUN R -e "install.packages('gridExtra', repos='http://cran.case.edu')"
RUN R -e "install.packages('openxlsx', repos='http://cran.case.edu')"
RUN R -e "install.packages('yaml', repos='http://cran.case.edu')"
RUN R -e "install.packages('data.table', repos='http://cran.case.edu')"
RUN R -e "install.packages('RMariaDB', repos='http://cran.case.edu')"
RUN R -e "install.packages('reldist', repos='http://cran.case.edu')"
RUN R -e "install.packages('wordcloud', repos='http://cran.case.edu')"
RUN R -e "install.packages('gtools', repos='http://cran.case.edu')"
RUN R -e "install.packages('kableExtra', repos='http://cran.case.edu')"
RUN R -e "install.packages('gridExtra', repos='http://cran.case.edu')"
RUN R -e "install.packages('RColorBrewer', repos='http://cran.case.edu')"
RUN R -e "install.packages('geneRxCluster', repos='http://cran.case.edu')"
RUN R -e "install.packages('ggforce', repos='http://cran.case.edu')"

COPY entrypoint.R /usr/local/bin/entrypoint.R
RUN chmod +x /usr/local/bin/entrypoint.R

ENTRYPOINT ["/usr/local/bin/entrypoint.R"]
