

# Download base docker with r-base already included

https://hub.docker.com/_/r-base/

```sh
sudo docker pull r-base
```

Check for succesfull install

Download base docker image with R

Test the base image works

This is for learning

```sh
sudo docker run -v /data:/ -ti --rm r-base

sudo docker run -ti --rm r-base
```

# Build image

If the dockerfile and the R and conda install scripts are all set up, then build the image with the code below.

```sh
cd /data/AAVengeR/docker

sudo docker build --no-cache -t aavenger_docker_v43 . 

# check install using interactive mode
sudo docker run -it --mount type=bind,source=/home/ubuntu/,dst=/home/ubuntu aavenger_docker bash
```

# Check aavenger_base passses the AAVengeR unit test

Code to test docker image/container. Runs through the AAvengeR unit test. Necessary to update bashrc with paths to bioinformatic dependencies and to provide an alias for python2. 

```sh
sudo docker run -it --mount type=bind,source=/home/ubuntu/,dst=/home/ubuntu aavenger_docker bash

echo "export PATH=$PATH:/opt/conda/envs/bioinformatic_dependencies/bin:/opt/conda/envs/python2_env/bin" >> /root/.bashrc

echo "alias python2=/opt/conda/envs/python2_env/bin/python" >> /root/.bashrc 

source /root/.bashrc

source activate bioinformatic_dependencies

Rscript /home/mount/AAVengeR/data/aavenger_unit_test/run_aavenger_unit_test.R build /home/mount/AAVengeR /home/mount/aavenger_unit_test_results2

Rscript /home/mount/AAVengeR/data/aavenger_unit_test/run_aavenger_unit_test.R test /home/mount/AAVengeR /home/mount/aavenger_unit_test_results2

Rscript aavenger.R /home/data/test_synthetic/config.yml

```

# Manual installation/checking

All sections below go into creating the necessary scripts and configuration files so that the the image can be built from the dockerfile. 


# System level depencencies

All system-level dependencies are installed in this chunk. In future apt-get update needs to be removed for reproducibility.

```sh
apt-get update
apt-get -y install libcurl4-openssl-dev
apt-get -y install libssl-dev
apt-get -y install libxml2-dev
apt-get -y install libgmp3-dev
apt-get -y install libmariadb-dev
#apt-get -y install mariadb-server mariadb-client # untested
```

# Python and bioinformatic dependencies

## Making environment files manually

This chunk manually installs each package in environments bioinformatic_dependencies or python2_env. 

Afterwards a yml config file is created for each.

```sh
# run in interactive
# <!-- sudo docker run -it --mount type=bind,source=/home/ubuntu/,dst=/home/data aavenger bash -->
# sudo docker run -it --mount type=bind,source=/home/ubuntu/,dst=/home/mount aavenger bash #CORRECT 
sudo docker run -it --mount type=bind,source=/data,dst=/data/ aavenger bash #CORRECT 

conda install -c conda-forge mamba -y

conda create -n bioinformatic_dependencies python=3.8 -y

source activate bioinformatic_dependencies

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# mamba install -c bioconda cutadapt=3.5 -y

pip install cutadapt==3.5

# mamba install cutadapt -y

# pip install cutadapt==3.5

# mamba install -c bioconda BLAT=36 -y 

# mamba install -c bioconda BLAT=35 -y 

mamba install -c bioconda muscle=3.8.31 -y

mamba install -c bioconda blast=2.12.0 -y

mamba install -c bioconda hmmer=3.3.2 -y

mamba install -c bioconda mafft=7.310 -y

RUN mamba install -c bioconda bwa-mem2=2.2.1 -y 

# # 4.8.1
RUN mamba install -c bioconda cd-hit=4.8.1 -y

conda deactivate

conda create -n python2_env python=2.7 -y

source activate python2_env

conda install -c bioconda biopython

pip install numpy

conda deactivate


source activate bioinformatic_dependencies

conda env export > bioinformatic_dependencies.yml

head -n 1000 bioinformatic_dependencies.yml

conda deactivate

source activate python2_env

conda env export > python2_env.yml

head -n 1000 python2_env.yml


conda deactivate


conda create -n blat_dependency

source activate blat_dependency

mamba install -c conda-forge libuuid=2.38.1  -y

mamba install -c conda-forge libstdcxx-ng=13.2.0 -y

mamba install -c conda-forge libzlib=1.2.13 -y

mamba install -c conda-forge zlib=1.3 -y

mamba install -c bioconda ucsc-blat=445 -y

conda env export > blat_dependency.yml

```

## Install packages from environment yml

If yml files are already available, this chunk installs package from those environment files.


```sh
conda env create -f python2_env.yml

conda env create -f bioinformatic_dependencies.yml
```

# Set up mariaDB connection

Code necessary for mariadb to be set up in the docker container. Necessary if you want to connect AAVengeR to SQL database of integration sites.


```sh
mariadb -u admin -p -h 174.129.238.44 AAVengeR
#[AAVengeR]
#user=admin
#password=iAmAdmin+1
#host=174.129.238.44
#port=3306
#database=AAVengeR
touch ~/.my.cnf 
vim ~/.my.cnf
# copy and paste above hashed info
:wq
```

# R packages dependencies

All R packages that are required by AAVengeR are here. After 

```R
R # enter R

as.data.frame(installed.packages())

install.packages('remotes',repos='https://cloud.r-project.org/')
library(remotes)
install.packages('devtools',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages("lubridate", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("data.table", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("ggplot2", repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages("scales", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("tidyverse", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("RColorBrewer", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("gridExtra", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("png", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("igraph", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("RMySQL", repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages("readr", repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('BiocManager', repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages('RcppAlgos', repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages('stringdist', repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages('openxlsx', repos='https://cloud.r-project.org/', quiet=FALSE)

BiocManager::install(c("GenomicRanges", "ShortRead",'GenomicRanges','Biostrings'), update = FALSE)

install.packages("RMariaDB", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages('pander', repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages('dtplyr', repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages('argparse', repos='https://cloud.r-project.org/', quiet=FALSE)


q() #quit R

Rscript R_dependencies.R # test that all packages were installed successfully
```









