
# Overview

Dockerfile, blat_dependency.yml, bionformatic_dependencies.yml, python2_env.yml, R_dependencies.R, and check_R_library_installation.R are updated in this file and then used to build the latest version of the AAVengeR docker image

# Latest

```sh
cd /data/AAVengeR/docker

sudo docker run -it --mount type=bind,source=/home/ubuntu/,dst=/home/ubuntu aavenger_docker_v2_ bash
```

## aavenger_docker_v2_2

* Alex McFarland

* 08/26/2024

### New system libraries

In dockerfile (except apt-get -y update) :

```
apt-get -y update

apt-get -y install mariadb-server # mariadb # doesn't work at the moment

apt-get -y install procps # free command

apt-get -y install nano

apt-get -y install htop

apt-get -y install ssh

apt-get -y install rsync
```

### New bioinformatic libraries

**PEAR**

```sh
source activate bioinformatic_dependencies

mamba install bioconda::pear
```

In bionformatic_dependencies.yml :

```
- pear=0.9.6=h9d449c0_10
```


### Notes

These are the system commands that need to be in the default search path:

```
blat
cutadapt    (3.5 or greater)
makeblastdb  (ncbi-blast-2.12.0+ or greater)
blastn                 (ncbi-blast-2.12.0+ or greater)
xz
unxz
wget
gzip
gunzip
tar
bwa-mem2
pear                 (0.9.11 or greater)
faToTwoBit
cd-hit-est     (4.8.1 or greater)
free
which
python
python2
find
rsync
```

```R
These are the R libraries currently being used:
> unique(gsub('\\(|\\)', '', c(gsub('library\\(', '', stringr::str_extract(system('find ./ -name "*.R" | xargs grep library', intern = TRUE), 'library\\(\\S+')),
+        gsub('::$', '', stringr::str_extract(system('find ./ -name "*.R" | xargs grep ::', intern = TRUE), '[a-zA-Z]+::')))))
 
   "yaml"          "dplyr"         "lubridate"     "pander"        "parallel"      "data.table"    "Biostrings"    "RMariaDB"      "ShortRead"     "dtplyr"        
   "stringdist"    "remotes"       "readr"         "GenomicRanges"
   "igraph"        "ggplot2"       "stringr"       "scales"        "gridExtra"     "png"           "RcppAlgos"     "openxlsx"      "rtracklayer"   "tidyr"         
   "optparse"      "RMySQL"        "IRanges"       "argparse"     
   "testthat"      "base"          "table"         "stringi"       "BiocManager"      "Vectors"       "Matrix"        "grDevices"     "path"          "tibble"        
   "grid" 
```

### Build

```sh
cd /data/AAVengeR/docker

sudo docker build --no-cache -t aavenger_docker_v2_2 . 

sudo docker save aavenger_docker_v2_2:latest | gzip > aavenger_docker_v2_2.tar.gz

scp ./aavenger_docker_v2_2.tar.gz agmcfarland@microb120.med.upenn.edu:/home/agmcfarland

agmcfarland@microb191.med.upenn.edu:/home/agmcfarland/flu_project/scripts


rm aavenger_docker_v2_2.tar.gz
```

## aavenger_docker_v2_1

* Alex McFarland

* 04/18/2024

* Enforced pulling from r-base:4.3.2

* Added this.path R library to R_dependencies.R

* Added optparse R library to R_dependencies.R

* Fixed iGraph R library not installing

* Updated check_r_library_installation.R with new libraries

```sh
cd /data/AAVengeR/docker

sudo docker build --no-cache -t aavenger_docker_v2_1 . 
```

## aavenger_docker_v2

2024/02/23

* Added RMariaDb R library and system libraries

* Added rtracklayer R library via BiocManager

* Added faToTwoBit via conda

```sh
/data/AAVengeR/docker

sudo docker build --no-cache -t aavenger_docker_v2 . 
```

## aavenger_docker_v43

* Added Blat v37

```sh
/data/AAVengeR/docker

sudo docker build --no-cache -t aavenger_docker_v43 . 
```