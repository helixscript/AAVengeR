#!/usr/bin/Rscript
source('./lib.R')
args = commandArgs(trailingOnly=TRUE)
f <- tmpFile()
system(paste0('wget -O ', f, ' http://bushmanlab.org/data/AAVengeR/genomeData/', args, '.tar'))
system(paste0('tar xvf ', f))
system(paste0('rsync -a ', args, '/ data/'))
system(paste0('rm -rf ', args, ' ', f))
