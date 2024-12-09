#!/usr/bin/Rscript

# AAVengeR/anchorReadStartSeqs
# John K. Everett, Ph.D.
# 
# This script determines the frequencies of anchor read sub-sequence windows.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))

# Read in the configuration file and perform basic sanity checks.
set.seed(1)
args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) stop('Expected at least one command line argument')
source(file.path(yaml::read_yaml(args[1])$softwareDir, 'lib.R'))
opt <- startModule(args)

dir.create(file.path(opt$outputDir, opt$anchorReadStartSeqs_outputDir))

reads <- tidyr::separate(readRDS(file.path(opt$outputDir, opt$anchorReadStartSeqs_inputFile)), 
                         uniqueSample, c('trial', 'subject', 'sample', 'replicate'), sep = '~')

reads$i <- group_by(reads, trial, subject, sample) %>% group_indices

cluster <- makeCluster(opt$anchorReadStartSeqs_CPUs)
clusterExport(cluster, c('opt'))

r <- bind_rows(parLapply(cluster, split(reads, reads$i), function(x){
       library(dplyr)
       z <- bind_rows(lapply(as.integer(unlist(strsplit(opt$anchorReadStartSeqs_windows, '\\s*,\\s*'))), function(n){
              seqs <- x[nchar(x$anchorReadSeq) >= n,]$anchorReadSeq
              if(length(seqs) == 0) return(tibble())
              o <- arrange(data.frame(table(substr(seqs, 1, n))), desc(Freq))
              o$n <- n
              rows <- ifelse(nrow(o) >= opt$anchorReadStartSeqs_nStartSeqs, opt$anchorReadStartSeqs_nStartSeqs, nrow(o))
              names(o) <- c('seq', 'count', 'window')
              o$freq <- o$count / length(seqs)
              o[1:rows,]
          }))
  
       z$trial   <- x$trial[1]
       z$subject <- x$subject[1]
       z$sample  <- x$sample[1]
       select(z, trial, subject, sample, window, count, freq, seq) %>% arrange(trial, subject, sample, window, desc(freq))
    }))

stopCluster(cluster)
saveRDS(r, file.path(opt$outputDir, opt$anchorReadStartSeqs_outputDir, 'startSeqs.rds'))

q()
