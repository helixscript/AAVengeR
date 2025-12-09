#!/usr/bin/Rscript
options(scipen = 999, useFancyQuotes = FALSE) 

# AAVengeR/anchorReadStartSeqs
# John K. Everett, Ph.D.
# 
# This script determines the frequencies of anchor read sub-sequence windows.

for (p in c('parallel', 'dplyr')) suppressPackageStartupMessages(library(p, character.only = TRUE))

# Read in the configuration file and perform basic sanity checks.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) stop('Expected at least one command line argument')
source(file.path(yaml::read_yaml(args[1])$softwareDir, 'lib', 'main.R'))
opt <- startModule(args)

dir.create(file.path(opt$outputDir, opt$anchorReadStartSeqs_outputDir), showWarnings = FALSE)

# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, opt$anchorReadStartSeqs_outputDir, 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)
write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)

updateLog('Starting anchorReadStartSeqs.')
updateMasterLog()

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
updateLog('anchorReadStartSeqs completed.')
updateMasterLog()
closeAllConnections()

q()
