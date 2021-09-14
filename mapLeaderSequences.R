library(dplyr)
library(yaml)
library(parallel)

opt <- read_yaml('config.yml')
source(file.path(opt$softwareDir, 'lib.R'))

frags <- readRDS(file.path(opt$outputDir, 'standardizeFragments', 'frags.rds'))

samples <- readr::read_tsv(opt$sampleConfigFile)
samples$uniqueSample <- paste0(samples$subject, '~', samples$sample, '~', samples$replicate)

samples <- subset(samples, uniqueSample %in% frags$uniqueSample)

frags$leaderSeqTestSeqs <- samples[match(frags$uniqueSample, samples$uniqueSample),]$anchorRead.identification
