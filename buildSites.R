library(dplyr)
library(yaml)
library(parallel)

sites <- readRDS('/home/everett/projects/AAVengeR/data/testData/sites.rds')

opt <- read_yaml('config.yml')
source(file.path(opt$softwareDir, 'lib.R'))

frags <- readRDS(file.path(opt$outputDir, 'mapLeaderSequences', 'frags.rds'))

frags$posid <- paste0(frags$seqnames, frags$strand, ifelse(frags$strand == '+', frags$start, frags$end))

# # Collapse replicate replicate fragments to sample replicates.
# o <- lapply(split(frags, paste(frags$subject, frags$sample, frags$seqnames, frags$start, frags$end)), function(x){
#    if(nrow(x) > 1)browser()
# })