library(dplyr)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

frags <- readRDS(file.path(opt$outputDir, opt$buildReplicateSiteSummary_inputFile))

sites$replicate <- as.integer(sites$replicate)

minReplicate <- min(sites$replicate)
maxReplicate <- max(sites$replicate)

lapply(split(sites, paste(sites$trial, sites$subject, sites$subject, sites$posid)), function(x){
  if(n_distinct(x$replicate) == 3) browser()

  return()
  
  r <- representativeSeq(x$repLeaderSeq)
  
  # Skip replicates that have a distincatly different LTR/ITR remant to the concensus sequences.
  replicatesToSkip <- NA
  if(r[[1]] > opt$buildStdFragments_maxLeaderSeqDiffScore){
    i <- as.vector(stringdist::stringdistmatrix(x$repLeaderSeq, r[[2]]) / nchar(x$repLeaderSeq) <= opt$buildStdFragments_maxLeaderSeqDiffScore)
    replicatesToSkip <- x[! i,]$replicate
  }
  
  bind_cols(lapply(minReplicate:maxReplicate, function(r){
    o <- subset(x, replicate == r)
        
    if(nrow(o) == 1){
      t <- tibble(x1 = o$estAbund, x2 = o$linkerMolCodes, x3 = o$reads, x4 = o$repLeaderSeq)
    } else {
      t <- tibble(x1 = NA, x2 = NA, x3 = NA, x4 = NA)
    }
    
    names(t) <- c(paste0('rep', o$replicate[1], '-estAbund'), 
                  paste0('rep', o$replicate[1], '-linkerMolCodes'), 
                  paste0('rep', o$replicate[1], '-reads'), 
                  paste0('rep', o$replicate[1], '-repLeaderSeq'))
    t
  }))
  
})