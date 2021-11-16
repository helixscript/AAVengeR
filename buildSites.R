library(dplyr)

opt <- yaml::read_yaml('config.yml')
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$buildSites_outputDir))

frags <- readRDS(file.path(opt$outputDir, opt$buildSites_inputFile))

if('buildSites_excludeSites' %in% names(opt)){
  p <- unlist(lapply(strsplit(unlist(strsplit(opt$buildSites_excludeSites, '\\|')), ','), function(x){
    if(length(x) != 3) return()
    subset(frags, seqnames == x[1] & start >= x[2] & start <= x[3])$posid
  }))
  
  frags <- subset(frags, ! posid %in% p)
}

samples <- loadSamples()

frags$fragWidth <- frags$fragEnd - frags$fragStart + 1
frags$posid <- paste0(frags$chromosome, frags$strand, ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd))


sites <- bind_rows(lapply(split(frags, paste(frags$trial, frags$subject, frags$sample, frags$posid)), function(x){
  if(nrow(x) > 1){
    i <- rep(TRUE, nrow(x))
    
    # Check that leader sequences are all similar now that we are combining fragments from different replicates.
    r <- representativeSeq(x$repLeaderSeq)
    
    if(r[[1]] > opt$buildStdFragments_maxLeaderSeqDiffScore){
      # There is a conflict, one or more fragments have a markedly different leader sequences then the other fragments.
      # Attempt to salvage this site by retaining the majority of fragments with similiar leader sequences.
      
      i <- as.vector(stringdist::stringdistmatrix(x$repLeaderSeq, r[[2]]) / nchar(x$repLeaderSeq) < opt$buildFragments_maxLeaderSeqDiffScore)
      if(sum(i)/nrow(x) >= opt$buildSites_assemblyConflictResolution){
        x <- x[i,]
      } else {
        return(tibble())
      }
    }
    
    r <- representativeSeq(x$repLeaderSeq)
    
    if('flags' %in% names(samples)){
      x$flags <- paste0(unique(subset(samples, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1])$flags), collapse = ',')
    } else {
      x$flags <- NA
    }
    
    return(dplyr::mutate(x, estAbund = n_distinct(fragWidth), position = ifelse(strand[1] == '+', fragStart[1], fragEnd[1]), 
                  reads = sum(reads), repLeaderSeq = r[[2]], fragmentsRemoved = sum(!i)) %>%
           dplyr::select(subject, sample, chromosome, strand, position, posid, estAbund, reads, fragmentsRemoved, repLeaderSeq, flags) %>%
           dplyr::slice(1))
    }else{
    
      if('flags' %in% names(samples)){
        x$flags <- paste0(unique(subset(samples, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1])$flags), collapse = ',')
      } else {
        x$flags <- NA
      }
      
      return(dplyr::mutate(x, estAbund = n_distinct(fragWidth), position = ifelse(strand[1] == '+', fragStart[1], fragEnd[1]), fragmentsRemoved = 0) %>%
             dplyr::select(subject, sample, chromosome, strand, position, posid, estAbund, reads, fragmentsRemoved, repLeaderSeq, flags))
  }
}))


saveRDS(sites, file.path(opt$outputDir, opt$buildSites_outputDir, opt$buildSites_outputFile))

