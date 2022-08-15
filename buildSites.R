library(dplyr)
library(lubridate)
library(Biostrings)
library(sonicLength)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$buildSites_outputDir))

# Read in Standardized fragments.
frags <- readRDS(file.path(opt$outputDir, opt$buildSites_inputFile))

randomIDs <- tibble()
if(opt$demultiplex_captureRandomLinkerSeq){
  randomIDs <- Reduce('append', lapply(list.files(file.path(opt$outputDir, opt$demultiplex_outputDir), pattern = 'randomIDs', full.names = TRUE), readDNAStringSet))
  randomIDs <- tibble(id = names(randomIDs), seq = as.character(randomIDs))
}

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

if(opt$buildSites_level == 'replicate'){
  frags$s <- paste(frags$trial, frags$subject, frags$sample, frags$replicate, frags$repLeaderSeqGroup, frags$posid)
} else if (opt$buildSites_level == 'sample'){
  frags$s <- paste(frags$trial, frags$subject, frags$sample, frags$repLeaderSeqGroup, frags$posid)
} else if (opt$buildSites_level == 'subject'){
  frags$s <- paste(frags$trial, frags$subject, frags$repLeaderSeqGroup, frags$posid)
} else {
  stop('Error: buildSites_level not defined with replicate, sample, or subject')
}

frags$uniqueSample <- paste0(frags$trial, '~', frags$subject, '~', frags$sample, '~', frags$replicate)

sonicLengths <- tibble()

if(opt$buildSites_sonicLengthAbund){
  sonicLengths <- bind_rows(lapply(split(frags, frags$s), function(x){
    
   #message('a ', x$uniqueSample[1], ' - ', nrow(x))
   #if(x$uniqueSample[1] == 'Encoded~p2003_L~GTSP5112~4') browser()
    
    # Sonic length fails on 1 fragment.
    if(nrow(x) == 1){
      z <- list()
      z$theta <- 1
    } else {
      z <- estAbund(x$posid, x$fragWidth) 
    }
    
    message('b ', x$uniqueSample[1])
    
    tibble(trial = x$trial[1], subject = x$subject[1], sample = x$sample[1], 
          replicate = x$replicate[1], uniqueSample = x$uniqueSample[1],
          repLeaderSeqGroup = x$repLeaderSeqGroup[1],
          posid = x$posid[1], estAbund = floor(z$theta))
    }))
  
  if(opt$buildSites_level == 'replicate'){
    sonicLengths$s <- paste(sonicLengths$trial, sonicLengths$subject, sonicLengths$sample, sonicLengths$replicate, sonicLengths$repLeaderSeqGroup, sonicLengths$posid)
  } else if (opt$buildSites_level == 'sample'){
    sonicLengths$s <- paste(sonicLengths$trial, sonicLengths$subject, sonicLengths$sample, sonicLengths$repLeaderSeqGroup, sonicLengths$posid)
  } else if (opt$buildSites_level == 'subject'){
    sonicLengths$s <- paste(sonicLengths$trial, sonicLengths$subject, sonicLengths$repLeaderSeqGroup, sonicLengths$posid)
  } else {
    stop('Error: buildSites_level not defined with replicate, sample, or subject')
  }
}

sites <- bind_rows(lapply(split(frags, frags$s), function(x){
  
  #message(x$posid[1])
  #if(x$posid[1] == 'chr2+46422396') browser()
  
  if(opt$buildStdFragments_categorize_anchorRead_remnants) x$posid <- paste0(x$posid, '.', x$repLeaderSeqGroup)
  
  if(opt$buildSites_sonicLengthAbund){
    estAbund <- sonicLengths[sonicLengths$s == x$s[1],]$estAbund
  }else{
    estAbund <- sum(unlist(lapply(split(x, x$uniqueSample), function(x) n_distinct(x$fragWidth))))
  }
  
  molCodes <- NA
  if(opt$demultiplex_captureRandomLinkerSeq){
    molCodes <- sum(unlist(lapply(split(x, x$uniqueSample), function(x){
      molCodes <- randomIDs[match(unlist(x$readIDlist), randomIDs$id),]$seq
      n_distinct(conformMinorSeqDiffs(molCodes, editDist = opt$buildSites_molCodes_minEditDist, abundSeqMinCount = opt$buildSites_molCodes_abundSeqMinCount))
    })))
  }
  
  if('flags' %in% names(samples)){
    x$flags <- paste0(unique(subset(samples, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1])$flags), collapse = ',')
  } else {
    x$flags <- NA
  }
  
  if(nrow(x) > 1){
    i <- rep(TRUE, nrow(x))
    
    # Check that leader sequences are all similar now that we are combining fragments from different replicates.
    r <- representativeSeq(x$repLeaderSeq)
    
    if(r[[1]] > opt$buildStdFragments_maxLeaderSeqDiffScore){
      # There is a conflict, one or more fragments have a markedly different leader sequences then the other fragments.
      # Attempt to salvage this site by retaining the majority of fragments with similar leader sequences.
      
      i <- as.vector(stringdist::stringdistmatrix(x$repLeaderSeq, r[[2]]) / nchar(x$repLeaderSeq) <= opt$buildStdFragments_maxLeaderSeqDiffScore)
      
      if(sum(i)/nrow(x) >= opt$buildSites_assemblyConflictResolution){
        x <- x[i,]
        r <- representativeSeq(x$repLeaderSeq)
      } else {
        return(tibble())
      }
    }
    
    # EstAbund and molCodes need to be recalculated since fragments were removed.
    if(opt$buildSites_sonicLengthAbund){
      estAbund <- sonicLengths[sonicLengths$s == x$s[1],]$estAbund
    }else{
      estAbund <- sum(unlist(lapply(split(x, x$uniqueSample), function(x) n_distinct(x$fragWidth))))
    }
    
    molCodes <- NA
    if(opt$demultiplex_captureRandomLinkerSeq){
      molCodes <- sum(unlist(lapply(split(x, x$uniqueSample), function(x){
        molCodes <- randomIDs[match(unlist(x$readIDlist), randomIDs$id),]$seq
        n_distinct(conformMinorSeqDiffs(molCodes, editDist = opt$buildSites_molCodes_minEditDist, abundSeqMinCount = opt$buildSites_molCodes_abundSeqMinCount))
      })))
    }
    
    return(dplyr::mutate(x, estAbund = estAbund, linkerMolCodes = molCodes, position = ifelse(strand[1] == '+', fragStart[1], fragEnd[1]), 
                         reads = sum(reads), repLeaderSeq = r[[2]], fragmentsRemoved = sum(!i)) %>%
             dplyr::select(trial, subject, sample, replicate, chromosome, strand, position, posid, estAbund, linkerMolCodes, reads, fragmentsRemoved, repLeaderSeq, flags) %>%
             dplyr::slice(1))
  }else{
    return(dplyr::mutate(x, estAbund = estAbund, linkerMolCodes = molCodes, position = ifelse(strand[1] == '+', fragStart[1], fragEnd[1]), fragmentsRemoved = 0) %>%
             dplyr::select(trial, subject, sample, replicate, chromosome, strand, position, posid, estAbund, linkerMolCodes, reads, fragmentsRemoved, repLeaderSeq, flags))
  }
}))

if(! opt$demultiplex_captureRandomLinkerSeq){
  sites <- dplyr::select(sites, -linkerMolCodes)
}

if (opt$buildSites_level == 'sample'){
  sites <- dplyr::select(sites, -replicate)
} 

if (opt$buildSites_level == 'subject'){
  sites <- dplyr::select(sites, -replicate, -sample)
} 

saveRDS(sites, file.path(opt$outputDir, opt$buildSites_outputDir, opt$buildSites_outputFile))

q(save = 'no', status = 0, runLast = FALSE) 
