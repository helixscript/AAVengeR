library(dplyr)
library(lubridate)
library(parallel)
library(data.table)
library(GenomicRanges)
library(Biostrings)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)

source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$buildStdFragments_outputDir))

write(c(paste(now(), 'Reading in fragment file(s).')), file = file.path(opt$outputDir, 'log'), append = TRUE)

samples <- loadSamples()

if('databaseGroup' %in% names(opt)){
  library(RMariaDB)
  library(DBI)
  con <- dbConnect(RMariaDB::MariaDB(), group = opt$databaseGroup)
  
  write(c(paste(now(), '  Reading in fragment data from database.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
  # Retrieve all fragments for unique trial / subject pairs.
  # This may retrieve fragments for samples outside of the current sample list but 
  # this is intended since we need to standardize sites across all samples.
  o <- distinct(dplyr::select(samples, trial, subject))
  frags <- bind_rows(lapply(split(o, 1:nrow(o)), function(x){
    r <- dbGetQuery(con, paste0("select * from fragments where subject='", x$subject, "' and trial='", x$trial, "'"))
    r2 <- tibble()
    if(nrow(r) > 0){
      r2 <- bind_rows(lapply(split(r, 1:nrow(r)), function(x2){
         d <- unserialize(x2$data[[1]])
         d$refGenome <- x2$refGenome
         d
      }))
    }
    r2
  }))
  
  dbDisconnect(con)
} else {
  opt$buildStdFragments_inputFiles <- gsub('\\s*,\\s*', ',', opt$buildStdFragments_inputFiles)

  frags <- bind_rows(lapply(unlist(strsplit(opt$buildStdFragments_inputFiles, ',')), function(x){
             write(c(paste(now(), '  reading in local data file ', x)), file = file.path(opt$outputDir, 'log'), append = TRUE)
             readRDS(file.path(opt$outputDir, x))
           }))
}

frags$uniqueSample <- paste0(frags$trial, '~', frags$subject, '~', frags$sample, '~', frags$replicate)
frags$fragID <- paste0(frags$trial, ':', frags$subject, ':', frags$sample, ':', frags$replicate, ':', frags$chromosome, ':', frags$strand, ':', frags$fragStart, ':', frags$fragEnd)


cluster <- makeCluster(opt$buildStdFragments_CPUs)
clusterExport(cluster, 'opt')

if(opt$buildStdFragments_categorize_anchorRead_remnants){
  write(c(paste(now(), 'Categorizing leader sequences.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  n <- 1
  frags$repLeaderSeqGroup <- 0

  while(any(frags$repLeaderSeqGroup == 0)){
      o <- dplyr::filter(frags, repLeaderSeqGroup == 0) %>% top_n(1, wt = reads) %>% dplyr::slice(1)
     
      i1 <- stringdist::stringdist(o$repLeaderSeq, frags$repLeaderSeq) <= ceiling(nchar(o$repLeaderSeq)*opt$buildStdFragments_categorize_anchorRead_remnants_percentDiff)
      i2 <- frags$repLeaderSeqGroup == 0
      frags[i1 & i2,]$repLeaderSeqGroup <- n
      n <- n+1
  }
} else {
  frags$repLeaderSeqGroup <- 0
}


# Standardize fragments.
# MN01490:66:000H3T5LT:1:22103:4969:7952
write(c(paste(now(), 'Standardizing fragment boundaries.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
frags_std <- standardizedFragments(frags, opt, cluster)

if(opt$buildStdFragments_standardizeSites & opt$buildStdFragments_standardizeBreaks){
  frags_std <- subset(frags_std, intSiteRefined == TRUE & breakPointRefined == TRUE) %>% dplyr::select(-intSiteRefined, -breakPointRefined, -width)
} else if(opt$buildStdFragments_standardizeSites == TRUE & opt$buildStdFragments_standardizeBreaks == FALSE){
  frags_std <- subset(frags_std, intSiteRefined == TRUE & breakPointRefined == FALSE) %>% dplyr::select(-intSiteRefined, -breakPointRefined, -width)
} else if(opt$buildStdFragments_standardizeSites == FALSE & opt$buildStdFragments_standardizeBreaks == TRUE){
  frags_std <- subset(frags_std, intSiteRefined == FALSE & breakPointRefined == TRUE) %>% dplyr::select(-intSiteRefined, -breakPointRefined, -width)
} else {
  # Something went wrong, eg. parameter missing or not set to a Boolean.
  q(save = 'no', status = 1, runLast = FALSE) 
}


# Join the standardized tallied fragments back onto the read level fragments in order to look for multi-hits.
write(c(paste(now(), 'Applying standardized boundaries.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
frags <- left_join(frags, dplyr::select(frags_std, fragID, start, end), by = 'fragID')

frags$fragStart <- frags$start
frags$fragEnd <- frags$end
frags$start <- NULL
frags$end <- NULL
frags$fragID <- paste0(frags$trial, ':', frags$subject, ':', frags$sample, ':', frags$replicate, ':', frags$chromosome, ':', frags$strand, ':', frags$fragStart, ':', frags$fragEnd)
frags$posid <-paste0(frags$chromosome, frags$strand, ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd))

multiHitFrags <- tibble()

# save.image('~/buildStdFrags.RData')

if(opt$buildStdFragments_salvageMultiHits){
  write(c(paste(now(), 'Salavaging multi-hit fragments.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
  # Expand fragments to read level.
  z <- tidyr::unnest(frags, readIDlist)

  # Identify reads which map to more than position id and define these as multi-hit reads.
  o <- group_by(z, readIDlist) %>%
       summarise(n = n_distinct(fragID)) %>%
       ungroup() %>%
       dplyr::filter(n > 1)

  if(nrow(o) > 0){
    multiHitFrags <- subset(z, readIDlist %in% o$readIDlist)  # Reads which map to more than one intSite positions.
    frags <- subset(z, ! readIDlist %in% o$readIDlist)        # Reads which map to a single intSite positions.
    
    browser()
  
    # There may be multi-hit reads where one of their position ids is the same as those in the frags
    # where the reads aligned uniquely. If only one of the possible alignment positions in multiHitFrags 
    # is in frags (unique alignments), then move those reads to frags.
    #
    # What if a read maps to more than fragment associatd with the same posid???
  
    if(any(multiHitFrags$posid %in% frags$posid)){
      multiHitFrags$returnToFrags <- FALSE
    
      # Create a table of valid posids for each subject.
      frags_sites <- dplyr::select(frags, trial, subject, posid) %>% dplyr::distinct()
      clusterExport(cluster, 'frags_sites')
    
      # Split multi-hit frags by read id and test if reads can be salvaged.
      # multiHitFrags <- bind_rows(parLapply(cluster, split(multiHitFrags, multiHitFrags$readIDlist), function(x){
      multiHitFrags <- bind_rows(lapply(split(multiHitFrags, multiHitFrags$readIDlist), function(x){
        
                                  # x will contain two or more frag ids supported by a SINGLE read.
                                  # Create a list of non-ambiguous sites from the subject this read came from.
                                  posidList <- subset(frags_sites, trial == x$trial[1] & subject == x$subject[1])$posid
                                
                                  # Return sites to the non-ambiguous read list if they map to a single site from the 
                                  # non-ambiguous site list.
                                  if(sum(unique(x$posid) %in% posidList) == 1){
                                    x[x$posid %in% posidList,]$returnToFrags <- TRUE
                                    
                                    if(length(unique(x[x$posid %in% posidList,]$fragID)) > 1) x[x$posid %in% posidList,]$returnToFrags <- FALSE 
                                  }
                                  x
                                }))
    
      if(any(multiHitFrags$returnToFrags == TRUE)){
        
        # Isolate frag reads that should be returned.
        m <- subset(multiHitFrags, returnToFrags == TRUE)
        write(paste0(sprintf("%.2f%%", (nrow(m) / nrow(multiHitFrags))*100), ' of multihit reads recovered because one potential site was in the unabiguous site list.'), file = file.path(opt$outputDir, 'log'), append = TRUE)

        # Remove salvaged reads from multi-hit reads.
        multiHitFrags <- subset(multiHitFrags, ! readIDlist %in% m$readIDlist)
        
        m$returnToFrags <- NULL
        multiHitFrags$returnToFrags <- NULL
        
        # Return salvaged reads to frags.
        frags <- bind_rows(frags, m)
      }
    }
    
    # Write out multiHitFrags here...
    
    # Regroup reads.
    frags <- group_by(frags, uniqueSample, chromosome, strand, fragStart, fragEnd) %>% 
             mutate(reads = n_distinct(readIDlist), readIDlist = list(readIDlist)) %>% 
             dplyr::slice(1) %>% ungroup()
  }
} else {
  
  # Expand fragments to the read level.
  z <- tidyr::unnest(frags, readIDlist)
  
  # Find reads that map to a single standardized fragment.
  o <- group_by(z, readIDlist) %>%
       summarise(n = n_distinct(fragID)) %>%
       ungroup() %>%
       dplyr::filter(n == 1)
  
  if(nrow(o) == 0){
    stop('Error -- no reads remain after removing reads that map to multiple std fragments.')
  }
  
  frags <- group_by(subset(z, readIDlist %in% o$readIDlist), uniqueSample, chromosome, strand, fragStart, fragEnd) %>% 
           mutate(reads = n_distinct(readIDlist), readIDlist = list(readIDlist)) %>% 
           dplyr::slice(1) %>% ungroup()
}
  

# Clear out the tmp/ directory.
invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))


# Redefine the cluster since it goes away sometimes.
stopCluster(cluster)
gc()

cluster <- makeCluster(opt$buildStdFragments_CPUs)
clusterExport(cluster, 'opt')

f <- as.data.table(frags)
s <- split(f, f$fragID)
 
write(c(paste(now(), 'Regrouping standardized fragments.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

frags <- bind_rows(parLapply(cluster, s, function(a){ 
#frags <- bind_rows(lapply(s, function(a){   
  library(dplyr)
  source(file.path(opt$softwareDir, 'lib.R'))
  
  r <- representativeSeq(a$repLeaderSeq)
  
  # Exclude reads where the leaderSeq is not similar to the representative sequence. 
  i <- stringdist::stringdist(r[[2]], a$repLeaderSeq) / nchar(r[[2]]) <= opt$buildStdFragments_maxLeaderSeqDiffScore
  if(all(! i)) return(tibble::tibble())
  
  a <- a[i,]
  
  a$repLeaderSeq <- r[[2]]
  a$reads <- sum(a$reads)
  
  reads <- unique(unlist(a$readIDlist))
  a <- a[1,]
  a$readIDlist <- list(reads)
  a
}))


if(nrow(multiHitFrags) > 0){
  write(c(paste(now(), 'Preparing mult-hit data files.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
  multiHitReads <- group_by(multiHitFrags, uniqueSample, readIDlist) %>%
    summarise(posids = list(posid)) %>%
    ungroup() 
  
  multiHitFrags <- bind_rows(parLapply(cluster, split(multiHitFrags, multiHitFrags$fragID), function(a){
  #multiHitFrags <- bind_rows(lapply(split(multiHitFrags, multiHitFrags$fragID), function(a){
    library(dplyr)
    library(data.table)
    source(file.path(opt$softwareDir, 'lib.R'))
    
    r <- representativeSeq(a$repLeaderSeq)
    
    # Exclude reads where the leaderSeq is not similar to the representative sequence. 
    i <- stringdist::stringdist(r[[2]], a$repLeaderSeq) / nchar(r[[2]]) <= opt$buildStdFragments_maxLeaderSeqDiffScore
    if(all(! i)) return(tibble())
    
    a <- a[i,]
    
    a$repLeaderSeq <- r[[2]]
    
    reads <- n_distinct(a$readIDlist)
    readIDlist <- unique(a$readIDlist)
    
    a <- a[1,]
    a$readIDlist <- list(readIDlist)
    a
  }))
  
  saveRDS(multiHitReads, file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitReads.rds'))
  saveRDS(multiHitFrags, file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitFrags.rds'))
}


# Clear out the tmp/ directory.
invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

saveRDS(data.frame(frags), file.path(opt$outputDir, opt$buildStdFragments_outputDir, opt$buildStdFragments_outputFile))

q(save = 'no', status = 0, runLast = FALSE) 