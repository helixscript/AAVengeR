library(dplyr)
library(parallel)
library(data.table)
library(GenomicRanges)

opt <- yaml::read_yaml('config.yml')
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$buildStdFragments_outputDir))

dups <- tibble()

if('buildStdFragments_duplicateReadFile' %in% names(opt)){
  dups <- readRDS(file.path(opt$outputDir, opt$buildStdFragments_duplicateReadFile))
  dups <- data.table(dplyr::distinct(dplyr::select(dups, id, n)))
}

# Read in read level fragments records.
frags <- readRDS(file.path(opt$outputDir, opt$buildFragments_outputDir, opt$buildFragments_outputFile))


cluster <- makeCluster(opt$buildFragments_CPUs)
clusterExport(cluster, c('opt', 'dups'))

frags$fragID <- paste0(frags$trial, ':', frags$subject, ':', frags$sample, ':', frags$replicate, ':', frags$strand, ':', frags$fragStart, ':', frags$fragEnd)

# Tally the number of reads for each fragment which is needed for standardization.
f1 <- data.table(dplyr::select(frags, fragID, readID))
f2 <- bind_rows(parLapply(cluster, split(f1, f1$fragID), function(a){
  library(data.table)
  if(nrow(a) == 1){
    a$reads = 1
  } else {
    a$reads <- nrow(a) 
    o <- dups[dups$id %in% a$readID]
    if(nrow(o) > 0) a$reads <- a$reads + sum(o$n)
  }
  
  a[1, c(1, 3)]
}))


frags2 <- left_join(frags, f2, by = 'fragID') %>% dplyr::select(-readID) %>% dplyr::distinct()



cluster <- makeCluster(opt$buildStdFragments_CPUs)
clusterExport(cluster, 'opt')

# Standardize fragments.
frags_std <- standardizedFragments(frags2, opt, cluster)
frags_std <- subset(frags_std, intSiteRefined == TRUE & breakPointRefined == TRUE) %>% dplyr::select(-intSiteRefined, -breakPointRefined, -width)

# Join the standardized tallied fragments back onto the read level fragments in order to look for multi-hits.
frags <- left_join(frags, dplyr::select(frags_std, fragID, start, end), by = 'fragID')

frags$fragStart <- frags$start
frags$fragEnd <- frags$end
frags$start <- NULL
frags$end <- NULL
frags$fragID <- paste0(frags$trial, ':', frags$subject, ':', frags$sample, ':', frags$replicate, ':', frags$strand, ':', frags$fragStart, ':', frags$fragEnd)
frags$posid <-paste0(frags$chromosome, frags$strand, ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd))


# Identify reads which map to more than position id and define these as multihit reads.
o <- group_by(frags, readID) %>%
  summarise(n = n_distinct(posid)) %>%
  ungroup() %>%
  dplyr::filter(n > 1)

multiHitFrags <- tibble()

if(nrow(o) > 0){
  multiHitFrags <- subset(frags, readID %in% o$readID)  # Reads which map to more than one intSite positions.
  frags <- subset(frags, ! readID %in% o$readID)        # Reads which map to a single intSite positions.
  
  # There may be multi hit reads where one of their position ids is the same as those in the frags
  # where the reads aligned uniquely. If only one of the possible alignment positions in multiHitFrags 
  # is in frags (unique alignments), then move those reads to frags.
  
  
  if(any(multiHitFrags$posid %in% frags$posid)){
    multiHitFrags$returnToFrags <- FALSE
    
    frags_sites <- dplyr::select(frags, trial, subject, posid) %>% dplyr::distinct()
    clusterExport(cluster, 'frags_sites')
    
    multiHitFrags <- bind_rows(parLapply(cluster, split(multiHitFrags, multiHitFrags$readID), function(x){
                                o <- subset(frags_sites, trial == x$trial[1] & subject == x$subject[1]) # List of posids from the trials and subject these reads came from.
                                
                                if(sum(unique(x$posid) %in% o$posid) == 1){
                                  x[x$posid %in% o$posid,]$returnToFrags <- TRUE
                                }
                                x
                              }))
    
    if(any(multiHitFrags$returnToFrags == TRUE)){
      m <- subset(multiHitFrags, returnToFrags == TRUE)
      write(paste0(sprintf("%.2f%%", (nrow(m) / nrow(multiHitFrags))*100), ' of multihit reads recovered because one potential site was in the unabiguous site list.'), file = file.path(opt$outputDir, 'log'), append = TRUE)
      multiHitFrags <- subset(multiHitFrags, ! readID %in% m$readID)
      m$returnToFrags <- NULL
      multiHitFrags$returnToFrags <- NULL
      frags <- bind_rows(frags, m)
    }
  }
}


# Clear out the tmp/ directory.
invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))


# Redefine the cluster since it goes away sometimes.
stopCluster(cluster)
gc()


cluster <- makeCluster(opt$buildStdFragments_CPUs)
clusterExport(cluster, c('opt', 'dups'))

# s <- split(frags, frags$fragID)
library(data.table)
f <- as.data.table(frags)
s <- split(f, f$fragID)
 
#frags <- bind_rows(lapply(s, function(a){   
frags <- bind_rows(parLapply(cluster, s, function(a){     
  library(dplyr)
  library(data.table)  # Needed for dup table.
  source(file.path(opt$softwareDir, 'lib.R'))
  
  r <- representativeSeq(a$leaderSeq.anchorReads)
  
  # Exclude reads where the leaderSeq is not similiar to the representative sequence. 
  i <- stringdist::stringdist(r[[2]], a$leaderSeq.anchorReads) / nchar(r[[2]]) <= opt$buildStdFragments_maxLeaderSeqDiffScore
  if(all(! i)) return(tibble::tibble())
  
  a <- a[i,]
  
  a$repLeaderSeq <- r[[2]]
  
  a$reads <- dplyr::n_distinct(a$readID)
  o <- dups[dups$id %in% a$readID]
  
  if(nrow(o) > 0) a$reads <- a$reads + sum(o$n)

  a[1, c(7:10, 2:5, 13:14)]
}))



if(nrow(multiHitFrags) > 0){
  
  multiHitReads <- group_by(multiHitFrags, trial, subject, sample, readID) %>%
    summarise(posids = I(list(posid))) %>%
    ungroup() %>%
    tidyr::unnest(posids) %>%
    dplyr::distinct()
  
  multiHitFrags <- bind_rows(parLapply(cluster, split(multiHitFrags, multiHitFrags$fragID), function(a){
  #multiHitFrags <- bind_rows(lapply(split(multiHitFrags, multiHitFrags$fragID), function(a){
    library(dplyr)
    library(data.table)
    source(file.path(opt$softwareDir, 'lib.R'))
    #message(a$fragID, '\n')
    
    r <- representativeSeq(a$leaderSeq.anchorReads)
    
    # Exclude reads where the leaderSeq is not similiar to the representative sequence. 
    i <- stringdist::stringdist(r[[2]], a$leaderSeq.anchorReads) / nchar(r[[2]]) <= opt$buildStdFragments_maxLeaderSeqDiffScore
    if(all(! i)) return(tibble())
    
    a <- a[i,]
    a$repLeaderSeq <- r[[2]]
    
    a$reads <- dplyr::n_distinct(a$readID)
    
    o <- dups[dups$id %in% a$readID]
    if(nrow(o) > 0) a$reads <- a$reads + sum(o$n)
    
    a[1, c(7:10, 2:5, 13:14)]
  }))
  
  saveRDS(multiHitReads, file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitReads.rds'))
  saveRDS(multiHitFrags, file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitFrags.rds'))
}


# Clear out the tmp/ directory.
invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

saveRDS(frags, file.path(opt$outputDir, opt$buildStdFragments_outputDir, opt$buildStdFragments_outputFile))

q(save = 'no', status = 0, runLast = FALSE) 