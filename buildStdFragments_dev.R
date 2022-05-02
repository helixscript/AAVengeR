library(dplyr)
library(lubridate)
library(parallel)
library(data.table)
library(GenomicRanges)
library(Biostrings)

configFile <- '/data/project/Encoded/220222_M03249_0243_000000000-JHTY6/config.yml'
opt <- yaml::read_yaml(configFile)

source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$buildStdFragments_outputDir))

# Read in read level fragments records.

opt$buildStdFragments_inputFiles <- 'buildFragments.org/fragments.rds'

frags <- bind_rows(lapply(unlist(strsplit(opt$buildStdFragments_inputFiles, ',')), function(x){
           message('reading in ', x)
           readRDS(file.path(opt$outputDir, x))
         }))

r <- readLines('chr17+27974433_readIDS')

frags$posid <- paste0(frags$chromosome, frags$strand, ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd))

z <- subset(frags, posid == 'chr17+27974433') %>% tidyr::unnest(readIDlist)

f <- tidyr::unnest(frags, readIDlist)
o <- subset(f,  readIDlist %in% z$readIDlist)
sort(table(o$posid), decreasing = TRUE)[1:50]

# The reads supporting chr17+27974433 also support hundred of other posids.





frags <- tidyr::unnest(frags, readIDlist)
frags <- subset(frags, readIDlist %in% r)

frags <- tidyr::nest(frags, readIDlist = c(readIDlist))


cluster <- makeCluster(opt$buildStdFragments_CPUs)
clusterExport(cluster, 'opt')

if(opt$buildStdFragments_categorize_anchorRead_remnants){
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
frags_std <- standardizedFragments(frags, opt, cluster)

if(opt$buildStdFragments_standardizeSites & opt$buildStdFragments_standardizeBreaks){
  frags_std <- subset(frags_std, intSiteRefined == TRUE & breakPointRefined == TRUE) %>% dplyr::select(-intSiteRefined, -breakPointRefined, -width)
} else if(opt$buildStdFragments_standardizeSites == TRUE & opt$buildStdFragments_standardizeBreaks == FALSE){
  frags_std <- subset(frags_std, intSiteRefined == TRUE & breakPointRefined == FALSE) %>% dplyr::select(-intSiteRefined, -breakPointRefined, -width)
} else if(opt$buildStdFragments_standardizeSites == FALSE & opt$buildStdFragments_standardizeBreaks == TRUE){
  frags_std <- subset(frags_std, intSiteRefined == FALSE & breakPointRefined == TRUE) %>% dplyr::select(-intSiteRefined, -breakPointRefined, -width)
} else {
  # Something went wrong, eg. parameter missing or not set to a boolean.
  q(save = 'no', status = 1, runLast = FALSE) 
}


# Join the standardized tallied fragments back onto the read level fragments in order to look for multi-hits.
frags <- left_join(frags, dplyr::select(frags_std, fragID, start, end), by = 'fragID')

frags$fragStart <- frags$start
frags$fragEnd <- frags$end
frags$start <- NULL
frags$end <- NULL
frags$fragID <- paste0(frags$trial, ':', frags$subject, ':', frags$sample, ':', frags$replicate, ':', frags$strand, ':', frags$fragStart, ':', frags$fragEnd)
frags$posid <-paste0(frags$chromosome, frags$strand, ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd))

# Expand fragments to read level.
z <- tidyr::unnest(frags, readIDlist)

# Identify reads which map to more than position id and define these as multihit reads.
o <- group_by(z, readIDlist) %>%
  summarise(n = n_distinct(posid)) %>%
  ungroup() %>%
  dplyr::filter(n > 1)

multiHitFrags <- tibble()

if(nrow(o) > 0){
  multiHitFrags <- subset(z, readIDlist %in% o$readIDlist)  # Reads which map to more than one intSite positions.
  frags <- subset(z, ! readIDlist %in% o$readIDlist)        # Reads which map to a single intSite positions.
  
  # There may be multihit reads where one of their position ids is the same as those in the frags
  # where the reads aligned uniquely. If only one of the possible alignment positions in multiHitFrags 
  # is in frags (unique alignments), then move those reads to frags.
  
  
  if(any(multiHitFrags$posid %in% frags$posid)){
    multiHitFrags$returnToFrags <- FALSE
    
    frags_sites <- dplyr::select(frags, trial, subject, posid) %>% dplyr::distinct()
    clusterExport(cluster, 'frags_sites')
    
    # Split multi hit frags by read id.
    multiHitFrags <- bind_rows(parLapply(cluster, split(multiHitFrags, multiHitFrags$readIDlist), function(x){
    #multiHitFrags <- bind_rows(lapply(split(multiHitFrags, multiHitFrags$readIDlist), function(x){
                                # x will contain two or more frag ids.
                                # Create a list of non-ambiguous sites from the subject this read came from.
                                posidList <- subset(frags_sites, trial == x$trial[1] & subject == x$subject[1])$posid
                                
                                
                                # Return sites to the non-ambigious read list if they map to a single site from the 
                                # non-ambiguous site list.
                                if(sum(unique(x$posid) %in% posidList) == 1){
                                  x[x$posid %in% posidList,]$returnToFrags <- TRUE
                                }
                                x
                              }))
    
    if(any(multiHitFrags$returnToFrags == TRUE)){
      m <- subset(multiHitFrags, returnToFrags == TRUE)
      write(paste0(sprintf("%.2f%%", (nrow(m) / nrow(multiHitFrags))*100), ' of multihit reads recovered because one potential site was in the unabiguous site list.'), file = file.path(opt$outputDir, 'log'), append = TRUE)
      multiHitFrags <- subset(multiHitFrags, ! readIDlist %in% m$readIDlist)
      m$returnToFrags <- NULL
      multiHitFrags$returnToFrags <- NULL
      frags <- bind_rows(frags, m)
    }
  }
}

frags <- tidyr::nest(frags, data = c(readIDlist))
frags$data <- NULL

# if(nrow(multiHitFrags) > 0){
#   multiHitFrags <- tidyr::nest(multiHitFrags, data = c(readIDlist))
#   multiHitFrags$data <- NULL
# }


# Clear out the tmp/ directory.
invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))


# Redefine the cluster since it goes away sometimes.
stopCluster(cluster)
gc()


cluster <- makeCluster(opt$buildStdFragments_CPUs)
clusterExport(cluster, 'opt')

f <- as.data.table(frags)
s <- split(f, f$fragID)
 
frags <- bind_rows(parLapply(cluster, s, function(a){ 
#frags <- bind_rows(lapply(s, function(a){   
  library(dplyr)
  source(file.path(opt$softwareDir, 'lib.R'))
  
  r <- representativeSeq(a$repLeaderSeq)
  
  # Exclude reads where the leaderSeq is not similiar to the representative sequence. 
  i <- stringdist::stringdist(r[[2]], a$repLeaderSeq) / nchar(r[[2]]) <= opt$buildStdFragments_maxLeaderSeqDiffScore
  if(all(! i)) return(tibble::tibble())
  
  a <- a[i,]
  
  a$repLeaderSeq <- r[[2]]
  a$reads <- sum(a$reads)
  
  a[1, c(6:10, 2:5, 13:14, 12, 15)]
}))



if(nrow(multiHitFrags) > 0){
  multiHitReads <- group_by(multiHitFrags, trial, subject, sample, readIDlist) %>%
    summarise(posids = I(list(posid))) %>%
    ungroup() %>%
    tidyr::unnest(posids) %>%
    dplyr::distinct()
  
  multiHitFrags <- bind_rows(parLapply(cluster, split(multiHitFrags, multiHitFrags$fragID), function(a){
  #multiHitFrags <- bind_rows(lapply(split(multiHitFrags, multiHitFrags$fragID), function(a){
    library(dplyr)
    library(data.table)
    source(file.path(opt$softwareDir, 'lib.R'))
    
    r <- representativeSeq(a$repLeaderSeq)
    
    # Exclude reads where the leaderSeq is not similiar to the representative sequence. 
    i <- stringdist::stringdist(r[[2]], a$repLeaderSeq) / nchar(r[[2]]) <= opt$buildStdFragments_maxLeaderSeqDiffScore
    if(all(! i)) return(tibble())
    
    a <- a[i,]
    
    a$repLeaderSeq <- r[[2]]
    a$reads <- sum(a$reads)
    
    a[1, c(6:10, 2:5, 13:14, 12, 15)]
  }))
  
  saveRDS(multiHitReads, file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitReads.rds'))
  saveRDS(multiHitFrags, file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitFrags.rds'))
}


# Clear out the tmp/ directory.
invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

saveRDS(data.frame(frags), file.path(opt$outputDir, opt$buildStdFragments_outputDir, opt$buildStdFragments_outputFile))

q(save = 'no', status = 0, runLast = FALSE) 