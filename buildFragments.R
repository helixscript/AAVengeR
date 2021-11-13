library(dplyr)
library(parallel)
library(data.table)
library(GenomicRanges)

opt <- yaml::read_yaml('config.yml')
source(file.path(opt$softwareDir, 'lib.R'))

dups <- tibble()
if('buildFragments_duplicateReadFile' %in% names(opt)){
  dups <- readRDS(file.path(opt$outputDir, opt$buildFragments_duplicateReadFile))
  dups <- data.table(dplyr::distinct(dplyr::select(dups, id, n)))
}

anchorReadAlignments <- readRDS(file.path(opt$outputDir, opt$buildFragments_anchorReadsAlignmentFile))
adriftReadAlignments <- readRDS(file.path(opt$outputDir, opt$buildFragments_adriftReadsAlignmentFile))

dir.create(file.path(opt$outputDir, opt$buildFragments_outputDir))

anchorReadAlignments <- select(anchorReadAlignments, sample, qName, tName, strand, tStart, tEnd, leaderSeq)
adriftReadAlignments <- select(adriftReadAlignments, sample, qName, tName, strand, tStart, tEnd)

# Some reads will have so many alignments after filtering that it will break the joining of the anchor and adrift reads
# tables and should be abandoned. Setting buildFragments_maxAlignmentsPerRead to a high value will include 
# these values but buildFragments_idGroup_size should be lowered to protect against breaking the joining 
# function by exceeding R's built in max table size limit.

a <- table(anchorReadAlignments$qName)
b <- table(adriftReadAlignments$qName)
x <- unique(c(names(a[a > opt$buildFragments_maxAlignmentsPerRead]), names(b[b > opt$buildFragments_maxAlignmentsPerRead])))

write(x, file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'highAlignmentReads'))

anchorReadAlignments <- subset(anchorReadAlignments, ! qName %in% x)
adriftReadAlignments <- subset(adriftReadAlignments, ! qName %in% x)

names(anchorReadAlignments) <- paste0(names(anchorReadAlignments), '.anchorReads')
names(adriftReadAlignments) <- paste0(names(adriftReadAlignments), '.adriftReads')

ids <- unique(anchorReadAlignments$qName.anchorReads)
id_groups <- split(ids, dplyr::ntile(1:length(ids), round(length(ids)/opt$buildFragments_idGroup_size)))

cluster <- makeCluster(opt$buildFragments_CPUs)
clusterExport(cluster, c('opt', 'anchorReadAlignments', 'adriftReadAlignments'))

# Build intial fragments.
frags <- bind_rows(parLapply(cluster, id_groups, function(id_group){
  library(dplyr)
  a <- subset(anchorReadAlignments, qName.anchorReads %in% id_group)
  b <- subset(adriftReadAlignments, qName.adriftReads %in% id_group)
  
  if(nrow(a) == 0 | nrow(b) == 0) return(data.frame())
  
  # Join adrift reads alignemnts to anchor read alignments to create potential read pairs.
  frags <- left_join(a, b, by = c('qName.anchorReads' = 'qName.adriftReads')) %>% tidyr::drop_na()
  
  # Remove combinations not found on the same chromosome. 
  i <- which(frags$tName.anchorReads != frags$tName.adriftReads)
  if(length(i) > 0) frags <- frags[-i,]
  if(nrow(frags) == 0) return(data.frame())
  
  # Remove combinations which have the same strand since fragment reads are expected to have opposite strands.
  i <- which(frags$strand.anchorReads == frags$strand.adriftReads)
  if(length(i) > 0) frags <- frags[-i,]
  if(nrow(frags) == 0) return(data.frame())
  
  # Determine the start and end of fragments based on their alignment strands
  # and perform some sanity tests then filter on fragment size. 
  mutate(frags, 
         fragStart  = ifelse(strand.anchorReads == '+', tStart.anchorReads + 1, tStart.adriftReads + 1),
         fragEnd    = ifelse(strand.anchorReads == '+', tEnd.adriftReads + 1,   tEnd.anchorReads + 1),
         strand     = ifelse(strand.anchorReads == '+', '+', '-'),
         chromosome = tName.anchorReads,  
         fragTest  = ifelse(strand.anchorReads == '+', tStart.anchorReads < tEnd.adriftReads, tStart.adriftReads < tEnd.anchorReads),  
         fragWidth = (fragEnd - fragStart) + 1) %>%
         filter(fragTest == TRUE, 
                fragWidth <= opt$buildFragments_maxFragLength,
                fragWidth >= opt$buildFragments_minFragLength) %>%
                mutate(uniqueSample = sample.anchorReads, readID = qName.anchorReads, fragID = paste0(uniqueSample, ':', chromosome, ':', strand, ':', fragStart, ':', fragEnd)) %>%
                select(fragID, uniqueSample, readID, chromosome, strand, fragStart, fragEnd, leaderSeq.anchorReads)
}))


# Stop and rebuild cluster objects to release alignment objects from memory.
stopCluster(cluster)
gc()
cluster <- makeCluster(opt$buildFragments_CPUs)
clusterExport(cluster, c('opt', 'dups'))


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

# Unpack ids and build a data frame for standardizaing fragments.
o <- stringr::str_split(f2$fragID, ':')
i <- stringr::str_split(f2$fragID, '[~:]')
d <- tibble(uniqueSample = unlist(lapply(o, '[[', 1)),
            subject      = unlist(lapply(i, '[[', 1)),
            sample       = unlist(lapply(i, '[[', 2)),
            replicate    = unlist(lapply(i, '[[', 3)),
            seqnames     = unlist(lapply(o, '[[', 2)),
            strand       = unlist(lapply(o, '[[', 3)),
            start        = unlist(lapply(o, '[[', 4)),
            end          = unlist(lapply(o, '[[', 5)),
            fragID       = f2$fragID,
            reads        = f2$reads)


# Standardize fragments.
frags_std <- standardizedFragments(d, opt, cluster)
frags_std <- subset(frags_std, intSiteRefined == TRUE & breakPointRefined == TRUE) %>% dplyr::select(-intSiteRefined, -breakPointRefined, -width)


# Expand standardized fragments back to the read level.
frags_std$fragID2 <- paste0(frags_std$uniqueSample, ':', frags_std$seqnames, ':', frags_std$strand, ':', frags_std$start, ':', frags_std$end)
frags <- left_join(frags, dplyr::select(frags_std, fragID, fragID2, start, end), by = 'fragID')


# Over-write previous fragment boundaries with standardized boundary paositions.
frags$fragID <- frags$fragID2
frags$fragStart <- frags$start
frags$fragEnd <- frags$end
frags$fragID2 <- NULL
frags$start <- NULL
frags$end <- NULL

# Create fragment position ids.
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
    clusterExport(cluster, 'frags')
    multiHitFrags <- bind_rows(parLapply(cluster, split(multiHitFrags, multiHitFrags$readID), function(x){

      if(sum(x$posid %in% frags$posid) == 1){
        x[x$posid %in% frags$posid,]$returnToFrags <- TRUE
      }
      
      x
    }))
    
    if(any(multiHitFrags$returnToFrags == TRUE)){
      m <- subset(multiHitFrags, returnToFrags == TRUE)
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
cluster <- makeCluster(opt$buildFragments_CPUs)
clusterExport(cluster, c('opt'))

# Collapse read level fragment records and determine the concensus leader sequences.
frags <- bind_rows(parLapply(cluster, split(frags, frags$fragID), function(a){
       library(dplyr)
       source(file.path(opt$softwareDir, 'lib.R'))

       r <- representativeSeq(a$leaderSeq.anchorReads)
    
       # Exclude reads where the leaderSeq is not similiar to the representative sequence. 
       i <- stringdist::stringdist(r[[2]], a$leaderSeq.anchorReads) / nchar(r[[2]]) <= opt$buildFragments_maxLeaderSeqDiffScore
       if(all(! i)) return(tibble())
    
       a <- a[i,]
       a$repLeaderSeq <- r[[2]]
    
       a$reads <- n_distinct(a$readID)
       
       bind_cols(unpackUniqueSampleID(a[1,2]),  a[1, c(4:7, 11, 9, 10)])
   }))


if(nrow(multiHitFrags) > 0){
  
  multiHitReads <- group_by(unpackUniqueSampleID(multiHitFrags), subject, sample, readID) %>%
                   summarise(posids = I(list(posid))) %>%
                   ungroup() %>%
                   tidyr::unnest(posids) %>%
                   dplyr::distinct()
  
   #multiHitFrags <- bind_rows(parLapply(cluster, split(multiHitFrags, multiHitFrags$fragID), function(a){
   multiHitFrags <- bind_rows(lapply(split(multiHitFrags, multiHitFrags$fragID), function(a){
     library(dplyr)
     source(file.path(opt$softwareDir, 'lib.R'))
     #message(a$fragID, '\n')
  
     r <- representativeSeq(a$leaderSeq.anchorReads)
  
     # Exclude reads where the leaderSeq is not similiar to the representative sequence. 
     i <- stringdist::stringdist(r[[2]], a$leaderSeq.anchorReads) / nchar(r[[2]]) <= opt$buildFragments_maxLeaderSeqDiffScore
     if(all(! i)) return(tibble())
  
     a <- a[i,]
     a$repLeaderSeq <- r[[2]]
  
     a$reads <- n_distinct(a$readID)
  
    bind_cols(unpackUniqueSampleID(a[1,2]),  a[1, c(4:7, 11, 9, 10)])
   }))
   
   saveRDS(multiHitReads, file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'multiHitReads.rds'))
   saveRDS(multiHitFrags, file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'multiHitFrags.rds'))
}


# Clear out the tmp/ directory.
invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

saveRDS(frags, file.path(opt$outputDir, opt$buildFragments_outputDir, opt$buildFragments_outputFile))

