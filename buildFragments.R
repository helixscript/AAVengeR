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
                mutate(uniqueSample = sample.anchorReads, readID = qName.anchorReads) %>%
                select(uniqueSample, readID, chromosome, strand, fragStart, fragEnd, leaderSeq.anchorReads)
}))


# Random ids.
if(opt$demultiplex_captureRandomLinkerSeq){
  startingReadCount <- n_distinct(frags$readID)
  
  # Collect the random IDs isolated during demultiplexing.
  files <- list.files(file.path(opt$outputDir, opt$demultiplex_outputDir), pattern = 'randomIDs', full.names = TRUE)
  randomIDs <- Reduce('append', lapply(files,  readDNAStringSet))
  r <- tibble(readID = names(randomIDs), randomLinkerSeq = as.character(randomIDs))
  
  frags <- left_join(frags, r, by = 'readID')
  frags <- unpackUniqueSampleID(frags)
  
  # Random linker ids should be sample specific.
  frags$s <- paste(frags$trial, frags$subject, frags$sample)
  o <- group_by(frags, randomLinkerSeq) %>%
       summarise(nSamples = n_distinct(s), reads = n_distinct(readID)) %>%
       ungroup() %>%
       filter(nSamples > 1)
  
  # Remove reads where the random id was seen in two or more samples and there are 
  # too few reads were associated with the random id to attemept to segreagate. 
  z <- subset(o, reads < opt$buildFragments_randomLinkerID_minReadCountToSegreagate)$randomLinkerSeq
  if(length(z) > 0){
    frags <- subset(frags, ! randomLinkerSeq %in% z)
  }
  
  z <- subset(o, reads >= opt$buildFragments_randomLinkerID_minReadCountToSegreagate)$randomLinkerSeq
  
  if(length(z) > 0){
    a <- subset(frags, ! randomLinkerSeq %in% z)  # Undisputed reads
    b <- subset(frags, randomLinkerSeq %in% z)    # Reads which need to be examined.
    
    c <- bind_rows(lapply(split(b, b$randomLinkerSeq), function(x){
           tab <- data.frame(table(x$s))
           tab$percentReads <- (tab$Freq / n_distinct(x$readID))*100 
           topSample <- dplyr::arrange(tab, desc(percentReads)) %>% dplyr::slice(1)
           
           if(topSample$percentReads >= opt$buildFragments_randomLinkerID_minSingleSampleMajorityPercent){
             return(subset(x, s %in% topSample$Var1))
           } else{
             return(tibble())
           }
         }))
    
    if(nrow(c) > 0){
      frags <- bind_rows(a, c)
    } else {
      frags <- a
    }
  }
  
  frags <- dplyr::select(frags, -randomLinkerSeq, -trial, -subject, -sample, -replicate, -s)
  msg <- paste0(startingReadCount - n_distinct(frags$readID), ' reads, ', sprintf("%.2f%%", (1 - n_distinct(frags$readID) / startingReadCount)*100), ' of total reads, removed due to random linker ID conflicts.')
  write(msg, file = file.path(opt$outputDir, 'log'), append = TRUE)
}


# Stop and rebuild cluster objects to release alignment objects from memory.
stopCluster(cluster)
gc()

frags <- unpackUniqueSampleID(frags)

samples <- loadSamples()
frags <- left_join(frags, select(samples, uniqueSample, flags), by = 'uniqueSample')

dups <- tibble()

if('buildStdFragments_duplicateReadFile' %in% names(opt)){
  dups <- readRDS(file.path(opt$outputDir, opt$buildStdFragments_duplicateReadFile))
  dups <- data.table(dplyr::distinct(dplyr::select(dups, id, n)))
}

cluster <- makeCluster(opt$buildFragments_CPUs)
clusterExport(cluster, c('opt', 'dups'))

frags$fragID <- paste0(frags$trial, ':', frags$subject, ':', frags$sample, ':', frags$replicate, ':', frags$strand, ':', frags$fragStart, ':', frags$fragEnd)

# Combine read level fragments into unique fragments with read counts.
f1 <- data.table(frags)
f2 <- split(f1, f1$fragID)

frags <- parLapply(cluster, f2, function(x){
#frags <- (lapply(f2, function(x){  
           library(dplyr)
           library(data.table)
           source(file.path(opt$softwareDir, 'lib.R'))
  
           if(nrow(x) == 1){
             x$reads <- nrow(x)
             x$repLeaderSeq <- x$leaderSeq.anchorReads
             x$readIDlist <- x$readID
             return(select(as_tibble(x), -readID, -leaderSeq.anchorReads))
           }
  
           r <- representativeSeq(x$leaderSeq.anchorReads)
  
           # Exclude reads where the leaderSeq is not similiar to the representative sequence. 
           i <- stringdist::stringdist(r[[2]], x$leaderSeq.anchorReads) / nchar(r[[2]]) <= opt$buildStdFragments_maxLeaderSeqDiffScore
           if(all(! i)) return(data.table::data.table())
  
           x <- x[i,]
          
           x$reads <- nrow(x)
           x$repLeaderSeq <- r[[2]]
  
           o <- dups[dups$id %in% x$readID]
           if(nrow(o) > 0){
             x$reads <- nrow(x) + sum(o$n)
           }
           readIDs <- x$readID
           x <- as_tibble(x[1,])
           x$readIDlist <- list(readIDs)
           
           return(select(x[1,], -readID, -leaderSeq.anchorReads))
})


frags <- bind_rows(lapply(frags, function(x){
           if(! is.list(x$readIDlist)) x$readIDlist <- list(x$readIDlist)
           x
        }))

stopCluster(cluster)

saveRDS(distinct(frags), file.path(opt$outputDir, opt$buildFragments_outputDir, opt$buildFragments_outputFile))

q(save = 'no', status = 0, runLast = FALSE) 