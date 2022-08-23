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

write(c(paste(now(), 'Reading in alignment results.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
anchorReadAlignments <- readRDS(file.path(opt$outputDir, opt$buildFragments_anchorReadsAlignmentFile))
adriftReadAlignments <- readRDS(file.path(opt$outputDir, opt$buildFragments_adriftReadsAlignmentFile))

dir.create(file.path(opt$outputDir, opt$buildFragments_outputDir))

anchorReadAlignments <- select(anchorReadAlignments, sample, qName, tName, strand, tStart, tEnd, leaderSeq)
adriftReadAlignments <- select(adriftReadAlignments, sample, qName, tName, strand, tStart, tEnd)

names(anchorReadAlignments) <- paste0(names(anchorReadAlignments), '.anchorReads')
names(adriftReadAlignments) <- paste0(names(adriftReadAlignments), '.adriftReads')

# Filter read pairs such that one must mate must not be too 'jumpy'.
# Setting buildFragments_min_num_alignments_for_multiHit to one would require 
# one mate to map uniquely while the other mate could have multiple alignments. 

write(c(paste(now(), 'Limiting multiple alignments.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

a <- as.data.frame(table(anchorReadAlignments$qName.anchorReads))
names(a) <- c('readID', 'anchorReadAlnFreq')
b <- as.data.frame(table(adriftReadAlignments$qName.adriftReads))
names(b) <- c('readID', 'adriftReadAlnFreq')
tab <- left_join(a, b, by = 'readID')

tab$anchorReadAlnFreq <- ifelse(tab$anchorReadAlnFreq >= opt$buildFragments_min_num_alignments_for_multiHit, 1, 0)
tab$adriftReadAlnFreq <- ifelse(tab$adriftReadAlnFreq >= opt$buildFragments_min_num_alignments_for_multiHit, 1, 0)
tab$s <- rowSums(tab[ , c(2,3)], na.rm=TRUE)

if(! opt$buildFragments_readMates_allowed_to_be_multiHit %in% c('none', 'one', 'both')){
  write(c(paste(now(), 'Errror -- buildFragments_readMates_allowed_to_be_multiHit must be set on none, one, or both.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}

if(opt$buildFragments_readMates_allowed_to_be_multiHit == 'none'){
  n <- 0
} else if (opt$buildFragments_readMates_allowed_to_be_multiHit == 'one'){
  n <- 1
} else {
  n <- 2
}

readsToRemove <- tab[tab$s > n,]$readID

anchorReadAlignments <- subset(anchorReadAlignments, ! qName.anchorReads %in% readsToRemove)
adriftReadAlignments <- subset(adriftReadAlignments, ! qName.adriftReads %in% readsToRemove)

write(c(paste(now(), sprintf("%.2f%%", (n_distinct(readsToRemove) / n_distinct(tab$readID))*100), ' reads removed due to multiple alignments.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

ids <- unique(anchorReadAlignments$qName.anchorReads)
id_groups <- split(ids, dplyr::ntile(1:length(ids), ceiling(length(ids)/opt$buildFragments_idGroup_size)))

write(c(paste(now(), 'Building initial fragments.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

frags <- bind_rows(lapply(id_groups, function(id_group){
  library(dplyr)
  
  a <- subset(anchorReadAlignments, qName.anchorReads %in% id_group)
  b <- subset(adriftReadAlignments, qName.adriftReads %in% id_group)
  
  if(nrow(a) == 0 | nrow(b) == 0) return(data.frame())
  
  # Join adrift reads alignments to anchor read alignments to create potential read pairs.
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
  r <- mutate(frags, 
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
   r
}))

# The join between anchor and adrift reads can assign a single read to multiple fragments.

# save.image('~/buildFragments1.RData')

# Row duplication is occuring somewhere
frags <- dplyr::distinct(frags)


# Random ids.
if(opt$demultiplex_captureRandomLinkerSeq){
  write(c(paste(now(), 'Capturing random IDs from adrift read linkers.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
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
  # too few reads were associated with the random id to attempt to segregate. 
  z <- subset(o, reads < opt$buildFragments_randomLinkerID_minReadCountToSegreagate)$randomLinkerSeq
  if(length(z) > 0){
    frags <- subset(frags, ! randomLinkerSeq %in% z)
  }
  
  
  # Define random ids with two or more samples and enough reads to try to salvage,
  # ie. retain reads where most are assigned to a single sample.
  
  z2 <- subset(o, reads >= opt$buildFragments_randomLinkerID_minReadCountToSegreagate)$randomLinkerSeq
  
  if(length(z2) > 0){
    a <- subset(frags, ! randomLinkerSeq %in% z2)  # Undisputed reads
    b <- subset(frags, randomLinkerSeq %in% z2)    # Reads which need to be examined.
    
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

  write(c(paste(now(), msg)), file = file.path(opt$outputDir, 'log'), append = TRUE)
}

gc()

frags <- unpackUniqueSampleID(frags)

samples <- loadSamples()
frags <- left_join(frags, select(samples, uniqueSample, flags), by = 'uniqueSample')

dups <- tibble()

if('buildFragments_duplicateReadFile' %in% names(opt)){
  write(c(paste(now(), 'Reading duplicate read file created by prepReads.R.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  dups <- readRDS(file.path(opt$outputDir, opt$buildFragments_duplicateReadFile))
  dups <- data.table(dplyr::distinct(dplyr::select(dups, id, n)))
}

cluster <- makeCluster(opt$buildFragments_CPUs)
clusterExport(cluster, c('opt', 'dups'))

frags$fragID <- paste0(frags$trial, ':', frags$subject, ':', frags$sample, ':', frags$replicate, ':', frags$chromosome, ':', frags$strand, ':', frags$fragStart, ':', frags$fragEnd)

# Combine read level fragments into unique fragments with read counts.
# Duplicate reads are tallied if a duplicate read file is provided.
# Using data.table objects since they are split more efficiently.

f <- data.table(frags)

write(c(paste(now(), 'Refining fragments.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

frags <- parLapply(cluster, split(f, f$fragID), function(x){
           library(dplyr)
           library(data.table)
           source(file.path(opt$softwareDir, 'lib.R'))

           if(nrow(x) == 1){
             x$reads <- 1
             x$repLeaderSeq <- x$leaderSeq.anchorReads
             x$readIDlist   <- x$readID
             return(select(as_tibble(x), -readID, -leaderSeq.anchorReads))
           }

           r <- representativeSeq(x$leaderSeq.anchorReads)

           # Exclude reads where the leaderSeq is not similar to the representative sequence.
           i <- stringdist::stringdist(r[[2]], x$leaderSeq.anchorReads) / nchar(r[[2]]) <= opt$buildStdFragments_maxLeaderSeqDiffScore
           if(all(! i)) return(data.table::data.table())

           x <- x[i,]

           x$reads <- nrow(x)
           x$repLeaderSeq <- r[[2]]

           o <- dups[dups$id %in% x$readID]
           if(nrow(o) > 0){
             x$reads <- nrow(x) + sum(o$n)
           }
           
           # Duplicate reads IDs are not stored.
           readIDs <- x$readID
           x <- as_tibble(x[1,])
           x$readIDlist <- list(readIDs)

           return(select(x[1,], -readID, -leaderSeq.anchorReads))
})

# Force single read ids into a list.
frags <- bind_rows(lapply(frags, function(x){
           if(! is.list(x$readIDlist)) x$readIDlist <- list(x$readIDlist)
           x
        }))

stopCluster(cluster)

frags <- dplyr::select(frags, trial, subject, sample, replicate, chromosome, strand, 
                       fragStart, fragEnd, reads, repLeaderSeq, flags, readIDlist)

frags$compDate <- as.character(lubridate::today())
frags$dataLabel <- opt$dataLabel
frags$uniqueSample <- paste0(frags$trial, '~', frags$subject,  '~', frags$sample,  '~', frags$replicate)
frags <- left_join(frags, dplyr::select(samples, uniqueSample, refGenome.id), by = 'uniqueSample')

saveRDS(distinct(frags), file.path(opt$outputDir, opt$buildFragments_outputDir, opt$buildFragments_outputFile), compress = 'xz')

if('databaseGroup' %in% names(opt)){
  library(RMariaDB)
  library(DBI)
  con <- dbConnect(RMariaDB::MariaDB(), group = opt$databaseGroup)
  
  invisible(lapply(split(frags, paste(frags$uniqueSample, frags$refGenome.id)), function(x){

    dbExecute(con, paste0("delete from fragments where trial='", x$trial[1], "' and subject='", x$subject[1], 
           "' and sample='", x$sample[1], "' and replicate='", x$replicate[1], "' and refGenome='", x$refGenome.id[1], "'"))
    
    r <- dbExecute(con,
              "insert into fragments values (?, ?, ?, ?, ?, ?, ?, ?)",
              params = list(x$trial[1], x$subject[1], x$sample[1], x$replicate[1], x$refGenome.id[1],
                            x$dataLabel[1], x$compDate[1], list(serialize(dplyr::select(x, -uniqueSample, -refGenome.id, -dataLabel, -compDate), NULL))))
    
    if(r == 0){
        write(c(paste(now(), 'Errror -- could not upload fragment data for ', x$uniqueSample[1], ' to the database.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
        q(save = 'no', status = 1, runLast = FALSE) 
      }
  }))
  
  dbDisconnect(con)
}

q(save = 'no', status = 0, runLast = FALSE) 
