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

dir.create(file.path(opt$outputDir, opt$buildFragments_outputDir))
dir.create(file.path(opt$outputDir, opt$buildFragments_outputDir, 'randomIDexcludedReads'))

# Read in adrift alignment reads.
adriftReadAlignments <- readRDS(file.path(opt$outputDir, opt$buildFragments_adriftReadsAlignmentFile))

cluster <- makeCluster(opt$buildFragments_CPUs)
clusterExport(cluster, c('opt'))

# Random ids.
if(opt$demultiplex_captureRandomLinkerSeq){
  
  # Collect the random IDs isolated during demultiplexing.
  files <- list.files(file.path(opt$outputDir, opt$demultiplex_outputDir), pattern = 'randomIDs', full.names = TRUE)
  randomIDs <- Reduce('append', lapply(files,  readDNAStringSet))
  r <- tibble(readID = names(randomIDs), randomLinkerSeq = as.character(randomIDs))
  
  # Limit random ids to those found in the adrift read alignments.
  r <- subset(r, readID %in% adriftReadAlignments$qName)
  message(n_distinct(r$randomLinkerSeq), ' unqiue random ids found in adrift read alignments.')
  
  # Add sample ids to read alignments.
  message('Adding sample names to adrift alignments')
  r <- dplyr::distinct(left_join(r, dplyr::select(adriftReadAlignments, qName, sample), by = c('readID' = 'qName')))
  
  # Remove replicate identifiers from sample names to evaluate on the sample level.
  r$sample <- sub('~\\d+$', '', r$sample)
  
  # Correct random ids to the most abundant within samples.
  # Uncorrectable codes are returned as NNNN and removed.
  message('Correcting minor differences in random ids.')
  #r <- bind_rows(lapply(split(r, r$sample), function(x){
  r <- bind_rows(parLapply(cluster, split(r, r$sample), function(x){
         source(file.path(opt$softwareDir, 'lib.R'))
         library(dplyr)
         x$randomLinkerSeqPre <- x$randomLinkerSeq
         x$randomLinkerSeq    <- conformMinorSeqDiffs(x$randomLinkerSeq)
         
         o <- x[grepl('N', x$randomLinkerSeq),]
         if(nrow(o) > 0) readr::write_tsv(subset(x, readID %in% o$readID), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'randomIDexcludedReads', paste0('uncorrectableIDs~', x$sample[1], '.tsv')))
         
         x[! grepl('N', x$randomLinkerSeq),]
       }))
  
  message('Determining which random ids span multiple samples.')
  o <- group_by(r, randomLinkerSeq) %>%
       summarise(nSamples = n_distinct(sample)) %>%
       ungroup() %>%
       filter(nSamples > 1)
  
  r$remove <- FALSE
  
  message(sprintf("%.2f%%", (n_distinct(o$randomLinkerSeq)/n_distinct(r$randomLinkerSeq))*100), ' random ids seen across two or more samples')
  
  if(nrow(o) > 0){
   message('Cleaning up instances where random ids are seen across samples')
   invisible(lapply(split(o, 1:nrow(o)), function(a){
     b <- subset(r, randomLinkerSeq == a$randomLinkerSeq)
     
     if(nrow(b) >= opt$buildFragments_randomLinkerID_minReadCountToSegreagate){
       # Yes, there are enough reads to attempt to segregate.
       t <- sort(table(b$sample), decreasing = TRUE)
       
       if((t[1]/sum(t))*100 >= opt$buildFragments_randomLinkerID_minSingleSampleMajorityPercent){
         # Yes, there is a predominant sample for this code.
         r[r$randomLinkerSeq == a$randomLinkerSeq & ! r$sample == names(t[1]),]$remove <<- TRUE
       } else {
         # No, no predominant sample, remove all codes.
         r[r$randomLinkerSeq == a$randomLinkerSeq,]$remove <<- TRUE
       }
     } else {
       r[r$randomLinkerSeq == a$randomLinkerSeq,]$remove <<- TRUE
     } 
   })) 
  }
  
  o <- subset(r, remove == TRUE)
  if(nrow(o) > 0)  readr::write_tsv(subset(r, remove == TRUE), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'randomIDexcludedReads', 'read_IDs_shared_between_samples.tsv'))
    
  saveRDS(subset(r, remove == TRUE), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'readsRemoved_randomID_sample_conflicts.rds'))
  
  r <- subset(r, remove != TRUE)
  
  adriftReadAlignments <- subset(adriftReadAlignments, qName %in% r$readID)
  adriftReadAlignments <- left_join(adriftReadAlignments, 
                                    distinct(dplyr::select(r, readID, randomLinkerSeq)), by = c('qName' = 'readID'))
  rm(r, randomIDs)
  gc()
  
} else {
  adriftReadAlignments$randomLinkerSeq <- 'NNNNNNNNNNNN'
}


write(c(paste(now(), 'Reading in alignment results.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

anchorReadAlignments <- readRDS(file.path(opt$outputDir, opt$buildFragments_anchorReadsAlignmentFile))

readID2reference <- distinct(tibble(readID = anchorReadAlignments$qName, refGenome = anchorReadAlignments$refGenome.id))

anchorReadAlignments <- subset(anchorReadAlignments, qName %in% adriftReadAlignments$qName)

dir.create(file.path(opt$outputDir, opt$buildFragments_outputDir))

anchorReadAlignments <- select(anchorReadAlignments, sample, qName, tName, strand, tStart, tEnd, leaderSeq)
adriftReadAlignments <- select(adriftReadAlignments, sample, qName, tName, strand, tStart, tEnd, randomLinkerSeq)

names(anchorReadAlignments) <- paste0(names(anchorReadAlignments), '.anchorReads')
names(adriftReadAlignments) <- paste0(names(adriftReadAlignments), '.adriftReads')

#----------------------------------

# Create a table where each row is a read pair id and the columns are the number
# of anchor and adrift read alignments. These values are converted to 0 if the 
# value is less than opt$buildFragments_maxReadAlignments and converted to 1 if 
# the value is greater. Row sums of 2 tell us that both mates have more than 
# alignments and should be pruned.

a <- as.data.frame(table(anchorReadAlignments$qName.anchorReads))
names(a) <- c('readID', 'anchorReadAlnFreq')
b <- as.data.frame(table(adriftReadAlignments$qName.adriftReads))
names(b) <- c('readID', 'adriftReadAlnFreq')
tab <- left_join(a, b, by = 'readID')

tab$anchorReadAlnFreq <- ifelse(tab$anchorReadAlnFreq > opt$buildFragments_maxReadAlignments, 1, 0)
tab$adriftReadAlnFreq <- ifelse(tab$adriftReadAlnFreq > opt$buildFragments_maxReadAlignments, 1, 0)
tab$s <- rowSums(tab[ , c(2,3)], na.rm=TRUE)

# Find read pairs where both anchor and adrift mates have more than opt$buildFragments_maxReadAlignments alignments.
z <- as.character(tab[tab$s == 2,]$readID)

if(length(z) > 0){
  
  write(c(paste(now(), paste0(length(z), ' read pairs, ', sprintf("%.2f%%", (length(z)/nrow(tab))*100),    
                              ' of reads, have more than ', 
                              opt$buildFragments_maxReadAlignments, ' alignments.'))), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
  # Find the anchor read alignments from reads with more than the maximum number of alignments.
  a <- filter(anchorReadAlignments, qName.anchorReads %in% z)
  
  # Sample a number of alignments from each overly aligned anchor read.
  set.seed(1)
  a <- group_by(a, qName.anchorReads) %>%
       sample_n(opt$buildFragments_maxReadAlignments) %>%
       ungroup()
  
  # Find the subset of adrift reads alignments which corespond to the randomly selected anchor read alignments.
  b <- filter(adriftReadAlignments, qName.adriftReads %in% a$qName.anchorReads)
  
  # Split the alignments by read id.
  a <- split(a, a$qName.anchorReads)
  b <- split(b, b$qName.adriftReads)
  if(! all(names(a) == names(b))) stop('Error -- excess read name alignment error.')
  
  # Split the anchor read and adrift reads lists into another list whose elements can be processed in parallel.
  k <- lapply(split(1:length(a), ntile(1:length(a), opt$buildFragments_CPUs)), function(x) list(a[x], b[x]))

  o <- parLapply(cluster, k, function(x){
  #o <- lapply(k, function(x){
         library(dplyr)
         library(GenomicRanges)
    
         a <- x[[1]]
         b <- x[[2]]

         j <- lapply(1:length(a), function(x2){
                a <- a[[x2]]
                b <- b[[x2]]
         
                a.pos <- subset(a, strand.anchorReads == '+')
                a.neg <- subset(a, strand.anchorReads == '-')
         
                a.pos.gr <- GenomicRanges::makeGRangesFromDataFrame(tibble(
                              seqnames = a.pos$tName.anchorReads, 
                              strand = '+',
                              start = a.pos$tStart.anchorReads,
                              end = a.pos$tEnd.anchorReads))
         
                a.neg.gr <- GenomicRanges::makeGRangesFromDataFrame(tibble(
                              seqnames = a.neg$tName.anchorReads, 
                              strand = '-',
                              start = a.neg$tStart.anchorReads,
                              end = a.neg$tEnd.anchorReads))
                
                b.pos <- subset(b, b$strand.adriftReads == '+')
                b.neg <- subset(b, b$strand.adriftReads == '-')
         
                b.pos.gr <- GenomicRanges::makeGRangesFromDataFrame(tibble(
                              seqnames = b.pos$tName.adriftReads, 
                              strand = '+',
                              start = b.pos$tStart.adriftReads-1000,
                              end = b.pos$tEnd.adriftReads+1000))
         
               b.neg.gr <- GenomicRanges::makeGRangesFromDataFrame(tibble(
                             seqnames = b.neg$tName.adriftReads, 
                             strand = '-',
                             start = b.neg$tStart.adriftReads-1000,
                             end = b.neg$tEnd.adriftReads+1000))
         
               negAlignmentsToReturn <- tibble()
               posAlignmentsToReturn <- tibble()
         
               i <- vector()
               if(length(a.pos.gr) > 0 & length(b.neg.gr) > 0) i <- findOverlaps(a.pos.gr, b.neg.gr, ignore.strand = TRUE)
               if(length(i) > 0) negAlignmentsToReturn <- b.neg[unique(as.integer(subjectHits(i))),]
         
               i <- vector()
               if(length(a.neg.gr) > 0 & length(b.pos.gr) > 0) i <- findOverlaps(a.neg.gr, b.pos.gr, ignore.strand = TRUE)
               if(length(i) > 0) posAlignmentsToReturn <- b.pos[unique(as.integer(subjectHits(i))),]
         
              list(bind_rows(a.pos, a.neg), bind_rows(negAlignmentsToReturn, posAlignmentsToReturn))
          })
  
         list(bind_rows(lapply(j, '[[', 1)), (bind_rows(lapply(j, '[[', 2)))) 
       })
  
  # Recombine the selected alignments from the parallel processing.
  o2 <- list(bind_rows(lapply(o, '[[', 1)), (bind_rows(lapply(o, '[[', 2)))) 
  
  # Remove alignments with too many alignments from the alignment data frames.
  anchorReadAlignments <- subset(anchorReadAlignments, ! qName.anchorReads %in% z)
  adriftReadAlignments <- subset(adriftReadAlignments, ! qName.adriftReads %in% z)
  
  # Add back the sampled alignments.
  anchorReadAlignments <- bind_rows(anchorReadAlignments, o2[[1]])
  adriftReadAlignments <- bind_rows(adriftReadAlignments, o2[[2]])
  
  rm(z, a, b, tab, k, o, o2)
  gc()
}

ids <- unique(anchorReadAlignments$qName.anchorReads)
id_groups <- split(ids, dplyr::ntile(1:length(ids), ceiling(length(ids)/opt$buildFragments_idGroup_size)))

write(c(paste(now(), 'Building initial fragments.')), file = file.path(opt$outputDir, 'log'), append = TRUE)


# Convert from tibble to data.table for increased subsetting efficiency. 
anchorReadAlignments <- data.table(anchorReadAlignments)
adriftReadAlignments <- data.table(adriftReadAlignments)


o <- lapply(id_groups, function(id_group){
       list(anchorReadAlignments[qName.anchorReads %in% id_group],
            adriftReadAlignments[qName.adriftReads %in% id_group])
     })

rm(anchorReadAlignments, adriftReadAlignments)
gc()


frags <- bind_rows(lapply(o, function(z){
  a <- z[[1]]
  b <- z[[2]]
  
  if(nrow(a) == 0 | nrow(b) == 0) return(data.frame())
  
  # Consider using data.table nomenclature for joins.
  
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
                select(uniqueSample, readID, chromosome, strand, fragStart, fragEnd, leaderSeq.anchorReads, randomLinkerSeq.adriftReads)
   r
}))

rm(o, id_groups)
gc()


if('buildFragments_duplicateReadFile' %in% names(opt)){
  write(c(paste(now(), 'Reading duplicate read file created by prepReads.R.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  dups <- readRDS(file.path(opt$outputDir, opt$buildFragments_duplicateReadFile))
  dups <- data.table(dplyr::distinct(dplyr::select(dups, id, n)))
  dups <- subset(dups, id %in% frags$readID)
  
  frags <- left_join(frags, dups, by = c('readID' = 'id'))
  frags$n <- ifelse(is.na(frags$n), 0, frags$n)
} else {
  frags$n <- 0
}

stopCluster(cluster)

samples <- loadSamples()
frags <- left_join(frags, select(samples, uniqueSample, flags), by = 'uniqueSample')
frags <- left_join(frags, readID2reference, by = 'readID')

saveRDS(frags, file.path(opt$outputDir, opt$buildFragments_outputDir, opt$buildFragments_outputFile))

q(save = 'no', status = 0, runLast = FALSE) 



#--------------------------------------------------------------------------------










# # Add sample flags.
# # samples <- loadSamples()
# # frags <- left_join(frags, select(samples, uniqueSample, flags), by = 'uniqueSample')
# 
# 
# frags$fragID <- paste0(frags$trial, ':', frags$subject, ':', frags$sample, ':', frags$replicate, ':', frags$chromosome, ':', 
#                        frags$strand, ':', frags$fragStart, ':', frags$fragEnd, ':', frags$randomLinkerSeq.adriftReads)
# 
# 
# stopCluster(cluster)
# 
# frags <- dplyr::select(frags, trial, subject, sample, replicate, chromosome, strand, 
#                        fragStart, fragEnd, reads, repLeaderSeq, flags, readIDlist)
# 
# frags$compDate <- as.character(lubridate::today())
# frags$dataLabel <- opt$dataLabel
# frags$uniqueSample <- paste0(frags$trial, '~', frags$subject,  '~', frags$sample,  '~', frags$replicate)
# frags <- left_join(frags, dplyr::select(samples, uniqueSample, refGenome.id), by = 'uniqueSample')
# 
# saveRDS(distinct(frags), file.path(opt$outputDir, opt$buildFragments_outputDir, opt$buildFragments_outputFile), compress = 'xz')
# 
# if('databaseGroup' %in% names(opt)){
#   library(RMariaDB)
#   library(DBI)
#   con <- dbConnect(RMariaDB::MariaDB(), group = opt$databaseGroup)
#   
#   invisible(lapply(split(frags, paste(frags$uniqueSample, frags$refGenome.id)), function(x){
# 
#     dbExecute(con, paste0("delete from fragments where trial='", x$trial[1], "' and subject='", x$subject[1], 
#            "' and sample='", x$sample[1], "' and replicate='", x$replicate[1], "' and refGenome='", x$refGenome.id[1], "'"))
#     
#     r <- dbExecute(con,
#               "insert into fragments values (?, ?, ?, ?, ?, ?, ?, ?)",
#               params = list(x$trial[1], x$subject[1], x$sample[1], x$replicate[1], x$refGenome.id[1],
#                             x$dataLabel[1], x$compDate[1], list(serialize(dplyr::select(x, -uniqueSample, -refGenome.id, -dataLabel, -compDate), NULL))))
#     
#     if(r == 0){
#         write(c(paste(now(), 'Errror -- could not upload fragment data for ', x$uniqueSample[1], ' to the database.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
#         q(save = 'no', status = 1, runLast = FALSE) 
#       }
#   }))
#   
#   dbDisconnect(con)
# }
# 
# q(save = 'no', status = 0, runLast = FALSE) 
