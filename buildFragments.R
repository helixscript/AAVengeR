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
write(c(paste(now(), '   Reading in anchor and adrift read alignments.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
adriftReadAlignments <- readRDS(file.path(opt$outputDir, opt$buildFragments_adriftReadsAlignmentFile))
anchorReadAlignments <- readRDS(file.path(opt$outputDir, opt$buildFragments_anchorReadsAlignmentFile))

anchorReadAlignments$posid <- paste0(anchorReadAlignments$tName, anchorReadAlignments$strand, ifelse(anchorReadAlignments$strand == '+', anchorReadAlignments$tStart, anchorReadAlignments$tEnd))

r <- dplyr::select(readRDS(file.path(opt$outputDir, opt$prepReads_outputDir, 'reads.rds')), readID, adriftReadRandomID) %>% dplyr::rename(randomLinkerSeq = adriftReadRandomID) %>% data.table()
  
# Limit random ids to those found in the adrift read alignments.
r <- subset(r, readID %in% adriftReadAlignments$readID)

anchorReadAlignments$sample <- sub('~\\d+$', '', anchorReadAlignments$uniqueSample)
adriftReadAlignments$sample <- sub('~\\d+$', '', adriftReadAlignments$uniqueSample)
r <- left_join(r, distinct(dplyr::select(adriftReadAlignments, readID, sample)), by = 'readID')  

  
  
# Correct random ids to the most abundant within samples.
# Uncorrectable codes are returned as NNNN and removed.

write(c(paste(now(), '   Correcting minor differences in random ids.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

# anchorReadAlignments too large to export to worker nodes.
# stringDist uses all cores regardless of nthread option... best not to use parLappy().
ara <- dplyr::select(anchorReadAlignments, readID, posid)

# nrow(subset(anchorReadAlignments, posid == 'chr2+202308587')); nrow(subset(anchorReadAlignments, posid == 'chr2-202308470'))

r <- rbindlist(lapply(split(r, r$sample), function(x){
         message(x$sample[1])
         x$randomLinkerSeqPre <- x$randomLinkerSeq
         x$randomLinkerSeq    <- conformMinorSeqDiffs(x$randomLinkerSeq, nThreads = opt$buildFragments_CPUs)
         
         o <- x[grepl('N', x$randomLinkerSeq),]
         if(nrow(o) > 0){
           o <- dplyr::select(dplyr::mutate(o, randomLinkerSeq = randomLinkerSeqPre), -randomLinkerSeqPre)
           o <- left_join(o, ara, by = 'readID')
           readr::write_tsv(o, file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'randomIDexcludedReads', paste0('uncorrectableIDs~', x$sample[1], '.tsv')))
         }     
         
         data.table(x[! grepl('N', x$randomLinkerSeq),])
       }))
  
write(c(paste(now(), '   Determining which random ids span multiple samples.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

o <- group_by(r, randomLinkerSeq) %>%
       summarise(nSamples = n_distinct(sample)) %>%
       ungroup() %>%
       filter(nSamples > 1)

write(c(paste(now(), paste0('   ', sprintf("%.2f%%", (n_distinct(o$randomLinkerSeq)/n_distinct(r$randomLinkerSeq))*100), ' random ids seen across two or more samples'))), file = file.path(opt$outputDir, 'log'), append = TRUE)
  

# Need to remove ids without enough reads to resolve...

randomIDsNoIssue   <- dplyr::filter(r, ! randomLinkerSeq %in% o$randomLinkerSeq)

randomIDsToResolve <- dplyr::filter(r, randomLinkerSeq %in% o$randomLinkerSeq) %>% 
                      dplyr::group_by(randomLinkerSeq) %>% 
                      dplyr::mutate(n = n()) %>% 
                      dplyr::ungroup() %>% 
                      dplyr::filter(n >= opt$buildFragments_randomLinkerID_minReadCountToSegreagate) %>%
                      dplyr::pull(unique(randomLinkerSeq))

if(length(randomIDsToResolve) > 0){
   # Split data frame into CPU chunks while not breaking appart random ids.
   r2 <- subset(r, randomLinkerSeq %in% randomIDsToResolve)
   
   randomIDIssuesResolved <- rbindlist(lapply(split(r2, r2$randomLinkerSeq), function(y){
      y$remove <- TRUE
          
      t <- sort(table(y$sample), decreasing = TRUE)
             if((t[1]/sum(t))*100 >= opt$buildFragments_randomLinkerID_minSingleSampleMajorityPercent){
                 y[y$sample == names(t)[1],]$remove <- FALSE
             }
       
             y
   })) %>% dplyr::filter(remove == FALSE) %>% dplyr::select(-remove)
   
   o <- subset(r2, ! randomLinkerSeq %in% randomIDIssuesResolved$randomLinkerSeq)
   if(nrow(o) > 0){
     o <- left_join(o, distinct(dplyr::select(anchorReadAlignments, readID, posid)), by = 'readID')
     readr::write_tsv(o, file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'randomIDexcludedReads', 'read_IDs_shared_between_samples.tsv'))
   }
     
   r <- bind_rows(randomIDsNoIssue, randomIDIssuesResolved)
}
  

# Remove reads where random IDs were in conflict and could not be corrected.
adriftReadAlignments <- subset(adriftReadAlignments, readID %in% r$readID)
anchorReadAlignments <- subset(anchorReadAlignments, readID %in% r$readID)
  
adriftReadAlignments <- left_join(adriftReadAlignments, distinct(dplyr::select(r, readID, randomLinkerSeq)), by = 'readID')
rm(r, randomIDs)
gc()


write(c(paste(now(), '   Preparing alignment data for fragment generation.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

anchorReadAlignments <- subset(anchorReadAlignments, readID %in% adriftReadAlignments$readID)

dir.create(file.path(opt$outputDir, opt$buildFragments_outputDir))

anchorReadAlignments <- select(anchorReadAlignments, uniqueSample, sample, readID, tName, strand, tStart, tEnd, leaderSeq)
adriftReadAlignments <- select(adriftReadAlignments, sample, readID, tName, strand, tStart, tEnd, randomLinkerSeq)

names(anchorReadAlignments) <- paste0(names(anchorReadAlignments), '.anchorReads')
names(adriftReadAlignments) <- paste0(names(adriftReadAlignments), '.adriftReads')

# Breaking reads into id groups prevents joins from exceeding internal row limits
# before they can be filtered for correct pairings.

ids <- unique(anchorReadAlignments$readID.anchorReads)
id_groups <- split(ids, dplyr::ntile(1:length(ids), ceiling(length(ids)/opt$buildFragments_idGroup_size)))

write(c(paste(now(), '   Building initial fragments.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

# Convert from tibble to data.table for increased subsetting efficiency. 
anchorReadAlignments <- data.table(anchorReadAlignments)
adriftReadAlignments <- data.table(adriftReadAlignments)



# Identify read pairs where both mates have many alignments that may cause the system 
# to run out of memory during anchor and adrift read joins.

a <- as.data.frame(table(anchorReadAlignments$readID.anchorReads))
names(a) <- c('readID', 'anchorReadAlnFreq')
b <- as.data.frame(table(adriftReadAlignments$readID.adriftReads))
names(b) <- c('readID', 'adriftReadAlnFreq')
tab <- left_join(a, b, by = 'readID')

tab$anchorReadAlnFreq <- ifelse(tab$anchorReadAlnFreq > opt$buildFragments_maxReadAlignments, 1, 0)
tab$adriftReadAlnFreq <- ifelse(tab$adriftReadAlnFreq > opt$buildFragments_maxReadAlignments, 1, 0)
tab$s <- rowSums(tab[ , c(2,3)], na.rm=TRUE)

# Find read pairs where both anchor and adrift mates have more than opt$buildFragments_maxReadAlignments alignments.
z <- as.character(tab[tab$s == 2,]$readID)

# For read pairs where both mates have too many alignments, randomly select anchor reads and associated adrift mates 
# that have the potential to form rational fragments.

write(paste0(now(), '    Building rational fragments from alignment data.'), file = file.path(opt$outputDir, 'log'), append = TRUE)

if(length(z) > 0){
  write(c(paste0(now(), '    ', length(z), ' reads pairs have more than ', opt$buildFragments_maxReadAlignments, ' alignments for both mates.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  write(c(paste0(now(), '    For each reach read, ', opt$buildFragments_maxReadAlignments, ' alignments will be randomly selected.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  write(c(paste0(now(), '    Adrift alignments that have the potential to form rational fragments with the selected anchor reads will be selected as well.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
  a1 <- subset(anchorReadAlignments, ! readID.anchorReads %in% z)
  b1 <- subset(adriftReadAlignments, ! readID.adriftReads %in% z)
  
  a2 <- subset(anchorReadAlignments, readID.anchorReads %in% z)
  b2 <- subset(adriftReadAlignments, readID.adriftReads %in% z)
  
  o <- lapply(split(a2, a2$readID.anchorReads), function(x){
          set.seed(1)
          x <- dplyr::sample_n(x, opt$buildFragments_maxReadAlignments)
          y <- subset(b2, readID.adriftReads %in% x$readID.anchorReads)  # find corresponding sampled ids in adrift reads.
          
          expand <- 1000
          y.pos <- tibble()
          y.neg <- tibble()
          
          if(nrow(y)){
             xpos <- subset(x, strand.anchorReads == '+')
             xneg <- subset(x, strand.anchorReads == '-')
          
             if(nrow(xpos) > 0){
                y$start <- y$tEnd.adriftReads - opt$buildFragments_maxFragLength
                y$end   <- y$tEnd.adriftReads
               
                z <- subset(y, strand.adriftReads == '-')
                g1 <- GenomicRanges::makeGRangesFromDataFrame(xpos, seqnames.field = 'tName.anchorReads', start.field = 'tStart.anchorReads', end.field = 'tStart.anchorReads', strand.field = 'strand.anchorReads')
                g2 <- GenomicRanges::makeGRangesFromDataFrame(z, seqnames.field = 'tName.adriftReads', start.field = 'start', end.field = 'end', strand.field = 'strand.adriftReads')
                o <- suppressWarnings(GenomicRanges::findOverlaps(g1, g2, ignore.strand = TRUE))
                if(length(o) > 0) y.pos <- z[subjectHits(o),]
             }
             
             if(nrow(xneg) > 0){
               y$start <- y$tStart.adriftReads
               y$end   <- y$tStart.adriftReads + opt$buildFragments_maxFragLength
               
               z <- subset(y, strand.adriftReads == '+')
               g1 <- GenomicRanges::makeGRangesFromDataFrame(xneg, seqnames.field = 'tName.anchorReads', start.field = 'tStart.anchorReads', end.field = 'tStart.anchorReads', strand.field = 'strand.anchorReads')
               g2 <- GenomicRanges::makeGRangesFromDataFrame(z, seqnames.field = 'tName.adriftReads', start.field = 'start', end.field = 'end', strand.field = 'strand.adriftReads')
               o <- suppressWarnings(GenomicRanges::findOverlaps(g1, g2, ignore.strand = TRUE))
               if(length(o) > 0) y.neg <- z[subjectHits(o),]
             }
          }
          
          y2 <- bind_rows(y.pos, y.neg)
          list(x, y2)
  }) 
  
  a3 <- bind_rows(lapply(o, '[[', 1))
  b3 <- bind_rows(lapply(o, '[[', 2))
  
  anchorReadAlignments <- bind_rows(a1, a3)
  adriftReadAlignments <- bind_rows(b1, b3)
}



# (dev) parallelize this and use rbindlist
# Large memory usage here.
# %in% slow -- consider dplyr group indicies solution.

o <- lapply(id_groups, function(id_group){
       list(anchorReadAlignments[readID.anchorReads %in% id_group],
            adriftReadAlignments[readID.adriftReads %in% id_group])
     })

rm(anchorReadAlignments, adriftReadAlignments)
gc()

counter <- 1
total <- length(o)

frags <- bind_rows(lapply(o, function(z){
  a <- z[[1]]
  b <- z[[2]]
  
  if(nrow(a) == 0 | nrow(b) == 0) return(data.frame())
  
  write(paste0(now(), '    ', counter, '/', total, ': ', nrow(a), ' anchorRead alignments, ', nrow(b), ' adriftRead_alignments.'), file = file.path(opt$outputDir, 'log'), append = TRUE)
  counter <<- counter + 1
  
  # Join adrift reads alignments to anchor read alignments to create potential read pairs.
  frags <- left_join(a, b, by = c('readID.anchorReads' = 'readID.adriftReads')) %>% tidyr::drop_na()
  
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
                mutate(uniqueSample = uniqueSample.anchorReads, readID = readID.anchorReads) %>%
                select(uniqueSample, readID, chromosome, strand, fragStart, fragEnd, leaderSeq.anchorReads, randomLinkerSeq.adriftReads)
   r
}))

rm(o, id_groups)
gc()

write(c(paste(now(), '   Fragment generation complete.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

if('buildFragments_duplicateReadFile' %in% names(opt)){
  write(c(paste(now(), '   Reading duplicate read file created by prepReads.R.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  dups <- readRDS(file.path(opt$outputDir, opt$buildFragments_duplicateReadFile))
  dups <- data.table(dplyr::distinct(dplyr::select(dups, id, n)))
  dups <- subset(dups, id %in% frags$readID)
  
  frags <- left_join(frags, dups, by = c('readID' = 'id'))
  frags$n <- ifelse(is.na(frags$n), 0, frags$n)
} else {
  frags$n <- 0
}


write(c(paste(now(), '   Adding sample details to fragment data.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
samples <- loadSamples()
frags <- left_join(frags, distinct(select(samples, uniqueSample, refGenome.id, vectorFastaFile, flags)), by = 'uniqueSample')

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
