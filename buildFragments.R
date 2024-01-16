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
setMissingOptions()
setOptimalParameters()
set.seed(1)

dir.create(file.path(opt$outputDir, opt$buildFragments_outputDir))
dir.create(file.path(opt$outputDir, opt$buildFragments_outputDir, 'tmp'))
dir.create(file.path(opt$outputDir, opt$buildFragments_outputDir, 'randomIDexcludedReads'))

# Read in adrift alignment reads.
write(c(paste(now(), '   Reading in anchor and adrift read alignments.')), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'log'), append = FALSE)
adriftReadAlignments <- readRDS(file.path(opt$outputDir, opt$buildFragments_adriftReadsAlignmentFile))
anchorReadAlignments <- readRDS(file.path(opt$outputDir, opt$buildFragments_anchorReadsAlignmentFile))

if(nrow(anchorReadAlignments) == 0 | nrow(adriftReadAlignments) == 0){
  write(c(paste(now(), '   Error - anchor and/or adrift alignment data files were empty.')), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'log'), append = TRUE)
  if(opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()
  q(save = 'no', status = 1, runLast = FALSE) 
}

incomingSamples <- unique(anchorReadAlignments$uniqueSample)

# Shorten file paths to save memory since we no longer need the full paths.
anchorReadAlignments$refGenome <- sapply(anchorReadAlignments$refGenome, lpe)
anchorReadAlignments$vectorFastaFile <- sapply(anchorReadAlignments$vectorFastaFile, lpe)


# Create posid column.
anchorReadAlignments$posid <- paste0(anchorReadAlignments$tName, anchorReadAlignments$strand, ifelse(anchorReadAlignments$strand == '+', anchorReadAlignments$tStart, anchorReadAlignments$tEnd))

anchorReadAlignments$sample <- sub('~\\d+$', '', anchorReadAlignments$uniqueSample)
adriftReadAlignments$sample <- sub('~\\d+$', '', adriftReadAlignments$uniqueSample)

r <- dplyr::select(readRDS(file.path(opt$outputDir, opt$prepReads_outputDir, 'reads.rds')), readID, adriftReadRandomID) %>% 
     dplyr::rename(randomLinkerSeq = adriftReadRandomID) %>% 
     data.table()

r <- r[readID %in% adriftReadAlignments$readID]

adriftReadAlignments <- left_join(adriftReadAlignments, r, by = 'readID')  

write(c(paste(now(), '   Preparing alignment data for fragment generation.')), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'log'), append = TRUE)

anchorReadAlignments <- subset(anchorReadAlignments, readID %in% adriftReadAlignments$readID)

dir.create(file.path(opt$outputDir, opt$buildFragments_outputDir))

anchorReadAlignments <- select(anchorReadAlignments, uniqueSample, sample, readID, tName, strand, tStart, tEnd, leaderSeq, refGenome, seqRunID, flags, vectorFastaFile)
adriftReadAlignments <- select(adriftReadAlignments, sample, readID, tName, strand, tStart, tEnd, randomLinkerSeq)

names(anchorReadAlignments) <- paste0(names(anchorReadAlignments), '.anchorReads')
names(adriftReadAlignments) <- paste0(names(adriftReadAlignments), '.adriftReads')

# Breaking reads into id groups prevents joins from exceeding internal row limits
# before they can be filtered for correct pairings.

ids <- unique(anchorReadAlignments$readID.anchorReads)
id_groups <- split(ids, dplyr::ntile(1:length(ids), ceiling(length(ids)/opt$buildFragments_idGroup_size)))

write(c(paste(now(), '   Building initial fragments.')), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'log'), append = TRUE)

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

write(paste0(now(), '    Building rational fragments from alignment data.'), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'log'), append = TRUE)

if(length(z) > 0 & opt$buildFragments_salvageReadsBeyondMaxNumAlignments){
  write(c(paste0(now(), '    ', length(z), ' reads pairs have more than ', opt$buildFragments_maxReadAlignments, ' alignments for both mates.')), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'log'), append = TRUE)
  write(c(paste0(now(), '    For each reach read, ', opt$buildFragments_maxReadAlignments, ' alignments will be randomly selected.')), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'log'), append = TRUE)
  write(c(paste0(now(), '    Adrift alignments that have the potential to form rational fragments with the selected anchor reads will be selected as well.')), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'log'), append = TRUE)
  
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
          
          if(nrow(y.pos) > 0 | nrow(y.neg) > 0){
            y2 <- bind_rows(y.pos, y.neg) %>% select(-start, -end)
          } else {
            y2 <- tibble()
          }
          
          list(x, y2)
  })
  
  a3 <- bind_rows(lapply(o, '[[', 1))
  b3 <- bind_rows(lapply(o, '[[', 2))
  
  anchorReadAlignments <- bind_rows(a1, a3)
  adriftReadAlignments <- bind_rows(b1, b3)
} else {
  anchorReadAlignments <- subset(anchorReadAlignments, ! readID.anchorReads %in% z)
  adriftReadAlignments <- subset(adriftReadAlignments, ! readID.adriftReads %in% z)
}


o <- lapply(id_groups, function(id_group){
       list(anchorReadAlignments[readID.anchorReads %in% id_group],
            adriftReadAlignments[readID.adriftReads %in% id_group])
     })

rm(anchorReadAlignments, adriftReadAlignments)
gc()

counter <- 1
total <- length(o)

# Running with parLapply() results in high memory usage.

frags <- bind_rows(lapply(o, function(z){
  a <- z[[1]]
  b <- z[[2]]
  
  if(nrow(a) == 0 | nrow(b) == 0) return(data.frame())
  
  write(paste0(now(), '    ', counter, '/', total, ': ', nrow(a), ' anchorRead alignments, ', nrow(b), ' adriftRead_alignments.'), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'log'), append = TRUE)
  counter <<- counter + 1
  
  # browser()
  # Join adrift reads alignments to anchor read alignments to create potential read pairs.
  # frags <- left_join(a, b, by = c('readID.anchorReads' = 'readID.adriftReads'), relationship = 'many-to-many') %>% tidyr::drop_na()
  frags <-  na.omit(a[b, on=c(readID.anchorReads = "readID.adriftReads"), allow.cartesian=TRUE])
  
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
                select(uniqueSample, readID, chromosome, strand, fragStart, fragEnd, leaderSeq.anchorReads, randomLinkerSeq.adriftReads, 
                       refGenome.anchorReads, vectorFastaFile.anchorReads, seqRunID.anchorReads, flags.anchorReads)
   r
}))

# Odd blat calls can lead to duplicate alignment entries, the blat parser is likely leaving out additional information about the alignments.
frags <- distinct(frags)

if(nrow(frags) == 0){
  write(c(paste(now(), '   Error - no fragments were identified.')), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'log'), append = TRUE)
  if(opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()
  q(save = 'no', status = 1, runLast = FALSE) 
}


write(c(paste(now(), '   Fragment generation complete.')), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'log'), append = TRUE)

# Add duplicate read count column from demultiplex module.
r <- readRDS(file.path(opt$outputDir, opt$prepReads_outputDir, 'reads.rds'))
frags <- left_join(frags, select(r, readID, nDuplicateReads))

rm(o, r, id_groups)
gc()


frags <- dplyr::rename(frags, leaderSeq = leaderSeq.anchorReads, randomLinkerSeq = randomLinkerSeq.adriftReads, refGenome = refGenome.anchorReads, 
                              vectorFastaFile = vectorFastaFile.anchorReads, seqRunID = seqRunID.anchorReads, flags = flags.anchorReads)

if(any(! incomingSamples %in% frags$uniqueSample) & opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()

saveRDS(frags, file.path(opt$outputDir, opt$buildFragments_outputDir, 'fragments.rds'), compress = opt$compressDataFiles)
write(date(), file.path(opt$outputDir, opt$buildFragments_outputDir, 'fragments.done'))


if('databaseGroup' %in% names(opt)){
  library(RMariaDB)
  
  write(c(paste(now(), '   Writing fragment data to the database.')), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'log'), append = TRUE)
  
  conn <- tryCatch({
    dbConnect(RMariaDB::MariaDB(), group = opt$databaseGroup)
  },
  error=function(cond) {
    write(c(paste(now(), '   Error - could not connect to the database.')), file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'log'), append = TRUE)
    if(opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()
    q(save = 'no', status = 1, runLast = FALSE) 
  })
  
  invisible(lapply(split(frags, frags$uniqueSample), function(x){
    x <- tidyr::separate(x, uniqueSample, c('trial', 'subject', 'sample', 'replicate'), sep = '~')
  
    dbExecute(conn, paste0("delete from fragments where trial='", x$trial[1], "' and subject='", x$subject[1],
                          "' and sample='", x$sample[1], "' and replicate='", x$replicate[1], "' and refGenome='", x$refGenome[1], "'"))
    
    o <- unlist(lapply(1:nrow(x), function(a){
           r <- x[a,]
           
           dbExecute(conn,
                     "insert into fragments values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                     params = list(r$trial, r$subject, r$sample, r$replicate, r$refGenome,
                                 r$vectorFastaFile, r$flags, r$seqRunID, r$nDuplicateReads, r$readID, 
                                 r$chromosome, r$strand, r$fragStart, r$fragEnd, r$leaderSeq, r$randomLinkerSeq))
    }))
    
    if(any(o == 0)){
      write(c(paste(now(), '   Error - could not upload all fragment records to the database.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)
      q(save = 'no', status = 1, runLast = FALSE) 
    }
  }))
  
  dbDisconnect(conn)
}

q(save = 'no', status = 0, runLast = FALSE) 
