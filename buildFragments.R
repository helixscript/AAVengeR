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
  r <- left_join(r, dplyr::select(adriftReadAlignments, qName, sample), by = c('readID' = 'qName'))
  
  # Remove replicate identifiers from sample names to evaluate on the sample level.
  r$sample <- sub('~\\d+$', '', r$sample)
  
  # Correct random ids to the most abundant within samples.
  # Uncorrectable codes are returned as NNNN and removed.
  message('Correcting minor differences in random ids.')
  # r <- bind_rows(lapply(split(r, r$sample), function(x){
  r <- bind_rows(parLapply(cluster, split(r, r$sample), function(x){
         source(file.path(opt$softwareDir, 'lib.R'))
         library(dplyr)
         x$randomLinkerSeq <- conformMinorSeqDiffs(x$randomLinkerSeq)
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

readID2reference <- tibble(readID = anchorReadAlignments$qName, refGenome = anchorReadAlignments$refGenome.id)

anchorReadAlignments <- subset(anchorReadAlignments, qName %in% adriftReadAlignments$qName)

dir.create(file.path(opt$outputDir, opt$buildFragments_outputDir))

anchorReadAlignments <- select(anchorReadAlignments, sample, qName, tName, strand, tStart, tEnd, leaderSeq)
adriftReadAlignments <- select(adriftReadAlignments, sample, qName, tName, strand, tStart, tEnd, randomLinkerSeq)

names(anchorReadAlignments) <- paste0(names(anchorReadAlignments), '.anchorReads')
names(adriftReadAlignments) <- paste0(names(adriftReadAlignments), '.adriftReads')



# Limit the number of alignments per read.

k <- table(anchorReadAlignments$qName.anchorReads)
p <- sprintf("%.2f%%", (sum(k > opt$buildFragments_maxReadAlignments)/n_distinct(anchorReadAlignments$qName.anchorReads))*100)
write(c(paste(now(), paste0('   ', p, '  anchor reads have more than ', opt$buildFragments_maxReadAlignments, ' alignments.'))), file = file.path(opt$outputDir, 'log'), append = TRUE)
k <- k[k > opt$buildFragments_maxReadAlignments]

if(length(k) > 0){
  a <- anchorReadAlignments[anchorReadAlignments$qName.anchorReads %in% names(k),]
  b <- anchorReadAlignments[! anchorReadAlignments$qName.anchorReads %in% names(k),]
  
  a2 <- bind_rows(lapply(split(a, a$qName.anchorReads), function(x){
          set.seed(1)
          dplyr::slice_sample(x, n = opt$buildFragments_maxReadAlignments)
         }))
  
  anchorReadAlignments <- bind_rows(a2, b)
}


k <- table(adriftReadAlignments$qName.adriftReads)
p <- sprintf("%.2f%%", (sum(k > opt$buildFragments_maxReadAlignments)/n_distinct(adriftReadAlignments$qName.adriftReads))*100)
write(c(paste(now(), paste0('   ', p, '  adrift reads have more than ', opt$buildFragments_maxReadAlignments, ' alignments.'))), file = file.path(opt$outputDir, 'log'), append = TRUE)
k <- k[k > opt$buildFragments_maxReadAlignments]

if(length(k) > 0){
  a <- adriftReadAlignments[adriftReadAlignments$qName.adriftReads %in% names(k),]
  b <- adriftReadAlignments[! adriftReadAlignments$qName.adriftReads %in% names(k),]
  
  a2 <- bind_rows(lapply(split(a, a$qName.adriftReads), function(x){
    set.seed(1)
    dplyr::slice_sample(x, n = opt$buildFragments_maxReadAlignments)
  }))
  
  adriftReadAlignments <- bind_rows(a2, b)
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
