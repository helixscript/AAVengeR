# John K. Everett, PhD
# AAVengeR/anchorReadRearrangements.R

library(ShortRead)
library(dplyr)
library(parallel)
library(lubridate)
library(data.table)
library(stringr)

# Read in config file and source library.
configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))


# Create module directory structure.
dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir))
dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs'))
dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'))

write(c(paste(now(), '   Reading in demultiplexed reads.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
reads <- readRDS(file.path(opt$outputDir, opt$anchorReadRearrangements_readsTable))


# Create virtual cluster.
cluster <- makeCluster(opt$anchorReadRearrangements_CPUs)
clusterExport(cluster, 'opt')

# Trim anchor read over-reading with cutadapt using the RC of the common linker in adrift reads.
# The trim sequence is created by demultiplex.R.
write(c(paste(now(), '   Trimming anchor read over-reading.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

reads <- data.table::rbindlist(parLapply(cluster, split(reads, dplyr::ntile(1:nrow(reads), opt$anchorReadRearrangements_CPUs)), function(x){
  source(file.path(opt$softwareDir, 'lib.R'))
  library(dplyr)
  library(data.table)
  library(Biostrings)
  
  data.table::rbindlist(lapply(split(x, x$adriftReadTrimSeq), function(y){
    f <- file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', tmpFile())
    o <- DNAStringSet(y$anchorReadSeq)
    names(o) <- y$readID
    Biostrings::writeXStringSet(o, f)
    
    # Use cut adapt to trim anchor read over-reading.
    system(paste0('cutadapt -e 0.15 -a ', y$adriftReadTrimSeq[1], ' ',
                  f, ' > ', paste0(f, '.cutAdapt')), ignore.stderr = TRUE)
    
    t <- readDNAStringSet(paste0(f, '.cutAdapt'))
    invisible(file.remove(f, paste0(f, '.cutAdapt')))
    
    t <- t[width(t) >= opt$anchorReadRearrangements_minAnchorReadLength]
    if(length(t) == 0) return(data.table())
    
    # Add back trimmed anchor read sequences.
    o <- subset(y, readID %in% names(t))
    trimmed <- data.table(readID = names(t), anchorReadSeq2 = as.character(t))
    o <- left_join(o, trimmed, by = 'readID')
    
    # Trim off adrift read adapters.
    o$adriftReadSeq2 <- substr(o$adriftReadSeq, o$adriftLinkerSeqEnd+1, nchar(o$adriftReadSeq))
    o <- o[nchar(o$adriftReadSeq2) >= opt$anchorReadRearrangements_minAdriftReadLength,]
    
    data.table(dplyr::select(o, -adriftReadTrimSeq, -adriftLinkerSeqEnd))
  }))
}))

if(nrow(reads) == 0){
  write(c(paste(now(), "Error - no reads remaining after trimming.")), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}


# Report trimming stats.
trimReport <- mutate(reads, sample = sub('~\\d+$', '', uniqueSample)) %>% 
              group_by(sample) %>% 
              summarise(reads = n_distinct(readID), percentTrimmed = sprintf("%.2f%%", (sum(anchorReadSeq != anchorReadSeq2)/n())*100)) %>% 
              ungroup()

trimReport$sample <- paste0('                       ', trimReport$sample)
write(paste0(now(),         '    Anchor reads trimmed with over-read adapter (sample, reads, percent trimmed):'), file = file.path(opt$outputDir, 'log'), append = TRUE)
readr::write_tsv(trimReport, file = file.path(opt$outputDir, 'log'), append = TRUE, col_names = FALSE)


# Switch original sequences for trimmed sequences.
reads$anchorReadSeq <- NULL
reads$adriftReadSeq <- NULL
reads <- dplyr::rename(reads, anchorReadSeq = anchorReadSeq2, adriftReadSeq = adriftReadSeq2)


# Create additional identifiers from the start of trimmed adrift reads.
# This sequence ids will help identify PCR linker crossovers. 
reads$adriftReadSeqID <- substr(reads$adriftReadSeq, 1, 15)
reads$adriftReadSeq <- NULL


# Identify and remove identical read pairs.
write(c(paste(now(), '   Identifying duplicate read pairs.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
reads$i <- group_by(reads, uniqueSample, anchorReadSeq, adriftReadSeqID) %>% group_indices() 
reads <- group_by(reads, i) %>% mutate(n = n()) %>% ungroup()

a <- subset(reads, n == 1) # Non-duplicated read pairs.
b <- subset(reads, n > 1)  # Duplicated read pairs.

write(c(paste(now(), '   Removing duplicate read pairs.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

# Create a table of duplicate read pairs where one is chosen (id) to move forward and the others are logged (id2).
c <- group_by(b, i) %>%
     summarise(id = readID[1], n = n() - 1, id2 = list(readID[2:n()])) %>%
     ungroup() %>% select(-i) %>% tidyr::unnest(id2)

# Exclude duplicate read pairs.
reads <- data.table(bind_rows(a, subset(b, ! readID %in% c$id2)) %>% dplyr::select(-i, -n))

# Clean up.
stopCluster(cluster)
invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), full.names = TRUE)))
rm(a, b, c)
gc()


# Correct minor UMI sequencing error by comparing to abundant UMIs within samples.
cluster <- makeCluster(opt$anchorReadRearrangements_CPUs)
clusterExport(cluster, c('opt', 'conformMinorSeqDiffs', 'tmpFile', 'waitForFile'))

reads <- unpackUniqueSampleID(reads)
reads <- select(reads, readID, adriftReadRandomID, adriftReadSeqID, anchorReadSeq, trial, subject, sample, replicate, vectorFastaFile)

reads <- rbindlist(parLapply(cluster, split(reads, reads$sample), function(x){
           library(dplyr)
           library(stringdist)
  
           x$adriftReadRandomID <- conformMinorSeqDiffs(x$adriftReadRandomID, 
                                                        abundSeqMinCount = opt$anchorReadRearrangements_abundSeqMinCount, 
                                                        nThreads = 3)
           x
         }))

# Remove reads where there were conflicts in the UMI correction (returned as poly N).
reads <- reads[! grepl('NNN', reads$adriftReadRandomID),]

            
# Determine which UMIs have multiple reads then select a representative read for each UMI.
d <- unique(reads$adriftReadRandomID[duplicated(reads$adriftReadRandomID)])
a <- subset(reads, ! adriftReadRandomID %in% d)
b <- subset(reads, adriftReadRandomID %in% d)

# Create grouping indices for each UMI and then use those indicies to create 
# a second index for splitting the data across CPUs.
b$i <- group_by(b, adriftReadRandomID) %>% group_indices() 
o <- tibble(i2 = 1:n_distinct(b$i))
o$n <- ntile(1:nrow(o), opt$anchorReadRearrangements_CPUs)
b <- left_join(b, o, by = c('i' = 'i2'))


# For each UMI, the start of R1 should be the same unless their was PCR recombination.
# Here, for each UMI, we look for the major start sequence for R1, select reads with 
# that start sequecne, and then return the concensus sequence. 
concensusReadSeq <- function(x){
  if(length(x) > 1000) x <- dplyr::sample_n(x, 1000)
  
  o <- sort(table(x$adriftReadSeqID), decreasing = TRUE)
  
  if(o[1]/nrow(x) > 0.5){
    x$anchorReadSeq <- as.character(Biostrings::consensusString(Biostrings::DNAStringSet(subset(x, adriftReadSeqID == names(o[1]))$anchorReadSeq), 
                                                                threshold = 0.75, 
                                                                ambiguityMap = 'N'))
  } else {
    x$anchorReadSeq <- NA
  }
  
  x[1,]
}

clusterExport(cluster, 'concensusReadSeq')

b <- bind_rows(parLapply(cluster, split(b, b$n), function(x){
       library(dplyr)
       library(Biostrings)
  
       set.seed(1)
       bind_rows(lapply(split(x, x$adriftReadRandomID), concensusReadSeq))
     }))

# Remove UMIs for which PCR crossovers could not be resolved.
b <- b[! is.na(b$anchorReadSeq)]

# Recombine read data.
n1 <- n_distinct(reads$adriftReadRandomID)
reads <- bind_rows(a, b)
n2 <- n_distinct(reads$adriftReadRandomID)


write(c(paste(now(), '   ', sprintf("%.2f%%", (1 - n2/n1)*100), 'of UMIs removed because potential PCR cross-overs could not be resolved.')), file = file.path(opt$outputDir, 'log'), append = TRUE)


blastWorker <- function(y, wordSize = 6, evalue = 10){
  library(Biostrings)
  library(data.table)
  library(dplyr)
  
  s <- DNAStringSet(y$anchorReadSeq)
  names(s) <- y$readID
  
  f <- tmpFile()
  writeXStringSet(s,  file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.fasta')))
  
  system(paste0('blastn -dust no -soft_masking false -word_size ', wordSize, ' -evalue ', evalue,' -outfmt 6 -query ',
                file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.fasta')), ' -db ',
                file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs', 'd'),
                ' -num_threads 1 -out ', file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.blast'))),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  waitForFile(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.blast')))
  
  b <- tibble()
  
  if(file.info(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.blast')))$size > 0){
    b <- read.table(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.blast')), sep = '\t', header = FALSE)
    names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
    
    o <- s[names(s) %in% b$qname]
    d <- tibble(qname = names(o), qlen = width(o))
    b <- left_join(b, d, by = 'qname')
    
    b$alignmentLength <- b$qend - b$qstart + 1
    b <- subset(b, pident >= opt$anchorReadRearrangements_seqsMinPercentID & 
                  alignmentLength >= opt$anchorReadRearrangements_seqsMinAlignmentLength &
                  gapopen <= opt$anchorReadRearrangements_seqsMaxGaps)
  }
  
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), pattern = f, full.names = TRUE)))
  
  b
}

# Restarting the cluster because random seeds were previously set within the cores
# which will result in all tmp files having the same names.

stopCluster(cluster)
cluster <- makeCluster(opt$anchorReadRearrangements_CPUs)
clusterExport(cluster, c('opt', 'blastWorker', 'tmpFile', 'waitForFile', 'buildRearrangementModel'))
        

write(c(paste(now(), '    Aligning reads to vector files.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
                 
vectorHits <- rbindlist(lapply(split(reads, reads$vectorFastaFile), function(x){
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$vectorFilter_outputDir, 'dbs'), full.names = TRUE)))

  system(paste0('makeblastdb -in ', file.path(opt$softwareDir, 'data', 'vectors', x$vectorFastaFile[1]), ' -dbtype nucl -out ', file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
  waitForFile(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs', 'd.nin'))
  
  o <- split(x, ntile(1:nrow(x), opt$anchorReadRearrangements_CPUs))
  
  rbindlist(parLapply(cluster, o, blastWorker))
}))

write(c(paste(now(), '    Vector alignments completed.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

o <- group_by(vectorHits, qname) %>% summarise(n = any(qstart >= 25)) %>% ungroup()

write(c(paste(now(), '   ', sprintf("%.2f%%", (sum(o$n)/ n_distinct(o$qname))*100), 'of reads that align to the vector files have a start position >= 25')), file = file.path(opt$outputDir, 'log'), append = TRUE)


# buildModel <- function(b2){
#   r <- vector()
#   counter <- 0
#   while(nrow(b2) != 0){
#     message(counter)
#     counter <- counter + 1
#     b2 <- arrange(b2, qstart, evalue)
#     x <- b2[1,]
#     if((x$qend - x$qstart + 1) < opt$anchorReadRearrangements_seqsMinAlignmentLength) next
#     r <- paste0(r, ';', x$qstart, '..', x$qend, '[', x$sstart2, x$strand, x$send2, ']')
#     b2 <- b2[abs(b2$qstart - x$qend) < 5,] # Remove alignments that were considered during this loop.
#     if(counter == 1000) break
#   } 
#   
#   sub('^;', '', r)
# }
# 
# blast2rearangements_worker <- function(b){
#   library(dplyr)
#   library(IRanges)
#   library(data.table)
#   
#   rbindlist(lapply(split(b, b$qname), function(b2){
# 
#     # Sort BLAST results by query start position and evalue (low to high).
#     b2 <- arrange(b2, qstart, evalue)
#     b2$strand <- ifelse(b2$send < b2$sstart, '-', '+')
#     
#     # Alignment to the negative strand will result in the subject end to come before the start.
#     # Switch it back so that they are sequential.
#     b2$sstart2 <- ifelse(b2$sstart > b2$send, b2$send, b2$sstart)
#     b2$send2   <- ifelse(b2$sstart > b2$send, b2$sstart, b2$send)
#     
#     data.table(qname = b2$qname[1], rearrangement = buildModel(b2))
#   }))

# Create a splitting CPU splitting vector that will not split up read ids.
vectorHits$i <- group_by(vectorHits, qname) %>% group_indices() 
o <- tibble(i2 = 1:n_distinct(vectorHits$i))
o$n <- ntile(1:nrow(o), opt$anchorReadRearrangements_CPUs)
vectorHits <- left_join(vectorHits, o, by = c('i' = 'i2'))


write(c(paste(now(), '    Building rearrangement models for each read.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

# Create alignment rearrangement models for each read.
r <- rbindlist(parallel::parLapply(cluster, split(vectorHits, vectorHits$n), blast2rearangements))
# r <- rbindlist(lapply(split(vectorHits, vectorHits$n), blast2rearangements))

# Add annotation for missing alignments based on read length.

write(c(paste(now(), '    Adding default annotations for reads without alignments.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

reads$width <- nchar(reads$anchorReadSeq)
r <- left_join(r, select(reads, readID, width), by = c('qname' = 'readID'))

r$n <- ntile(1:nrow(r), opt$anchorReadRearrangements_CPUs)

r2 <- bind_rows(parLapply(cluster, split(r, r$n), function(y){
        library(dplyr)
        library(stringr)
  
        bind_rows(lapply(split(y, 1:nrow(y)), function(x){
          if(is.na(x$rearrangement)) return(NA)
  
          o <- unlist(strsplit(as.character(x$rearrangement), ';'))
          lastRangeEnd <- as.integer(sub('\\.\\.', '', str_extract(o[length(o)], '..(\\d+)')))
  
          # Add unknown segment to end of reads if last range is shorter than read length.
          if((x$width - lastRangeEnd) > opt$anchorReadRearrangements_minAllowableGap){
            x$rearrangement <- paste0(x$rearrangement, ';', lastRangeEnd+1, '..', x$width, '[x]')
          }
  
          o <- unlist(strsplit(as.character(x$rearrangement), ';'))
  
          if(length(o) > 1){
           
            invisible(lapply(1:(length(o)-1), function(n){
              lastRangeEnd <- as.integer(sub('\\.\\.', '', str_extract(o[[n]], '..(\\d+)')))
              nextFirstRangeEnd <- as.integer(sub('\\.\\.', '', str_extract(o[[n+1]], '(\\d+)..')))
      
              if((nextFirstRangeEnd - lastRangeEnd) > opt$anchorReadRearrangements_minAllowableGap ){
               o[[n]] <<- paste0(o[[n]], ';', lastRangeEnd+1, '..', nextFirstRangeEnd-1, '[x];')
              }
            }))
          }
  
          x$rearrangement <- gsub(';;', ';', paste0(o, collapse = ';'))
          x
        }))
}))

# Remove last unknown segments since they may be genomic sequences from integrated 
# vectors or poor base calls at the ends of reads.
if(opt$anchorReadRearrangements_removeTailingUnknownSegments) r2$rearrangement <- sub(';\\d+\\.\\.\\d+\\[x\\]$', '', r2$rearrangement)

# Add read metadata.
write(c(paste(now(), '    Adding metadata to rearrangement models.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

r2 <- left_join(r2, select(reads, trial, subject, sample, readID, adriftReadRandomID), by = c('qname' = 'readID'))

r3 <- bind_rows(lapply(split(r2, paste(r2$trial, r2$subject, r2$sample)), function(x){
        mutate(x, UMIs = n_distinct(adriftReadRandomID), percentUMIsRearranged = sprintf("%.2f%%", (sum(grepl(';', x$rearrangement)) / UMIs)*100)) %>%
        select(trial, subject, sample, UMIs, percentUMIsRearranged) %>% 
        dplyr::slice(1)
       }))

saveRDS(r2, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'readRearrangements.rds'))
openxlsx::write.xlsx(r2, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'readRearrangements.xlsx'))
readr::write_tsv(r2, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'readRearrangements.tsv'))

saveRDS(r3, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'result.rds'))
openxlsx::write.xlsx(r3, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'result.xlsx'))
readr::write_tsv(r3, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'result.tsv'))
