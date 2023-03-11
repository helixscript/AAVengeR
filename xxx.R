# John K. Everett, PhD
# AAVengeR/prepReads.R
#
# This module prepares reads for alignment to a reference genome by removing duplicate
# read pairs, trimming adapter sequences, removing reads with high sequence homology to 
# the vector and selecting anchor reads with the expected leader sequences. 

library(ShortRead)
library(dplyr)
library(parallel)
library(lubridate)
options(stringsAsFactors = FALSE)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)

invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$prepReads_outputDir))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs'))

write(c(paste(now(), '   Reading in demultiplexed reads.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
reads <- readRDS(file.path(opt$outputDir, opt$prepReads_readsTable))

cluster <- makeCluster(opt$prepReads_CPUs)
clusterExport(cluster, 'opt')

# Trim anchor read over-reading with cutadapt using the RC of the common linker in adrift reads.
# The trim sequence is created by demultiplex.R.

write(c(paste(now(), '   Trimming anchor read over-reading.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

reads <- data.table::rbindlist(parLapply(cluster, split(reads, dplyr::ntile(1:nrow(reads), opt$prepReads_CPUs)), function(x){
  source(file.path(opt$softwareDir, 'lib.R'))
  library(dplyr)
  library(data.table)
  library(Biostrings)
  
  data.table::rbindlist(lapply(split(x, x$adriftReadTrimSeq), function(y){
    f <- file.path(opt$outputDir, 'tmp',  tmpFile())
    o <- DNAStringSet(y$anchorReadSeq)
    names(o) <- y$readID
    Biostrings::writeXStringSet(o, f)
    
    # Use cut adapt to trim anchor read over-reading.
    system(paste0('cutadapt -e 0.15 -a ', y$adriftReadTrimSeq[1], ' ',
                  f, ' > ', paste0(f, '.cutAdapt')), ignore.stderr = TRUE)
    
    t <- readDNAStringSet(paste0(f, '.cutAdapt'))
    invisible(file.remove(f, paste0(f, '.cutAdapt')))
    
    t <- t[width(t) >= opt$prepReads_minAnchorReadLength]
    if(length(t) == 0) return(data.table())
    
    # Add back trimmed anchor read sequences.
    o <- subset(y, readID %in% names(t))
    trimmed <- data.table(readID = names(t), anchorReadSeq2 = as.character(t))
    o <- left_join(o, trimmed, by = 'readID')
    
    # Trim off adrift read adapters.
    o$adriftReadSeq2 <- substr(o$adriftReadSeq, o$adriftLinkerSeqEnd+1, nchar(o$adriftReadSeq))
    o <- o[nchar(o$adriftReadSeq2) >= opt$prepReads_minAdriftReadLength,]
    
    data.table(dplyr::select(o, -adriftReadTrimSeq, -adriftLinkerSeqEnd))
  }))
}))

stopCluster(cluster)

if(nrow(reads) == 0){
  write(c(paste(now(), "Error - no reads remaining after trimming.")), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}

trimReport <- mutate(reads, sample = sub('~\\d+$', '', uniqueSample)) %>% 
  group_by(sample) %>% 
  summarise(reads = n_distinct(readID), percentTrimmed = sprintf("%.2f%%", (sum(anchorReadSeq != anchorReadSeq2)/n())*100)) %>% 
  ungroup()

trimReport$sample <- paste0('                       ', trimReport$sample)
write(paste0(now(),         '    Anchor reads trimmed with over-read adapter (sample, reads, percent trimmed):'), file = file.path(opt$outputDir, 'log'), append = TRUE)
readr::write_tsv(trimReport, file = file.path(opt$outputDir, 'log'), append = TRUE, col_names = FALSE)

reads$anchorReadSeq <- NULL
reads$adriftReadSeq <- NULL
reads <- dplyr::rename(reads, anchorReadSeq = anchorReadSeq2, adriftReadSeq = adriftReadSeq2)


# Identify and remove identical read pairs.
write(c(paste(now(), '   Identifying duplicate read pairs.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
reads$i <- group_by(reads, uniqueSample, adriftReadRandomID, anchorReadSeq, adriftReadSeq) %>% group_indices() 
reads <- group_by(reads, i) %>% mutate(n = n()) %>% ungroup()

a <- subset(reads, n == 1) # Non-duplicated read pairs.
b <- subset(reads, n > 1)  # Duplicated read pairs.

write(c(paste(now(), '   Removing duplicate read pairs.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

# Create a table of duplicate read pairs where one is chosen (id) to move forward and the others are logged (id2).
c <- group_by(b, i) %>%
  summarise(id = readID[1], n = n() - 1, id2 = list(readID[2:n()])) %>%
  ungroup() %>% select(-i) %>% tidyr::unnest(id2)

saveRDS(c, file.path(opt$outputDir, opt$prepReads_outputDir, 'duplicateReads.rds'))


# Exclude duplicate read pairs.
reads <- data.table(bind_rows(a, subset(b, ! readID %in% c$id2)) %>% dplyr::select(-i, -n))
saveRDS(reads, file.path(opt$outputDir, opt$prepReads_outputDir, 'uniqueReadPairs.rds'))

# Clean up.
invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))
rm(a, b, c)
gc()