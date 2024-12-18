#!/usr/bin/Rscript

# AAVengeR/filterVectorAgainstSeqs.R
# John K. Everett, Ph.D.
#
# This script reads in repLeaderSeq strings from buildSites or subsequent modules
# and aligns the against their respective vector sequences to build recombinations 
# maps which are expressed using a custom shorthand.

suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(data.table))

# Read in the configuration file and perform basic sanity checks.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) stop('Expected at least one command line argument')
source(file.path(yaml::read_yaml(args[1])$softwareDir, 'lib.R'))
opt <- startModule(args)


createOuputDir()
if(! dir.exists(file.path(opt$outputDir, opt$prepReads_outputDir))) dir.create(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir))
if(! dir.exists(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs'))) dir.create(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'dbs'))
if(! dir.exists(file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'))) dir.create(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'tmp'))
invisible(file.remove(list.files(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'tmp'), full.names = TRUE)))

# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)
write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)

quitOnErorr <- function(msg){
  updateLog(msg)
  message(msg)
  message(paste0('See log for more details: ', opt$defaultLogFile))
  q(save = 'no', status = 1, runLast = FALSE) 
}

updateLog('Starting mapSiteLeaderSequences.')

if(! file.exists(file.path(opt$outputDir, opt$mapSiteLeaderSequences_inputFile))) quitOnErorr("Error -- the sites input file could not be found.")

sites <- readRDS(file.path(opt$outputDir, opt$mapSiteLeaderSequences_inputFile))

if(! opt$mapSiteLeaderSequences_addAfter %in% names(sites)) quitOnErorr(paste0('Error - ', opt$mapSiteLeaderSequences_addAfter, ' is not a column in your input data frame.'))

cluster <- makeCluster(opt$mapSiteLeaderSequences_CPUs)
clusterExport(cluster, c('opt', 'blast2rearangements', 'buildRearrangementModel'))

# Loop over chunks of sites that used the same vector.
m <- rbindlist(lapply(split(sites, sites$vector), function(x){
 
  updateLog(paste0('Processing maps for vector: ', x$vector[1], '.'))
  
  # Clean up tmp database directory.
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'dbs'), full.names = TRUE)))
  
  # Make a blast database for the vector sequence.
  system(paste0('makeblastdb -in ', file.path(opt$softwareDir, 'data', 'vectors', x$vector[1]), 
                ' -dbtype nucl -out ', file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
  
  waitForFile(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'dbs', 'd.nin'))
  
  # Create a data frame of repLeader sequences to map.
  s <- data.frame(seq = unique(x$repLeaderSeq))
  s$id <- paste0('s', 1:nrow(s))
  
  # Bin the repLeader sequences by length using bins incremented by 10 nt.
  s$cut <- cut(nchar(s$seq), c(-Inf, seq(0, max(nchar(s$seq)), by = 10), Inf), labels = FALSE)
  
  # Create a splitting vector that breaks across bins. 
  s <- group_by(s, cut) %>% 
       mutate(n = ntile(1:n(), opt$mapSiteLeaderSequences_CPUs)) %>% 
       ungroup()
  
  reads <- DNAStringSet(s$seq)
  names(reads) <- s$id
  
  # Chunk the repLeaderSeqs and process each chunk.
  # Here qName stands for a particular repLeaderSeq while in other parts 
  # of the software it typically denotes a single read.
  
  b <- bind_rows(lapply(split(reads, s$n), function(a){
         f <- tmpFile()
 
         # Write the chunk out to a tmp file.
         writeXStringSet(a,  file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'tmp', paste0(f, '.fasta')))
    
         # Align the chunk to the vector sequence.
         system(paste0('blastn -dust no -soft_masking false -word_size 5 -evalue 50 -outfmt 6 -query ',
                  file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'tmp', paste0(f, '.fasta')), ' -db ',
                  file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'dbs', 'd'),
                  ' -num_threads 1 -out ', file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'tmp', paste0(f, '.blast'))),
                ignore.stdout = TRUE, ignore.stderr = TRUE)
    
         waitForFile(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'tmp', paste0(f, '.blast')))
    
         # Catch instances where no alignments were returned.
         if(file.info(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'tmp', paste0(f, '.blast')))$size == 0) return(tibble())
    
         # Parse blastn result.
         b <- read.table(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'tmp', paste0(f, '.blast')), sep = '\t', header = FALSE)
         names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
         b$alignmentLength <- b$qend - b$qstart + 1
         b$strand <- ifelse(b$sstart > b$send, '-', '+')
          
         # Limit repLeaderSeqs in this chunk to those with one or more alignments.
         o <- a[names(a) %in% b$qname]
          
         # Add query lengths to the blastn alignments.
         d <- tibble(qname = names(o), qlen = width(o))
         b <- left_join(b, d, by = 'qname')
    
         # Limit alignments based on parameters in the configuration file.
         dplyr::filter(b, pident >= opt$mapSiteLeaderSequences_minAlignmentPercentID, alignmentLength >= opt$mapSiteLeaderSequences_minAlignmentLength)
       }))
  
  # Create a grouping index for read IDs.
  b$i <- group_by(b, qname) %>% group_indices()
  
  # Create a CPU splitting vector that would not separate the read ID indices. 
  o <- tibble(i2 = 1:n_distinct(b$i))
  o$n <- ntile(1:nrow(o), opt$prepReads_CPUs)
  
  # Bind the CPU splitting vector.
  b <- left_join(b, o, by = c('i' = 'i2'))
  
  # Pass alignment chunks to rearrangement function.
  
  #r <- rbindlist(lapply(split(b, b$n), blast2rearangements, opt$mapSiteLeaderSequences_maxMissingTailNTs))
  r <- rbindlist(parallel::parLapply(cluster, split(b, b$n), blast2rearangements, opt$mapSiteLeaderSequences_maxMissingTailNTs, opt$mapSiteLeaderSequences_minLocalAlignmentLength))
  
  z <- left_join(tibble(qname = names(reads), leaderSeq = as.character(reads)), dplyr::rename(r, repLeaderSeqMap = rearrangement), by = 'qname')
  z$vector <- x$vector[1]
  z
}))

stopCluster(cluster)

invisible(file.remove(list.files(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'tmp'), full.names = TRUE)))

sites <- bind_rows(lapply(split(sites, sites$vector), function(x){
      m2 <- subset(m, vector == x$vector[1])
      m2 <- m2[! duplicated(m2$leaderSeq),]
      left_join(x, select(m2, leaderSeq, repLeaderSeqMap), by = c('repLeaderSeq' = 'leaderSeq'))
}))

sites <- dplyr::relocate(sites, repLeaderSeqMap, .after = opt$mapSiteLeaderSequences_addAfter)

if(tolower(opt$databaseConfigGroup) != 'none'){
  suppressPackageStartupMessages(library(RMariaDB))
  uploadSitesToDB(sites)
}

saveRDS(sites, file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'sites.rds'))
openxlsx::write.xlsx(sites, file = file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'sites.xlsx'))
readr::write_tsv(sites, file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'sites.tsv.gz'))

updateLog('mapSiteLeaderSequences completed.')

q(save = 'no', status = 0, runLast = FALSE) 
