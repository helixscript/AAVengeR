# AAVengeR/alignReads.R
# John K. Everett, Ph.D.
# 
# This scripts accepts prepared reads from the prepReads module and aligns
# them to reference genomes. Alignments are filtered according to parameters
# found in the configuration file.

suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(RMariaDB))
options(stringsAsFactors = FALSE)

set.seed(1)

source(file.path(yaml::read_yaml(commandArgs(trailingOnly=TRUE)[1])$softwareDir, 'lib.R'))
opt <- loadConfig()
optionsSanityCheck()

createOuputDir()
if(! dir.exists(file.path(opt$outputDir, opt$alignReads_outputDir))) dir.create(file.path(opt$outputDir, opt$alignReads_outputDir))
if(! dir.exists(file.path(opt$outputDir, opt$alignReads_outputDir, 'tmp'))) dir.create(file.path(opt$outputDir, opt$alignReads_outputDir, 'tmp'))
if(! dir.exists(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1'))) dir.create(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1'))
if(! dir.exists(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2'))) dir.create(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2'))

# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, opt$alignReads_outputDir, 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)
write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)

quitOnErorr <- function(msg){
  if(opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()
  updateLog(msg)
  message(msg)
  message(paste0('See log for more details: ', opt$defaultLogFile))
  q(save = 'no', status = 1, runLast = FALSE) 
}

# Read in sequencing data.
updateLog('Reading in prepped reads.')

if(! file.exists(file.path(opt$outputDir, opt$alignReads_inputFile))) quitOnErorr('Error - the input file does not exist.')

reads <- readRDS(file.path(opt$outputDir, opt$alignReads_inputFile))

if(nrow(reads) == 0) quitOnErorr('Error - the input table does not contain any rows.')

if(previousSampleDatabaseCheck(tibble(uniqueSample = reads$uniqueSample, refGenome = reads$refGenome) %>% tidyr::separate(uniqueSample, c('trial', 'subject', 'sample', 'replicate'), sep = '~')  %>% distinct())) q(save = 'no', status = 1, runLast = FALSE) 

incomingSamples <- unique(reads$uniqueSample)

# Create a CPU cluster.
cluster <- makeCluster(opt$alignReads_CPUs)
clusterSetRNGStream(cluster, 1)
clusterExport(cluster, c('opt', 'tmpFile'))

alignReads <- function(r, refGenome, minPercentSeqID, maxQstart, dir){
  
  refGenomePath <-  file.path(opt$softwareDir, 'data', 'referenceGenomes', 'blat', paste0(refGenome, '.2bit'))
      
  # Create an OOC file if requested
  if(opt$alignReads_blatUseOocFile){
      system(paste0('blat ', refGenomePath, ' /dev/null /dev/null -repMatch=',
                    opt$alignReads_genomeAlignment_blatRepMatch, ' -makeOoc=',
                    file.path(opt$outputDir, opt$alignReads_outputDir, paste0(opt$alignReads_genomeAlignment_blatTileSize, '.ooc'))))
  }

  # Create a split vector that attempts to equally distributes different length sequences between CPUs.
  r$cut <- cut(nchar(r$seq), c(-Inf, seq(0, max(nchar(r$seq)), by = 10), Inf), labels = FALSE)
  r <- group_by(r, cut) %>% mutate(n = ntile(1:n(), opt$alignReads_CPUs)) %>% ungroup()
  invisible(parLapply(cluster, split(r, r$n), blat, refGenomePath, dir))
  
  b <- tibble()
  f <- list.files(dir, pattern = '*.psl', full.names = TRUE)
  
  if(length(f) > 0){
    b <- rbindlist(lapply(f, function(x){
           b <- data.table(parseBLAToutput(x))
           
           if(nrow(b) == 0) return(data.table())
         
           dplyr::filter(b, alignmentPercentID >= minPercentSeqID, tNumInsert <= 1, qNumInsert <= 1, 
                           tBaseInsert <= 2, qBaseInsert <= 2, qStart <= maxQstart) %>%
           dplyr::select(qName, strand, qSize, qStart, qEnd, tName, tStart, tEnd, alignmentPercentID)
        }))
  }

  # Files need to be removed otherwise the join below will keep joining files from previous genomes.
  invisible(file.remove(list.files(dir, full.names = TRUE)))

  b
}

# Align anchor reads.

waitForMemory(stepDesc = 'Align reads / anchor reads, replicate level jobs', minMem = opt$system_minMemThreshold, maxWaitSecs = opt$system_minMemWaitTime, sleepSecs = opt$system_minMemSleepTime)

updateLog(paste0('Aligning ', ppNum(n_distinct(reads$readID)), ' anchor reads to reference genomes.'))

anchorReadAlignments <- rbindlist(lapply(split(reads, reads$refGenome), function(x){
  s <- data.table(seq = unique(x$anchorReadSeq))
  s$id <- paste0('s', 1:nrow(s))
  
  b <- alignReads(s, x$refGenome[1], opt$alignReads_genomeAlignment_minPercentID, opt$alignReads_genomeAlignment_anchorRead_maxStartPos,  file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1'))
  
  if(nrow(b) > 0){
    # Use the unique sequences to map alignments back to sequences
    b <- left_join(b, distinct(dplyr::select(s, id, seq)), by = c('qName' = 'id'))
    
    # Select reads that have an alignment.
    x2 <- dplyr::select(x, uniqueSample, readID, refGenome, anchorReadSeq) %>% dplyr::filter(anchorReadSeq %in% b$seq)
    
    # Join alignments to all reads using aligned sequences as keys.
    b <- left_join(x2, b, by = c('anchorReadSeq' = 'seq'), relationship = 'many-to-many') %>% dplyr::select(-anchorReadSeq, -qName)
  }
  
  b
}))

if(nrow(anchorReadAlignments) == 0) quitOnErorr('Error - no anchor reads aligned to the reference genomes.')

updateLog(paste0(sprintf("%.2f%%", (n_distinct(anchorReadAlignments$readID)/n_distinct(reads$readID))*100), 
                 ' of prepped anchor reads aligned to the reference genome.'))

# Select anchor reads where the ends align to the genome.
i <- (anchorReadAlignments$qSize - anchorReadAlignments$qEnd) <= opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned

if(sum(i) == 0) quitOnErorr('Error - no anchor read alignments remain after alignReads_genomeAlignment_anchorReadEnd_maxUnaligned filter.')

alignedReadIDsBeforeFilter <- n_distinct(anchorReadAlignments$readID)
anchorReadAlignments <- anchorReadAlignments[i,]

updateLog(paste0(sprintf("%.2f%%", (1 - n_distinct(anchorReadAlignments$readID) / alignedReadIDsBeforeFilter)*100), 
                 ' of anchor reads removed because their alignments ended more than ',
                 opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned,
                 ' NTs from the end of reads.'))

# Subset reads to those with good anchor read alignments.
reads <- subset(reads, readID %in% anchorReadAlignments$readID)


# Align adrift reads.

waitForMemory(stepDesc = 'Align reads / adrift reads, replicate level jobs', minMem = opt$system_minMemThreshold, maxWaitSecs = opt$system_minMemWaitTime, sleepSecs = opt$system_minMemSleepTime)

updateLog(paste0('Aligning ', ppNum(n_distinct(reads$readID)), ' adrift reads to reference genomes.'))

adriftReadAlignments <- rbindlist(lapply(split(reads, reads$refGenome), function(x){
  s <- data.table(seq = unique(x$adriftReadSeq))
  s$id <- paste0('s', 1:nrow(s))
  
  b <- alignReads(s, x$refGenome[1], opt$alignReads_genomeAlignment_minPercentID, opt$alignReads_genomeAlignment_adriftRead_maxStartPos,  file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2'))
  
  if(nrow(b) > 0){
    # Use the unique sequences to map alignments back to sequences
    b <- left_join(b, distinct(dplyr::select(s, id, seq)), by = c('qName' = 'id'))
    
    # Select reads that have an alignment.
    x2 <- dplyr::select(x, uniqueSample, readID, refGenome, adriftReadSeq) %>% dplyr::filter(adriftReadSeq %in% b$seq)
    
    # Join alignments to all reads using aligned sequences as keys.
    b <- left_join(x2, b, by = c('adriftReadSeq' = 'seq'), relationship = 'many-to-many') %>% dplyr::select(-adriftReadSeq, -qName)
  }
  
  b
}))

if(nrow(adriftReadAlignments) == 0) quitOnErorr('Error - no adrift reads aligned to the reference genomes.')

updateLog(paste0(sprintf("%.2f%%", (n_distinct(adriftReadAlignments$readID)/n_distinct(reads$readID))*100), 
                 ' of prepped adrift reads aligned to the reference genome.'))

# Select adrift reads where the ends align to the genome.
i <- (adriftReadAlignments$qSize - adriftReadAlignments$qEnd) <= opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned

if(sum(i) == 0) quitOnErorr('Error - no adrift read alignments remain after alignReads_genomeAlignment_adriftReadEnd_maxUnaligned filter.')


alignedReadIDsBeforeFilter <- n_distinct(adriftReadAlignments$readID)
adriftReadAlignments <- adriftReadAlignments[i,]

updateLog(paste0(sprintf("%.2f%%", (1 - n_distinct(adriftReadAlignments$readID) / alignedReadIDsBeforeFilter)*100), 
               ' of adrift reads removed because their alignments ended more than ',
               opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned,
               ' NTs from the end of reads.'))

# Select adrift reads with alignments that start at the beginning of reads.
i <- adriftReadAlignments$qStart <= opt$alignReads_genomeAlignment_adriftRead_maxStartPos

if(sum(i) == 0) quitOnErorr('Error - no adrift read alignments remain after alignReads_genomeAlignment_adriftRead_maxStartPos filter.')

alignedReadIDsBeforeFilter <- n_distinct(adriftReadAlignments$readID)
adriftReadAlignments <- adriftReadAlignments[i,]


# Find the intersection between retained anchor and adrift reads and limit read pairs to pairs
# where both mates aligned within parameters.

updateLog('Finding read pairs where both anchor and adrift read pair mates aligned well.')

# Select read pairs that have both a good anchor and adrift alignment.
i <- base::intersect(anchorReadAlignments$readID, adriftReadAlignments$readID)

if(length(i) == 0) quitOnErorr('Error -- no alignments remaing after filtering for paired alignments.')

anchorReadAlignments <- anchorReadAlignments[anchorReadAlignments$readID %in% i,]
adriftReadAlignments <- adriftReadAlignments[adriftReadAlignments$readID %in% i,]


# Add refGenome, vector, and flags.
anchorReadAlignments <- left_join(anchorReadAlignments, distinct(select(reads, readID, vectorFastaFile, flags)), by = 'readID')


# Expand predicted leaderSeq sequences by extending with delayed alignment sequences. 
anchorReadAlignments <- left_join(anchorReadAlignments, select(reads, readID, anchorReadSeq, leaderSeq), by = 'readID')
anchorReadAlignments$leaderSeq <- paste0(anchorReadAlignments$leaderSeq, substr(anchorReadAlignments$anchorReadSeq, 1, anchorReadAlignments$qStart))
anchorReadAlignments <- dplyr::select(anchorReadAlignments, -anchorReadSeq)


# Save anchor and adrift read alignments.
saveRDS(dplyr::distinct(anchorReadAlignments), file.path(opt$outputDir, opt$alignReads_outputDir, 'anchorReadAlignments.rds'), compress = opt$compressDataFiles)  
saveRDS(dplyr::distinct(adriftReadAlignments), file.path(opt$outputDir, opt$alignReads_outputDir, 'adriftReadAlignments.rds'), compress = opt$compressDataFiles) 


# Remove temp alignment files.
if(! opt$core_keepIntermediateFiles) unlink(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1'), recursive = TRUE)
if(! opt$core_keepIntermediateFiles) unlink(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2'), recursive = TRUE)

updateLog('alignReads completed.')

if(any(! incomingSamples %in% anchorReadAlignments$uniqueSample) & opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()

q(save = 'no', status = 0, runLast = FALSE) 
