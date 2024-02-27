# AAVengeR/anchorReadRearrangements.R
# John K. Everett, Ph.D.
# 
# This script accepts the a demultiplexed read table from the demultiplex module
# and identifies instances of vector rearrangements. This analysis focus on anchor
# reads and it is often applied to sequencing experiments where anchor reads lengths
# are prioritized over adrift read lengths.

suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))

# Read the configuration file and set additional parameters.
configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')

opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

createOuputDir()
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir))
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs'))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs'))
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'))
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'repLeaderSeqStructs'))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'repLeaderSeqStructs'))
invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), full.names = TRUE)))

setMissingOptions()
setOptimalParameters()
set.seed(1)


# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)
write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)


quitOnErorr <- function(msg){
  updateLog(msg)
  message(msg)
  message(paste0('See log for more details: ', opt$defaultLogFile))
  q(save = 'no', status = 1, runLast = FALSE)
}

updateLog('Reading sample data.')

cluster <- makeCluster(opt$anchorReadRearrangements_CPUs)
clusterSetRNGStream(cluster, 1)

if(! file.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_inputFile))) quitOnErorr('Error - the input data file does not exist.')

# Read in the reads table.
updateLog('Reading in demultiplexed reads.')
reads <- readRDS(file.path(opt$outputDir, opt$anchorReadRearrangements_inputFile))

reads <- reads[nchar(reads$anchorReadSeq) >= opt$anchorReadRearrangements_minAnchorReadLength,]

# Hot fix for SampleSheet errors during demultiplexing.
#
o <- readr::read_csv('correctedMetaData.csv', col_names = TRUE)
o$uniqueSample <- paste0(o$trial, '~', o$subject, '~', o$sample, '~', o$replicate)
reads <- rbindlist(lapply(split(reads, reads$uniqueSample), function(x){
           x$vectorFastaFile <- subset(o, uniqueSample == x$uniqueSample[1])$vectorFastaFile
           x
         }))
reads$vectorFastaFile <- sub('BushmanAAVcontrols.fasta', 'BushmanAAVcontrolsPlasmidLargestRemnant.fasta', reads$vectorFastaFile)

if(nrow(reads) == 0) quitOnErorr('Error - the input table contained no rows.')

reads <- select(reads, uniqueSample, readID, anchorReadSeq, adriftReadTrimSeq, vectorFastaFile)

updateLog('Limiting reads to select anchor read start sequences.')

reads$sample <- sub('~\\d+$', '', reads$uniqueSample)

save.image('~/anchorReadRearrangements_save1.RData')

if(! all(unique(reads$vectorFastaFile) %in%  names(opt$anchorReadRearrangements_expectedSeqs))){
  updateLog('Error -- all the vector names found in the read data are not defined in the anchorReadRearrangements_expectedSeqs section of the config file.')
  updateLog(paste0('These vectors are missing: ', paste0(unique(reads$vectorFastaFile)[! unique(reads$vectorFastaFile) %in%  names(opt$anchorReadRearrangements_expectedSeqs)], collapse = ',')))
  q()
}

nReadsPreFilter <- nrow(reads)

reads <- rbindlist(lapply(split(reads, reads$vectorFastaFile), function(x){
           v <- unique(x$vectorFastaFile)
           updateLog(paste0('Limiting ', ppNum(nrow(x)), ' reads for vector:  ', v))

           testSeqs <- substr(opt[['anchorReadRearrangements_expectedSeqs']][[v]], 1, opt$anchorReadRearrangements_vectorLeaderSeqFilterLength)

           rbindlist(lapply(testSeqs, function(k){
             updateLog(paste0('  Limiting reads for this vector to those starting with ', k))

             r <- x[vcountPattern(k,
                                  DNAStringSet(substr(x$anchorReadSeq, 1, nchar(k))),
                                  max.mismatch = opt$anchorReadRearrangements_vectorLeaderSeqFilterMaxMismatch) == 1]
             updateLog(paste0('    ', ppNum(nrow(r)), ' reads matched the leader sequence.'))
             r
           }))
         }))

if(nrow(reads) == 0) quitOnErorr('Error - not reads matched any of the provided leader sequences.')

updateLog(paste0(sprintf("%.2f%%", (nrow(reads) / nReadsPreFilter)*100), ' reads retained.'))

# Trim anchor read over-reading with cutadapt using the RC of the common linker in adrift reads.
# The trim sequence is defined in the reads table created by demultiplex.R.

updateLog('Trimming anchor read over-reading.')

if(tolower(opt$anchorReadRearrangements_additionalAnchorReadOverReadingSeqs) != 'none'){
  opt$anchorReadRearrangements_additionalAnchorReadOverReadingSeqs <- paste0(' -a ', paste0(unlist(strsplit(opt$anchorReadRearrangements_additionalAnchorReadOverReadingSeqs, ',')), collapse = ' -a '))
} else{
  opt$anchorReadRearrangements_additionalAnchorReadOverReadingSeqs = ''
}

clusterExport(cluster, c('opt', 'buildRearrangementModel'))

reads <- data.table::rbindlist(parLapply(cluster, split(reads, dplyr::ntile(1:nrow(reads), opt$anchorReadRearrangements_CPUs)), function(x){
  source(file.path(opt$softwareDir, 'lib.R'))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(Biostrings))

  # Split reads in this CPU chunk by their adapter sequences.
  data.table::rbindlist(lapply(split(x, x$adriftReadTrimSeq), function(y){

    # Write out reads in this chunk to a tmp file.
    f <- file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp',  tmpFile())
    o <- DNAStringSet(y$anchorReadSeq)
    names(o) <- y$readID
    Biostrings::writeXStringSet(o, f)

    # Use cut adapt to trim anchor read over-reading.
    system(paste0('cutadapt -e ', opt$anchorReadRearrangements_cutAdaptErrorRate, ' -a ', y$adriftReadTrimSeq[1], ' ', opt$anchorReadRearrangements_additionalAnchorReadOverReadingSeqs, ' ', f, ' > ', paste0(f, '.cutAdapt')), ignore.stderr = TRUE)

    # Read in the cutadapt trimmed sequences.
    t <- readDNAStringSet(paste0(f, '.cutAdapt'))
    invisible(file.remove(f, paste0(f, '.cutAdapt')))

    # Remove reads that were trimmed too short.
    t <- t[width(t) >= opt$anchorReadRearrangements_minAnchorReadLength]
    if(length(t) == 0) return(data.table())

    # Limit incoming reads to those found in the cutAdapt output.
    o <- subset(y, readID %in% names(t))
    trimmed <- data.table(readID = names(t), anchorReadSeq2 = as.character(t))

    # Attach the trimmed sequences to the input sequence labeled as anchorReadSeq2.
    data.table(left_join(o, trimmed, by = 'readID'))
  }))
}))


# Report trimming stats.
trimReport <- mutate(reads, sample = sub('~\\d+$', '', uniqueSample)) %>%
  group_by(sample) %>%
  summarise(reads = n_distinct(readID), percentTrimmed = sprintf("%.2f%%", (sum(anchorReadSeq != anchorReadSeq2)/n())*100)) %>%
  ungroup() %>%
  readr::write_tsv(file = file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'commonLinkerTrimReport.tsv'), append = FALSE, col_names = TRUE)

# Switch original sequences for trimmed sequences.
reads$anchorReadSeq <- NULL
reads <- dplyr::rename(reads, anchorReadSeq = anchorReadSeq2)

nReadsPreFilter <- n_distinct(reads$readID)

reads <- bind_rows(lapply(split(reads, reads$vectorFastaFile), function(x){
           f <- tmpFile()
           x$lastNTs <- substr(x$anchorReadSeq, nchar(x$anchorReadSeq) - (opt$anchorReadRearrangements_readEndAlignmentTestLength - 1), nchar(x$anchorReadSeq))

           system(paste0('/home/ubuntu/software/faToTwoBit ',
                         file.path(opt$softwareDir, 'data', 'vectors', x$vectorFastaFile[1]), ' ',
                         file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', f)))

           d <- DNAStringSet(unique(x$lastNTs))
           names(d) <- as.character(d)
           writeXStringSet(d,  file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.fasta')))

           system(paste0('blat ', file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', f), ' ',
                         file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.fasta')), ' ',
                         file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.psl')),
                         ' -tileSize=6 -stepSize=3 -repMatch=3000 -out=psl -t=dna -q=dna -minScore=0 -minIdentity=0 -noHead -noTrimA'))

          o <- parseBLAToutput(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.psl')))
          o <- subset(o, matches >= opt$anchorReadRearrangements_readEndAlignmentTestMinMatch)

          invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), pattern = f, full.names = TRUE)))

          # add alignment failure test...
          subset(x, lastNTs %in% o$qName)
        }))

updateLog(paste0(sprintf("%.2f%%", (n_distinct(reads$readID) / nReadsPreFilter)*100), ' reads remain after removing reads whoes tails do not align to the vector.'))

save.image('~/anchorReadRearrangements_save2.RData')

# load('~/anchorReadRearrangements_save2.RData')
# opt <- yaml::read_yaml('config.yml')
# source(file.path(opt$softwareDir, 'lib.R'))
# opt$defaultLogFile <- file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'log')
# logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
# write(logo, opt$defaultLogFile, append = FALSE)
# write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)
# cluster <- makeCluster(opt$anchorReadRearrangements_CPUs)
# clusterSetRNGStream(cluster, 1)
# clusterExport(cluster, c('opt', 'buildRearrangementModel'))

samplesProcessed <- 1

r <- bind_rows(lapply(split(reads, paste(reads$vectorFastaFile, reads$sample)), function(x){
       updateLog(paste0('Processing sample ', samplesProcessed, ' / ', n_distinct(reads$sample), ' - sample name: ', x$sample[1]))
       samplesProcessed <<- samplesProcessed  + 1
   
       #browser()
       #x.save <- x
       #x <- x[1:12,]
       #x$anchorReadSeq <- as.character(tt)
       #x$readID <- names(tt)
       
       d <- DNAStringSet(x$anchorReadSeq)
       names(d) <- x$readID
       
       d2 <- DNAStringSet(opt[['anchorReadRearrangements_expectedSeqs']][[unique(x$vectorFastaFile)]])
       d2 <- subseq(d2, 1, max(nchar(x$anchorReadSeq)) + 5)
       names(d2) <- paste0('expectedSeq', 1:length(d2))
       
       message('vector: ', unique(x$vectorFastaFile), ' nSeqs: ', length(d2))
       
       writeXStringSet(c(d, d2), 'tmpDNAseq.fasta')
       
       system(paste0("cd-hit-est -i tmpDNAseq.fasta -o tmpDNAseq -T ", 
                     opt$anchorReadRearrangements_CPUs, ' ',
                     opt$anchorReadRearrangements_clusterParams))
       
       r <- paste0(readLines('tmpDNAseq.clstr'), collapse = '')
       
       invisible(file.remove(list.files(pattern = 'tmpDNAseq', full.names = TRUE)))
       
       o <- bind_rows(lapply(unlist(strsplit(r, '>Cluster')), function(x){
              e <- sub('>', '', unlist(stringr::str_extract_all(x, '>[^\\.]+')))
         
              if(length(e) > 0){
                repSeq <- sub('^>', '', stringr::str_extract(stringr::str_extract(x, '>[^\\.]+\\.+\\s+\\*'), '>[^\\.]+'))
                return(tibble(readID = e, rep = repSeq))
              } else {
                return(tibble())
              }
            }))
       
       # Remove expected sequence(s).
       o <- o[! grepl('expectedSeq', o$readID),]
       
       expectedSeqReadsIDs <- unique(o[grepl('expectedSeq', o$rep),]$readID)
       expectedSeqReadsNTs <- sum(nchar(subset(reads, readID %in% expectedSeqReadsIDs)$anchorReadSeq))
       
       unexpectedSeqReadsIDs <- unique(o[! grepl('expectedSeq', o$rep),]$readID)
       unexpectedSeqReadsNTs <- sum(nchar(subset(reads, readID %in% unexpectedSeqReadsIDs)$anchorReadSeq))
       
       allReadNTs <- sum(nchar(x$anchorReadSeq))
       allReadIDs <- unique(o$readID)
       
       altStructBreaks <- NA 
       
       o <- o[! grepl('expectedSeq', o$rep),]
       
       if(nrow(o) > 0){
         # Make a blast database for the vector sequence.
         files <- list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs'), full.names = TRUE)
         if(length(files) > 0) invisible(file.remove(files))
         system(paste0('makeblastdb -in ', file.path(opt$softwareDir, 'data', 'vectors', x$vector[1]), 
                       ' -dbtype nucl -out ', file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
       
         repReads <- DNAStringSet(x[x$readID %in% o$rep,]$anchorReadSeq)
         names(repReads) <- x[x$readID %in% o$rep,]$readID
       
         writeXStringSet(repReads, 'repReads.fasta')
       
         # Align the chunk to the vector sequence.  repReads.fasta
         #
         # (!) Here we increase the penaly from the default -3 to -4 to prevent runs of mismatches to come through if followed by a string of matches.
         #
         system(paste0('blastn -penalty -4 -max_target_seqs 1000 -gapopen 10 -gapextend 5 -dust no -soft_masking false -word_size 5 -evalue 100 -outfmt 6 -query repReads.fasta -db ',
                       file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs', 'd'),
                       ' -num_threads ', opt$anchorReadRearrangements_CPUs  ,' -out repReads.blast'), ignore.stdout = TRUE, ignore.stderr = TRUE)
       
         invisible(file.remove('repReads.fasta'))
       
         if(file.info('repReads.blast')$size != 0){
           # Parse blastn result.
           b <- read.table('repReads.blast', sep = '\t', header = FALSE)
           invisible(file.remove('repReads.blast'))
           names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
           b$alignmentLength <- b$qend - b$qstart + 1
           b$strand <- ifelse(b$sstart > b$send, '-', '+')
         
           w <- tibble(qname = names(d), qlen = width(d))
           b <- left_join(b, w, by = 'qname')
         
           b <- dplyr::filter(b, pident >= opt$anchorReadRearrangements_minAlignmentPercentID, alignmentLength >= opt$anchorReadRearrangements_minAlignmentLength)
         
           # Create a grouping index for read IDs.
           b$i <- group_by(b, qname) %>% group_indices()
         
           # Create a CPU splitting vector that would not separate the read ID indices. 
           k <- tibble(i2 = 1:n_distinct(b$i))
           k$n <- ntile(1:nrow(k), opt$anchorReadRearrangements_CPUs)
         
           # Bind the CPU splitting vector.
           b <- left_join(b, k, by = c('i' = 'i2'))
         
           r <- rbindlist(parallel::parLapply(cluster, split(b, b$n), blast2rearangements, opt$anchorReadRearrangements_maxMissingTailNTs, opt$anchorReadRearrangements_minAlignmentLength))
         
           r <- left_join(r, data.frame(table(o$rep)), by = c('qname' = 'Var1'))
           names(r) <- c('clusterRepReadID', 'rearrangement', 'readsInCluster')
           r <- arrange(r, desc(readsInCluster))
         
           r <- mutate(r, sample = x$sample[1], .before = 'clusterRepReadID')
           readr::write_tsv(r, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'repLeaderSeqStructs.tsv'), append = TRUE)
           
           altStructBreaks <- sum(stringr::str_count(r$rearrangement, ';'))
         }
       }   
       
       tibble(sample                   = x$sample[1], 
              totalReads               = n_distinct(allReadIDs), 
              totalNTs                 = allReadNTs,
              expectedStructReads      = n_distinct(expectedSeqReadsIDs),
              expectedStructNTs        = expectedSeqReadsNTs,
              altStructs               = n_distinct(o$rep), 
              altStructReads           = n_distinct(unexpectedSeqReadsIDs),
              altStructNTs             = unexpectedSeqReadsNTs,
              altStructBreaks          = altStructBreaks,
              percentAltStructReads    = (n_distinct(unexpectedSeqReadsIDs)/n_distinct(allReadIDs))*100,
              percentAltStructNTs      = (unexpectedSeqReadsNTs/allReadNTs)*100,
              altStructsPer1KB         = altStructs / (totalNTs/1000),
              altStructsPer10KB        = altStructs / (totalNTs/10000),
              altStructBreaksPer1KB    = altStructBreaks / (totalNTs/5000),    
              altStructBreaksPer10KB   = altStructBreaks / (totalNTs/10000))
}))

# r$n <- ifelse(grepl('NoTemp', r$sample), 0, ifelse(grepl('Positive', r$sample), 2, ifelse(grepl('Bogus', r$sample), 3, ifelse(grepl('TAK', r$sample), 4, 5))))

saveRDS(r, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'result.rds'))
openxlsx::write.xlsx(r, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'result.xlsx'))

