# John K. Everett, PhD
# AAVengeR/demultipex.R
#
# This module demultiplexes paired-end reads based on barcode sequences found
# in the sampleData configuration file pointed to by the AAVengeR configuration file.

suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')

opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

createOuputDir()
if(! dir.exists(file.path(opt$outputDir))) dir.create(file.path(opt$outputDir))
if(! dir.exists(file.path(opt$outputDir, opt$demultiplex_outputDir))) dir.create(file.path(opt$outputDir, opt$demultiplex_outputDir))
if(! dir.exists(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp')))  dir.create(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp'))
if(! dir.exists(file.path(opt$outputDir, opt$demultiplex_outputDir, 'logs'))) dir.create(file.path(opt$outputDir, opt$demultiplex_outputDir, 'logs'))

# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, opt$demultiplex_outputDir, 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)



setOptimalParameters()
set.seed(1)

quitOnErorr <- function(msg){
  updateLog(msg)
  message(msg)
  message(paste0('See log for more details: ', opt$defaultLogFile))
  q(save = 'no', status = 1, runLast = FALSE) 
}


if(! file.exists(opt$demultiplex_sampleDataFile)) quitOnErorr('Error - the sample configuration file could not be found.')


# Read in sample data.
updateLog('Loading sample data.')
samples <- loadSamples()

# Throw errors if expected files are missing.
if(! file.exists(opt$demultiplex_adriftReadsFile)) quitOnErorr('Error - the adrift reads file could not be found.')
if(! file.exists(opt$demultiplex_anchorReadsFile)) quitOnErorr('Error - the index reads file could not be found.')
if(! file.exists(opt$demultiplex_index1ReadsFile)) quitOnErorr('Error - the anchor reads file could not be found.')

# Read in the I1 fastq file which will be used to determine the length of the data set
# and to determine if the reverse compliment of I1 should be used.

I1 <- ShortRead::readFastq(opt$demultiplex_index1ReadsFile)
dataSetLength <- length(I1)

updateLog(paste0(ppNum(dataSetLength), ' NTs in paired end data set.'))

if(opt$demultiplex_RC_I1_barcodes_auto){
  updateLog('Determining if I1 barcodes should be switched to RC.')
  
  d <- data.table(select(samples, subject, sample, replicate, index1Seq))
  i <- as.character(I1@sread)
  
  o <- rbindlist(lapply(split(d, 1:nrow(d)), function(x){
         x$barcodePercent <- sum(i %in% x$index1Seq)/length(i) * 100
         x$barcodePercentRC <- sum(i %in% as.character(Biostrings::reverseComplement(Biostrings::DNAString(x$index1Seq))))/length(i) * 100
         x
       }))
  
  r <- ifelse(sum(o$barcodePercent) > sum(o$barcodePercentRC), FALSE, TRUE)

  updateLog(paste0('Setting demultiplex_RC_I1_barcodes to ', r, '.'))

  opt$demultiplex_RC_I1_barcodes <- r
}

# Reverse compliment index1 sequences if requested.
if(opt$demultiplex_RC_I1_barcodes) samples$index1Seq <- as.character(reverseComplement(DNAStringSet(samples$index1Seq)))


# Create CPU clusters.
cluster <- makeCluster(opt$demultiplex_CPUs)
clusterSetRNGStream(cluster, 1)
clusterExport(cluster, c('opt', 'samples'))


# Create FASTQ file streamer objects.
index1.strm <- FastqStreamer(opt$demultiplex_index1ReadsFile, n = as.integer(opt$demultiplex_sequenceChunkSize))
anchor.strm <- FastqStreamer(opt$demultiplex_anchorReadsFile, n = as.integer(opt$demultiplex_sequenceChunkSize))
adrift.strm <- FastqStreamer(opt$demultiplex_adriftReadsFile, n = as.integer(opt$demultiplex_sequenceChunkSize))


# Stream chunks of FASTQ for I1, R1, and R2, write chunks to disk, then process in parallel.
n <- 1
k <- 1
processedReads <- 0

options(scipen = 999)

updateLog('Starting demultipling.')

repeat {
  
  updateLog(paste0('   processing data chunk ', k, '.'))
  
  index1.fq <- yield(index1.strm)
  if(length(index1.fq) == 0) break
  anchor.fq <- yield(anchor.strm)
  adrift.fq <- yield(adrift.strm)
  
  id <- paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''))
  
  writeFastq(index1.fq, file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0('index1_fastqChunk', '.', id)), compress = FALSE)
  writeFastq(anchor.fq, file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0('anchor_fastqChunk', '.', id)), compress = FALSE)
  writeFastq(adrift.fq, file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0('adrift_fastqChunk', '.', id)), compress = FALSE)
  
  processedReads <- processedReads + length(index1.fq)

  updateLog(paste0('   ', ppNum(processedReads), ' reads processed.'))
  
  if(n == opt$demultiplex_CPUs | length(index1.fq) < opt$demultiplex_sequenceChunkSize){
   
    updateLog('   Batch read limit reached, calling demultiplex().')
    
    o <- data.table(file = list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp'), pattern = 'fastqChunk'))
    o$n <- unlist(lapply(stringr::str_split(o$file, '\\.'), '[', 2))
    
    invisible(parLapply(cluster, split(o, o$n), demultiplex))
    
    updateLog(paste0('   ', 'processed ', sprintf("%.2f%%", (processedReads / dataSetLength)*100), ' of all reads.'))
    
    n <- 0
  }
  
  n <- n + 1
  k <- k + 1
}

updateLog('Completed processing read batches.')

# Collate demultiplexed chunks into a single data table.

updateLog('Colating data files.')

### counter <- 1
### totalRepSamples <- n_distinct(samples$uniqueSample)

reads <-  rbindlist(lapply(unique(samples$uniqueSample), function(x){
  ### message(counter, '/', totalRepSamples)
  ### counter <<- counter + 1

  f1 <- list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp'), pattern = paste0(x, '\\.[^\\.]+\\.anchorReads'), full.names = TRUE)
  f2 <- list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp'), pattern = paste0(x, '\\.[^\\.]+\\.adriftReads'), full.names = TRUE)
  f3 <- list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp'), pattern = paste0(x, '\\.[^\\.]+\\.randomAdriftReadIDs'), full.names = TRUE)
  if(length(f1) == 0 | length(f2) == 0 | length(f1) != length(f2)) return()

  updateLog(paste0('Colating reads for ', x, '.'))
  
  anchorReads <- Reduce('append', lapply(f1, readDNAStringSet))
  adriftReads <- Reduce('append', lapply(f2, readDNAStringSet))
  randomIDs   <- Reduce('append', lapply(f3, readDNAStringSet))
  
  closeAllConnections()
  
  r <- subset(samples, uniqueSample == x)
  c <- substr(r$adriftReadLinkerSeq, max(stringr::str_locate_all(r$adriftReadLinkerSeq, 'NNN')[[1]][,2])+1, nchar(r$adriftReadLinkerSeq))
  t <- as.character(reverseComplement(DNAString(substr(c, nchar(c) - 14, nchar(c)))))

  data.table(uniqueSample = x, readID = names(anchorReads), anchorReadSeq = as.character(anchorReads), adriftReadSeq = as.character(adriftReads), 
             adriftReadRandomID = as.character(randomIDs), adriftReadTrimSeq = t, adriftLinkerSeqEnd = nchar(r$adriftReadLinkerSeq),
             vectorFastaFile = r$vectorFastaFile, refGenome = r$refGenome, flags = r$flags)
}))



if(nrow(reads) == 0) quitOnErorr('Error - no reads were demultiplexed for any sample.')

updateLog('Clearing tmp files.')

invisible(file.remove(list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp'), full.names = TRUE)))

# Collect all the logs from the different computational nodes and create a single report.

updateLog('Colating log files.')

logReport <- bind_rows(lapply(list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'logs'), pattern = '*.logReport$', full.names = TRUE), function(f){
  read.table(f, header = TRUE, sep = '\t')
}))

logReport <- bind_rows(lapply(split(logReport, logReport$sample), function(x){
  o <- data.frame(lapply(2:length(x), function(y){
    if(all(is.na(x[,y]))){
      return(NA)
    } else {
      return(sum(x[,y], na.rm = TRUE))
    }
  }))
  
  names(o) <- names(x)[2:length(x)]
  bind_cols(data.frame(sample = x[1,1]), o)
})) %>% dplyr::arrange(demultiplexedReads)


updateLog('Writing attrition table.')

invisible(unlink(file.path(opt$outputDir, opt$demultiplex_outputDir, 'logs'), recursive = TRUE))
write.table(logReport, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'readAttritionTbl.tsv'))



# Expand read table with additional columns when appropriate.

updateLog('Expanding read table with meta data from sample table.')

if('anchorReadStartSeq' %in% names(samples)){
  reads <- left_join(reads, select(samples, uniqueSample, anchorReadStartSeq), by = 'uniqueSample')
}

if('leaderSeqHMM' %in% names(samples)){
  reads <- left_join(reads, select(samples, uniqueSample, leaderSeqHMM), by = 'uniqueSample')
}

if(! opt$processAdriftReadLinkerUMIs){
 updateLog('Setting read UMI sequences to poly-A because processAdriftReadLinkerUMIs is set to FALSE.')  
 reads$adriftReadRandomID <- 'AAAAAAAAAAAA'
 reads <- left_join(reads, select(samples, uniqueSample, adriftRead.linkerRandomID.start, adriftRead.linkerRandomID.end), by = 'uniqueSample')
 substr(reads$adriftReadSeq, reads$adriftRead.linkerRandomID.start, reads$adriftRead.linkerRandomID.end) <- 'AAAAAAAAAAAA'
 reads <- select(reads, -adriftRead.linkerRandomID.start, -adriftRead.linkerRandomID.end)
}

# Filter reads based on the requested demultiplex level.
# Read clustering is performed with CD-HIT-EST using the parameters passed in with demultiplex_mergeSimilarReadPairsParams.

if(opt$demultiplex_level == 'all'){
  reads$nDuplicateReads <- 0
} else if(opt$demultiplex_level == 'unique'){
  
  updateLog('Identifying unique read pairs.')
  
  cluster <- makeCluster(opt$demultiplex_CPUs)
  clusterSetRNGStream(cluster, 1)
  clusterExport(cluster, 'opt')
  
  reads <- data.table(reads)
  o <- rbindlist(parLapply(cluster, split(reads, reads$uniqueSample), function(x){
              suppressPackageStartupMessages(library(data.table))
              rbindlist(lapply(split(x, paste(x$anchorReadSeq, x$adriftReadSeq)), function(x2){
                          if(nrow(x2) > 1 ) x2 <- x2[order(x2$readID)]
                          x2$duplicated <- TRUE
                          x2$duplicatedRepID <- x2[1]$readID
                          x2$nDuplicateReads <- nrow(x2) - 1
                          x2[1]$duplicated <- FALSE
                          x2
                        }))
              }))
  
  stopCluster(cluster)
  
  reads <- dplyr::filter(o, duplicated == FALSE) %>% dplyr::select(-duplicated, -duplicatedRepID)
  dupReadsTable <- dplyr::filter(o, duplicated == TRUE) %>% dplyr::select(duplicatedRepID, readID)
  saveRDS(dupReadsTable, file =  file.path(opt$outputDir, opt$demultiplex_outputDir, 'dupReadsTable.rds'), compress = opt$compressDataFiles)
  rm(o)
  
} else if(opt$demultiplex_level == 'clustered'){
  
  updateLog('Identifying read pair clusters.')
  
  o <- rbindlist(lapply(split(reads, reads$uniqueSample), function(x){
    o <- DNAStringSet(paste0(x$anchorReadSeq, substr(x$adriftReadSeq, x$adriftLinkerSeqEnd+1, nchar(x$adriftReadSeq))))
    names(o) <- x$readID
    f <- tmpFile()
    writeFasta(o, file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', f))
    
    system(paste0("cd-hit-est -i ", file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', f), 
                  " -o ",  file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0(f, '.out')), " -T ", opt$demultiplex_CPUs, 
                  " ", opt$demultiplex_mergeSimilarReadPairsParams))
    
    r <- paste0(readLines(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0(f, '.out.clstr'))), collapse = '')
    
    invisible(file.remove(c(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', f),
                            file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0(f, '.out')),
                            file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0(f, '.out.clstr')))))
    
    m <- rbindlist(lapply(unlist(strsplit(r, '>Cluster')), function(x){
      e <- sub('>', '', unlist(stringr::str_extract_all(x, '>[^\\.]+')))
      
      if(length(e) > 0){
        repSeq <- sub('^>', '', stringr::str_extract(stringr::str_extract(x, '>[^\\.]+\\.+\\s+\\*'), '>[^\\.]+'))
        return(data.table(readID = e, rep = repSeq))
      } else {
        return(data.table())
      }
    })) 
    
    m <- left_join(m, select(reads, readID, adriftReadRandomID), by = 'readID')
    
    rbindlist(lapply(split(m, m$rep), function(x){
      rbindlist(lapply(split(x, x$adriftReadRandomID), function(x2){
        data.table(readID = ifelse(x$rep[1] %in% x2$readID, x$rep[1], x2[1,]$readID), nDuplicateReads = nrow(x2) - 1)
      }))
    }))
  }))
  
  reads <- subset(reads, readID %in% o$readID)
  reads <- left_join(reads, o, by = 'readID')
} else {
  # error
}


# Align demultiplexed reads to the reference genome with bwa2-mem and only 
# retain reads that align. These alignments will not be as accurate and blat's
# alignments but removing reads that blat will likely not be able to align 
# greatly increase the speed of the pipeline.

# Here we are aligning reads before collapsing reads to unique or clustered groupings.
# This puts more of a burden on bwa2 but helps with tracking reads since running this 
# alignment after grouping the reads would cause groups of reads to fall out of the analysis 
# and disrupt read counts. 

reads$quickFilterStartPos <- NA

if(opt$demultiplex_quickAlignFilter){
  readsLengthPreFilter <- n_distinct(reads$readID)
  
  updateLog(paste0('Prefiltering', ppNum(readsLengthPreFilter), ' reads with bwa2.'))
  
  reads <- rbindlist(lapply(split(reads, reads$refGenome), function(x){
    
    o <- DNAStringSet(substr(x$adriftReadSeq, x$adriftLinkerSeqEnd+1, nchar(x$adriftReadSeq)))
    names(o) <- x$readID
    writeXStringSet(o, file.path(opt$outputDir, opt$demultiplex_outputDir, 'adriftReads.fasta'))
    
    system(paste0('bwa-mem2 mem -t ', opt$demultiplex_CPUs, ' -c 100000 -a ', file.path(opt$softwareDir, 'data', 'referenceGenomes', 'bwa2', x$refGenome[1]), ' ',
                  file.path(opt$outputDir, opt$demultiplex_outputDir, 'adriftReads.fasta'), ' > ',
                  file.path(opt$outputDir, opt$demultiplex_outputDir, 'adriftReads.sam')))
    
    system(paste0(file.path(opt$softwareDir, 'bin', 'sam2psl.py'), ' -i ', 
                  file.path(opt$outputDir, opt$demultiplex_outputDir, 'adriftReads.sam'), ' -o ',
                  file.path(opt$outputDir, opt$demultiplex_outputDir, 'adriftReads.psl')))
    
    a <- parseBLAToutput(file.path(opt$outputDir, opt$demultiplex_outputDir, 'adriftReads.psl'), convertToBlatPSL = TRUE)
    a <- dplyr::filter(a, alignmentPercentID >= opt$demultiplex_quickAlignFilter_minPercentID, tNumInsert <= 1, 
                       qNumInsert <= 1, tBaseInsert <= 2, qBaseInsert <= 2, qStart <= 5, matches >= opt$demultiplex_quickAlignFilter_minMatches)
    
    invisible(file.remove(file.path(opt$outputDir, opt$demultiplex_outputDir, 'adriftReads.sam'), 
                          file.path(opt$outputDir, opt$demultiplex_outputDir, 'adriftReads.fasta'),
                          file.path(opt$outputDir, opt$demultiplex_outputDir, 'adriftReads.psl')))
    
    x <- subset(x, readID %in% a$qName)
    
    o <- DNAStringSet(substr(x$anchorReadSeq, opt$demultiplex_quickAlignFilter_minEstLeaderSeqLength, nchar(x$anchorReadSeq)))
    names(o) <- x$readID
    writeXStringSet(o, file.path(opt$outputDir, opt$demultiplex_outputDir, 'anchorReads.fasta'))
    
    # -B -O flags makes bwa-mem2 more tolerant of mismatches near the ends of alignments.
    system(paste0('bwa-mem2 mem -B 2 -O 4 -t ', opt$demultiplex_CPUs, ' -c 100000 -a ', file.path(opt$softwareDir, 'data', 'referenceGenomes', 'bwa2', x$refGenome[1]), ' ',
                  file.path(opt$outputDir, opt$demultiplex_outputDir, 'anchorReads.fasta'), ' > ',
                  file.path(opt$outputDir, opt$demultiplex_outputDir, 'anchorReads.sam')))
    
    system(paste0(file.path(opt$softwareDir, 'bin', 'sam2psl.py'), ' -i ', 
                  file.path(opt$outputDir, opt$demultiplex_outputDir, 'anchorReads.sam'), ' -o ',
                  file.path(opt$outputDir, opt$demultiplex_outputDir, 'anchorReads.psl')))
    
    p <- readr::read_tsv(file.path(opt$outputDir, opt$demultiplex_outputDir, 'anchorReads.psl'), col_names = FALSE)
    invisible(file.remove(file.path(opt$outputDir, opt$demultiplex_outputDir, 'anchorReads.psl')))
    
    p$X12 <- ifelse(p$X9 == '-', p$X11 - p$X13, p$X12)
    p$X13 <- ifelse(p$X9 == '-', p$X11 - p$X12, p$X13)
    readr::write_tsv(p, file.path(opt$outputDir, opt$demultiplex_outputDir, 'anchorReads.psl'), col_names = FALSE)
    
    a <- parseBLAToutput(file.path(opt$outputDir, opt$demultiplex_outputDir, 'anchorReads.psl'))
    
    a <- dplyr::filter(a, alignmentPercentID >= opt$demultiplex_quickAlignFilter_minPercentID, tNumInsert <= 1, 
                       qNumInsert <= 1, tBaseInsert <= 2, qBaseInsert <= 2, matches >= opt$demultiplex_quickAlignFilter_minMatches)
    
    invisible(file.remove(file.path(opt$outputDir, opt$demultiplex_outputDir, 'anchorReads.sam'), 
                          file.path(opt$outputDir, opt$demultiplex_outputDir, 'anchorReads.fasta'),
                          file.path(opt$outputDir, opt$demultiplex_outputDir, 'anchorReads.psl')))
    
    x <- subset(x, readID %in% a$qName)
    
    minAlnStarts <- group_by(a, qName) %>% 
                    slice_max(matches, with_ties = TRUE) %>% 
                    slice_min(qStart, with_ties = FALSE) %>%
                    select(qName, qStart) %>% 
                    ungroup() %>%
                    mutate(quickFilterStartPos = qStart + opt$demultiplex_quickAlignFilter_minEstLeaderSeqLength) %>%
                    select(-qStart)
    
    x$quickFilterStartPos <- NULL
    left_join(x, minAlnStarts, by = c('readID' = 'qName'))
  }))
  
  d <- sprintf("%.2f%%", (1 - (n_distinct(reads$readID)/readsLengthPreFilter))*100)
  
  updateLog(paste0('Prefiltering done.', ppNum(d), 'reads removed.'))
}

# FASTQ export if requested in the configuration file.
if(opt$demultiplex_exportFASTQ){
  
  updateLog('Exporting FASTQ.')
  
  dir.create(file.path(opt$outputDir, opt$demultiplex_outputDir, 'fastq'))
  
  invisible(lapply(c('I1', 'R1', 'R2'), function(a){
    
    if(a == 'I1'){
      f <- opt$demultiplex_index1ReadsFile
    } else if (a == 'R1'){
      f <- opt$demultiplex_adriftReadsFile
    } else {
      f <- opt$demultiplex_anchorReadsFile
    }
    
    r <- readFastq(f)
    ids <- sub('\\s+.+$', '', as.character(r@id))
  
    invisible(lapply(split(reads, reads$uniqueSample), function(x){
      writeFastq(r[match(x$readID, ids)],
                 file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'fastq', paste0(x$uniqueSample[1], '.', a, '.fastq.gz')),
                 compress = TRUE)
    }))
  }))
}


if(opt$demultiplex_requirePostUmiLinker){
  # Common Bushman linker: CTCCGCTTAAGGGACT
  
  updateLog('Applying post UMI filter because demultiplex_requirePostUmiLinker is set to TRUE.')
  
  reads <- rbindlist(lapply(split(reads, reads$uniqueSample), function(x){
    d <- subset(samples, uniqueSample == x$uniqueSample[1])
    s <- substr(x$adriftReadSeq, d$adriftRead.linkerRandomID.end, nchar(d$adriftReadLinkerSeq) + 1)   # expand region +1 on both ends.
    e <- substr(d$adriftReadLinkerSeq, d$adriftRead.linkerRandomID.end + 1, nchar(d$adriftReadLinkerSeq))
    i <- vcountPattern(e, DNAStringSet(s), max.mismatch = opt$demultiplex_requirePostUmiLinker_maxMismatch) == 1
    
    msg <- paste0(x$uniqueSample[1], ' - ', sprintf("%.2f%%", (sum(i == FALSE)/nrow(x))*100), ' reads missing post UMI linker.')
    updateLog(msg)
   
    x[i]
  }))
}


# Replicate reassingment if requested.
if(file.exists(opt$demultiplex_replicateMergingInstructions)){
  updateLog('Found replicate merging instruction file.')
  
  f <- readr::read_tsv(opt$demultiplex_replicateMergingInstructions)
  names(f) <- c('uniqueSample', 'uniqueSample2')
  
  if(any(f$uniqueSample %in% reads$uniqueSample)){
    updateLog('Merging replicates.')
    reads <- left_join(reads, f, by = 'uniqueSample')
    reads$uniqueSample <- ifelse(is.na(reads$uniqueSample2), reads$uniqueSample, reads$uniqueSample2)
    reads$uniqueSample2 <- NULL
  } else {
    updateLog('None of the replicates in the merging instructions were found in the demultiplexed data.')
  }
}

updateLog('Writing output files.')

# Save demultiplexed reads and clean up.
reads$seqRunID <- opt$demultiplex_seqRunID
saveRDS(reads, file =  file.path(opt$outputDir, opt$demultiplex_outputDir, 'reads.rds'), compress = opt$compressDataFiles)
invisible(file.remove(list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp'), full.names = TRUE)))

updateLog('Demultiplex completed.')

q(save = 'no', status = 0, runLast = FALSE) 
