# John K. Everett, PhD
# AAVengeR/demultipex.R
#
# This module demultiplexes paired-end reads based on barcode sequences found
# in the sampleData configuration file pointed to by the AAVengeR configuration file.

library(ShortRead)
library(readr)
library(parallel)
library(lubridate)
library(dplyr)
library(data.table)
library(dtplyr)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))
setOptimalParameters()

if(! file.exists(opt$demultiplex_sampleDataFile)){
  write(c(paste(now(), '   Error - the sample configuration file could not be found')), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}

# Create required directory structure.
if(! dir.exists(file.path(opt$outputDir))) dir.create(file.path(opt$outputDir))
if(! dir.exists(file.path(opt$outputDir, opt$demultiplex_outputDir))) dir.create(file.path(opt$outputDir, opt$demultiplex_outputDir))
if(! dir.exists(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp')))  dir.create(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp'))
if(! dir.exists(file.path(opt$outputDir, opt$demultiplex_outputDir, 'logs'))) dir.create(file.path(opt$outputDir, opt$demultiplex_outputDir, 'logs'))

# Read in sample data.
write(c(paste(now(), '   Loading sample data')), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = FALSE)
samples <- loadSamples()

# Throw errors if expected files are missing.
if(! file.exists(opt$demultiplex_adriftReadsFile)){
  write(c(paste(now(), 'Error - the adrift reads file could not be found')), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}

if(! file.exists(opt$demultiplex_anchorReadsFile)){
    write(c(paste(now(), '   Error - the index reads file could not be found')), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE) 
}

if(! file.exists(opt$demultiplex_index1ReadsFile)){
  write(c(paste(now(), '   Error - the anchor reads file could not be found')), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
} 



# Read in the I1 fastq file which will be used to determine the length of the data set
# and to determine if the reverse compliment of I1 should be used.

I1 <- ShortRead::readFastq(opt$demultiplex_index1ReadsFile)
dataSetLength <- length(I1)

if(opt$demultiplex_RC_I1_barcodes_auto){
  write(c(paste(now(), '   Determining if I1 barcodes should be switched to RC.')), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
  d <- data.table(select(samples, subject, sample, replicate, index1Seq))
  i <- as.character(I1@sread)
  
  o <- rbindlist(lapply(split(d, 1:nrow(d)), function(x){
         x$barcodePercent <- sum(i %in% x$index1Seq)/length(i) * 100
         x$barcodePercentRC <- sum(i %in% as.character(Biostrings::reverseComplement(Biostrings::DNAString(x$index1Seq))))/length(i) * 100
         x
       }))
  
  r <- ifelse(sum(o$barcodePercent) > sum(o$barcodePercentRC), FALSE, TRUE)
  write(paste(now(), '   Setting demultiplex_RC_I1_barcodes to', r), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
  opt$demultiplex_RC_I1_barcodes <- r
}

# Reverse compliment index1 sequences if requested.
if(opt$demultiplex_RC_I1_barcodes) samples$index1Seq <- as.character(reverseComplement(DNAStringSet(samples$index1Seq)))


# Create CPU clusters.
cluster <- makeCluster(opt$demultiplex_CPUs)
clusterExport(cluster, c('opt', 'samples'))


# Create FASTQ file streamer objects.
index1.strm <- FastqStreamer(opt$demultiplex_index1ReadsFile, n = as.integer(opt$demultiplex_sequenceChunkSize))
anchor.strm <- FastqStreamer(opt$demultiplex_anchorReadsFile, n = as.integer(opt$demultiplex_sequenceChunkSize))
adrift.strm <- FastqStreamer(opt$demultiplex_adriftReadsFile, n = as.integer(opt$demultiplex_sequenceChunkSize))



# Stream chunks of FASTQ for I1, R1, and R2, write chunks to disk, then process in parallel.
n <- 1
k <- 1
processedReads <- 0

repeat {
  message('n: ', n, ' k: ', k)
  index1.fq <- yield(index1.strm)
  if(length(index1.fq) == 0) break
  anchor.fq <- yield(anchor.strm)
  adrift.fq <- yield(adrift.strm)
  
  id <- paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''))
  
  writeFastq(index1.fq, file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0('index1_fastqChunk', '.', id)), compress = FALSE)
  writeFastq(anchor.fq, file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0('anchor_fastqChunk', '.', id)), compress = FALSE)
  writeFastq(adrift.fq, file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0('adrift_fastqChunk', '.', id)), compress = FALSE)
  
  processedReads <- processedReads + length(index1.fq)
  
  if(n == opt$demultiplex_CPUs | length(index1.fq) < opt$demultiplex_sequenceChunkSize){
    o <- data.table(file = list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp'), pattern = 'fastqChunk'))
    o$n <- unlist(lapply(stringr::str_split(o$file, '\\.'), '[', 2))
    
    invisible(parLapply(cluster, split(o, o$n), demultiplex))
    #invisible(lapply(split(o, o$n), demultiplex))
    
    message('Processed ', sprintf("%.2f%%", (processedReads / dataSetLength)*100), ' reads')
    
    n <- 0
  }
  
  n <- n + 1
  k <- k + 1
}



# Collate demultiplexed chunks into a single data table.

write(paste(now(), '   Colating data files.'), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
counter <- 1
totalRepSamples <- n_distinct(samples$uniqueSample)

reads <-  rbindlist(lapply(unique(samples$uniqueSample), function(x){
  message(counter, '/', totalRepSamples)
  counter <<- counter + 1
  
  f1 <- list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp'), pattern = paste0(x, '\\.[^\\.]+\\.anchorReads'), full.names = TRUE)
  f2 <- list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp'), pattern = paste0(x, '\\.[^\\.]+\\.adriftReads'), full.names = TRUE)
  f3 <- list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp'), pattern = paste0(x, '\\.[^\\.]+\\.randomAdriftReadIDs'), full.names = TRUE)
  if(length(f1) == 0 | length(f2) == 0 | length(f1) != length(f2)) return()

  write(paste0(now(), '    Colating reads for ', x), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
  
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

if(nrow(reads) == 0){
  write(c(paste(now(), '   Error - no reads were demultiplexed for any sample.')), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}

write(paste(now(), '   Clearing tmp files.'), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
invisible(file.remove(list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp'), full.names = TRUE)))



# Collect all the logs from the different computational nodes and create a single report.

write(paste(now(), '   Colating log files.'), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)

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

write(paste(now(), '   Writing attrition table.'), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
invisible(unlink(file.path(opt$outputDir, opt$demultiplex_outputDir, 'logs'), recursive = TRUE))
write.table(logReport, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'readAttritionTbl.tsv'))



# Expand read table with additional columns when appropriate.

if('anchorReadStartSeq' %in% names(samples)){
  reads <- left_join(reads, select(samples, uniqueSample, anchorReadStartSeq), by = 'uniqueSample')
}

if('leaderSeqHMM' %in% names(samples)){
  reads <- left_join(reads, select(samples, uniqueSample, leaderSeqHMM), by = 'uniqueSample')
}



# If a database group is defined in the configuration file, make connect to the database 
# and save the contents of the sample details file save a record of linker and barcodes. 

if('databaseGroup' %in% names(opt)){
  library(RMariaDB)
  
  write(paste(now(), '   Uploading to database.'), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
  
  conn <- tryCatch({
    dbConnect(RMariaDB::MariaDB(), group = opt$databaseGroup)
  },
  error = function(cond) {
    write(c(paste(now(), '   Error - could not connect to the database.')), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE) 
  })
  
  invisible(lapply(split(samples, samples$uniqueSample), function(x){
    dbExecute(conn, paste0("delete from demultiplex where trial='", x$trial, "' and subject='", x$subject, "' and sample='", x$sample, "' and replicate='", x$replicate, "'"))
  }))
  
  r <- unlist(lapply(split(samples, 1:nrow(samples)), function(x){
        dbExecute(conn, "insert into demultiplex values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", 
                   params = list(x$trial, x$subject, x$sample, x$replicate, x$refGenome, lpe(x$vectorFastaFile), x$flags, 
                                 x$index1Seq, x$adriftReadLinkerSeq, ifelse('leaderSeqHMM' %in% names(x), lpe(x$leaderSeqHMM), 'NA'),
                                 opt$demultiplex_seqRunID))
        }))
  
  if(any(r == 0)){
    write(c(paste(now(), '   Error - could not upload all sample records to the database.')), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE) 
  }
  
  dbDisconnect(conn)
}



# If UMI processing is disabled, set all random UMI sequences to poly-A and updated
# adrift reads by replacing the random UMI sequences with poly-A.

if(! opt$demultiplex_processAdriftReadLinkerUMIs){
  reads$adriftReadRandomID <- 'AAAAAAAAAAAA'
  
  reads <- rbindlist(lapply(split(reads, reads$uniqueSample), function(x){
             start <- subset(samples, uniqueSample == x$uniqueSample[1])$adriftRead.linkerRandomID.start
             stop  <- subset(samples, uniqueSample == x$uniqueSample[1])$adriftRead.linkerRandomID.end
             substr(x$adriftReadSeq, start, stop) <- 'AAAAAAAAAAAA'
             x
            }))
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
  
  write(paste(now(), '   Prefiltering', readsLengthPreFilter, 'reads with bwa2.'), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
  
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
  
  write(paste(now(), '   Prefiltering done.', d, 'reads removed.'), file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), append = TRUE)
}





# Filter reads based on the requested demultiplex level.
# All reads, unique reads (R1 [minus linker] + R2), or clustered reads can be retained.
# Read clustering is performed with CD-HIT-EST using the parameters passed in with demultiplex_mergeSimilarReadPairsParams.

if(opt$demultiplex_level == 'all'){
  reads$nDuplicateReads <- 0
} else if(opt$demultiplex_level == 'unique'){
  reads$adriftReadSeq2 <- substr(reads$adriftReadSeq, reads$adriftLinkerSeqEnd+1, nchar(reads$adriftReadSeq))

  reads <- lazy_dt(reads, immutable = FALSE) %>% 
           dplyr::group_by(uniqueSample, adriftReadRandomID, anchorReadSeq, adriftReadSeq2) %>%
           dplyr::mutate(nDuplicateReads = n() - 1) %>% 
           dplyr::slice(1) %>% 
           dplyr::ungroup() %>%
           dplyr::select(-adriftReadSeq2) %>%
           as.data.table()
  
} else if(opt$demultiplex_level == 'clustered'){
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


# Save demultiplexed reads and clean up.
reads$seqRunID <- opt$demultiplex_seqRunID
saveRDS(reads, file =  file.path(opt$outputDir, opt$demultiplex_outputDir, 'reads.rds'), compress = opt$compressDataFiles)
invisible(file.remove(list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp'), full.names = TRUE)))


# FASTQ export if requested in the configuration file.
if(opt$demultiplex_exportFASTQ){
  
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

q(save = 'no', status = 0, runLast = FALSE) 
