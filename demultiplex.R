library(ShortRead)
library(readr)
library(parallel)
library(lubridate)
library(dplyr)
library(data.table)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

# The launch script creates (if missing) and writes to the output directory.
# File write permission issues should be caught before starting modules.
dir.create(file.path(opt$outputDir, opt$demultiplex_outputDir))
dir.create(file.path(opt$outputDir, opt$demultiplex_outputDir, 'seqChunks'))
dir.create(file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'))

if(! file.exists(opt$sampleConfigFile)){
  write(c(paste(now(), '   Error - the sample configuration file could not be found')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}

write(c(paste(now(), '   Loading sample data')), file = file.path(opt$outputDir, 'log'), append = TRUE)
samples <- loadSamples()

if(! file.exists(opt$demultiplex_adriftReadsFile)){
  write(c(paste(now(), 'Error - the adrift reads file could not be found')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}

if(! file.exists(opt$demultiplex_anchorReadsFile)){
    write(c(paste(now(), '   Error - the index reads file could not be found')), file = file.path(opt$outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE) 
}

if(! file.exists(opt$demultiplex_index1ReadsFile)){
  write(c(paste(now(), '   Error - the anchor reads file could not be found')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
} 

if(opt$demultiplex_RC_I1_barcodes_auto){
  write(c(paste(now(), '   Determining if I1 barcodes should be switched to RC.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  r <- determine_RC_I1()
  write(paste(now(), '   Setting demultiplex_RC_I1_barcodes to', r), file = file.path(opt$outputDir, 'log'), append = TRUE)
  opt$demultiplex_RC_I1_barcodes <- r
}

# Reverse compliment index1 sequences if requested.
if(opt$demultiplex_RC_I1_barcodes) samples$index1.seq <- as.character(reverseComplement(DNAStringSet(samples$index1.seq)))

write(c(paste(now(), '   Reading in index 1 sequencing data.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
index1Reads <- shortRead2DNAstringSet(readFastq(opt$demultiplex_index1ReadsFile))

cluster <- makeCluster(opt$demultiplex_CPUs)
clusterExport(cluster, c('opt', 'samples'))

# Correct Golay encoded barcodes if requested.
if(opt$demultiplex_correctGolayIndexReads){
  write(paste(now(), '   Starting Golay correction.'), file = file.path(opt$outputDir, 'log'), append = TRUE)
  index1Reads.org <- index1Reads
  index1Reads <- Reduce('append', parLapply(cluster, split(index1Reads, ntile(1:length(index1Reads), opt$demultiplex_CPUs)), golayCorrection))
  percentChanged <- (sum(! as.character(index1Reads.org) == as.character(index1Reads)) / length(index1Reads))*100
  write(paste(now(), '   Golay correction complete.', sprintf("%.2f", percentChanged), '% of reads updated via Golay correction.'), file = file.path(opt$outputDir, 'log'), append = TRUE)
  rm(index1Reads.org)
}


# Quality trim virus reads and break reads.
write(paste(now(), '   Trimming anchor and adrift reads.'), file = file.path(opt$outputDir, 'log'), append = TRUE)
invisible(parLapply(cluster,     
                    list(c(opt$demultiplex_adriftReadsFile,  opt$demultiplex_sequenceChunkSize, 'adriftReads',  file.path(opt$outputDir, opt$demultiplex_outputDir, 'seqChunks')),
                         c(opt$demultiplex_anchorReadsFile,  opt$demultiplex_sequenceChunkSize, 'anchorReads',  file.path(opt$outputDir, opt$demultiplex_outputDir, 'seqChunks'))), 
                    function(x){
                      library(ShortRead)
                      source(file.path(opt$softwareDir, 'lib.R'))
                      qualTrimReads(x[[1]], x[[2]], x[[3]], x[[4]])
                    }))

write(paste(now(), '   Syncing anchor and adrift reads post-trimming.'), file = file.path(opt$outputDir, 'log'), append = TRUE)
adriftReads <- Reduce('append', lapply(list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'seqChunks'), pattern = 'adriftReads', full.names = TRUE), function(x) shortRead2DNAstringSet(readFastq(x))))
anchorReads <- Reduce('append', lapply(list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'seqChunks'), pattern = 'anchorReads', full.names = TRUE), function(x) shortRead2DNAstringSet(readFastq(x))))

reads <- syncReads(shortRead2DNAstringSet(readFastq(opt$demultiplex_index1ReadsFile)), anchorReads, adriftReads)
index1Reads <- reads[[1]];  anchorReads  <- reads[[2]];  adriftReads  <- reads[[3]]

invisible(file.remove(list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'seqChunks'), full.names = TRUE)))
rm(reads)
gc()


# Split the trimmed reads into chunks for parallel processing.
write(paste(now(), '   Dividing sequencing data into chunks.'), file = file.path(opt$outputDir, 'log'), append = TRUE)
chunkNum <- 1
d <- tibble(i = ntile(1:length(index1Reads), opt$demultiplex_CPUs), n = 1:length(index1Reads))
invisible(lapply(split(d, d$i), function(x){
  index1Reads <- index1Reads[min(x$n):max(x$n)]
  anchorReads <- anchorReads[min(x$n):max(x$n)]
  adriftReads <- adriftReads[min(x$n):max(x$n)]
  save(index1Reads, anchorReads, adriftReads, file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'seqChunks', chunkNum))
  chunkNum <<- chunkNum + 1
}))


# Clean up and free up memory. 
rm(d, chunkNum, index1Reads, anchorReads, adriftReads)
gc()


write(paste(now(), '   Demultiplexing sequence chunks.'), file = file.path(opt$outputDir, 'log'), append = TRUE)
if('anchorRead.startSeq' %in% names(samples)) write(paste(now(), '   Anchor read start sequence filter enabled.'), file = file.path(opt$outputDir, 'log'), append = TRUE)


invisible(parLapply(cluster, list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'seqChunks'), full.names = TRUE), function(f){
#invisible(lapply(list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'seqChunks'), full.names = TRUE), function(f){
  library(ShortRead)
  library(dplyr)
  library(stringr)
  source(file.path(opt$softwareDir, 'lib.R'))
  
  load(f)
 
  # Capture the chunk identifier.
  chunk.n <- unlist(str_match_all(f, '(\\d+)$'))[2]
  
  # Loop through samples in sample data file to demultiplex and apply read specific filters.
  invisible(lapply(1:nrow(samples), function(r){
    r <- samples[r,]

    v0 <- rep(TRUE, length(anchorReads))
    if('anchorRead.startSeq' %in% names(r)){
      v0 <- vcountPattern(r$anchorRead.startSeq, subseq(anchorReads, 1, nchar(r$anchorRead.startSeq)), max.mismatch = opt$demultiplex_anchorRead.startSeq.maxMisMatch) == 1
    }
       
    # Create barcode demultiplexing vectors.
    v1 <- vcountPattern(r$index1.seq, index1Reads, max.mismatch = opt$demultiplex_index1ReadMaxMismatch) > 0
    
    log.report <- tibble(sample = r$uniqueSample, demultiplexedIndex1Reads = sum(v1))
    
    # Create break read linker barcode demultiplexing vector.
    v2 <- rep(TRUE, length(adriftReads))
    if(opt$demultiplex_useAdriftReadUniqueLinkers){
      #browser()
      testSeq <- substr(r$adriftRead.linker.seq, r$adriftRead.linkerBarcode.start, r$adriftRead.linkerBarcode.end)
      v2 <- vcountPattern(testSeq, subseq(adriftReads, r$adriftRead.linkerBarcode.start, r$adriftRead.linkerBarcode.end), max.mismatch = opt$demultiplex_adriftReadLinkerBarcodeMaxMismatch) > 0
      log.report$demultiplexedLinkerReads <- sum(v2)
    } else {
      log.report$demultiplexedLinkerReads <- NA
    }
    
    # Test to see if any reads demultiplex to this row of the sample table and then subset reads to this sample.
    i <- Reduce(base::intersect, list(which(v0), which(v1), which(v2)))
    if(length(i) == 0){
      log.report$demultiplexedReads <- 0
    } else {
      reads <- syncReads(index1Reads[i], anchorReads[i], adriftReads[i])
      index1Reads <- reads[[1]];  anchorReads  <- reads[[2]]; adriftReads  <- reads[[3]]
      
      if(length(index1Reads) == 0){
          log.report$demultiplexedReads <- 0
      } else {
          writeFasta(subseq(adriftReads, r$adriftRead.linkerRandomID.start, r$adriftRead.linkerRandomID.end), 
                     file.path(opt$outputDir, 'tmp', paste0(r$uniqueSample, '.', chunk.n, '.randomAdriftReadIDs.gz')), compress = TRUE)
        
        writeFasta(anchorReads, file.path(opt$outputDir, 'tmp', paste0(r$uniqueSample, '.', chunk.n, '.anchorReads.gz')), compress = TRUE)
        writeFasta(adriftReads, file.path(opt$outputDir, 'tmp', paste0(r$uniqueSample, '.', chunk.n, '.adriftReads.gz')), compress = TRUE)
        log.report$demultiplexedReads <- length(index1Reads)
      }
    }
  
    write.table(log.report, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'log', paste0(r$uniqueSample, '.', chunk.n, '.logReport')))
  }))
}))

stopCluster(cluster)

invisible(unlink(file.path(opt$outputDir, opt$demultiplex_outputDir, 'seqChunks'), recursive = TRUE))

# Collate chunked reads and write out sample read files.
write(paste(now(), '   Colating data files.'), file = file.path(opt$outputDir, 'log'), append = TRUE)

reads <-  rbindlist(lapply(unique(samples$uniqueSample), function(x){
  f1 <- list.files(file.path(opt$outputDir, 'tmp'), pattern = paste0(x, '\\.\\d+\\.anchorReads'), full.names = TRUE)
  f2 <- list.files(file.path(opt$outputDir, 'tmp'), pattern = paste0(x, '\\.\\d+\\.adriftReads'), full.names = TRUE)
  f3 <- list.files(file.path(opt$outputDir, 'tmp'), pattern = paste0(x, '\\.\\d+\\.randomAdriftReadIDs'), full.names = TRUE)
  if(length(f1) == 0 | length(f2) == 0 | length(f1) != length(f2)) return()

  write(paste0(now(), '    Colating reads for ', x), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
  anchorReads <- Reduce('append', lapply(f1, readDNAStringSet))
  adriftReads <- Reduce('append', lapply(f2, readDNAStringSet))
  randomIDs   <- Reduce('append', lapply(f3, readDNAStringSet))
  
  closeAllConnections()
  
  r <- subset(samples, uniqueSample == x)
  c <- substr(r$adriftRead.linker.seq, max(stringr::str_locate_all(r$adriftRead.linker.seq, 'NNN')[[1]][,2])+1, nchar(r$adriftRead.linker.seq))
  t <- as.character(reverseComplement(DNAString(substr(c, nchar(c) - 14, nchar(c)))))

  data.table(uniqueSample = x, readID = names(anchorReads), anchorReadSeq = as.character(anchorReads), adriftReadSeq = as.character(adriftReads), 
             adriftReadRandomID = as.character(randomIDs), adriftReadTrimSeq = t, adriftLinkerSeqEnd = nchar(r$adriftRead.linker.seq),
             vectorFastaFile = r$vectorFastaFile, refGenome = r$refGenome.id, flags = r$flags)
}))

if(nrow(reads) == 0){
  write(c(paste(now(), '   Error - no reads were demultiplexed for any sample.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}

write(paste(now(), '   Clearing tmp files.'), file = file.path(opt$outputDir, 'log'), append = TRUE)
invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

# Collect all the logs from the different computational nodes and create a single report.
write(paste(now(), '   Colating log files.'), file = file.path(opt$outputDir, 'log'), append = TRUE)

logReport <- bind_rows(lapply(list.files(file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), pattern = '*.logReport$', full.names = TRUE), function(f){
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


invisible(unlink(file.path(opt$outputDir, opt$demultiplex_outputDir, 'log'), recursive = TRUE))
write.table(logReport, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'readAttritionTbl.tsv'))

write(paste(now(), '   Writing outputs.'), file = file.path(opt$outputDir, 'log'), append = TRUE)

saveRDS(reads, file =  file.path(opt$outputDir, opt$demultiplex_outputDir, 'reads.rds'), compress = FALSE)
readr::write_csv(reads,  file.path(opt$outputDir, opt$demultiplex_outputDir, 'reads.csv.gz'))

q(save = 'no', status = 0, runLast = FALSE) 
