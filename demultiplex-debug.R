library(ShortRead)
library(readr)
library(yaml)
library(parallel)
library(dplyr)

opt <- read_yaml('config.yml')

if(! file.exists(opt$sampleConfigFile)) stop('Error - sample configuration file could not be found.')
if(! file.exists(opt$adriftReadsFile)) stop('Error - R1 seq file could not be found.')
if(! file.exists(opt$anchorReadsFile)) stop('Error - R2 seq file could not be found.')
if(! file.exists(opt$index1ReadsFile)) stop('Error - I1 seq file could not be found.')

source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, 'demultiplex'))
dir.create(file.path(opt$outputDir, 'demultiplex', 'log'))
dir.create(file.path(opt$outputDir, 'demultiplex', 'tmp'))
dir.create(file.path(opt$outputDir, 'demultiplex', 'seqChunks'))


index1Reads <- shortRead2DNAstringSet(readFastq(opt$index1ReadsFile))

# Correct Golay encoded barcodes if requested.
if(opt$correctGolayIndexReads){
  cluster <- makeCluster(opt$correctGolayIndexReads_CPUs)
  clusterExport(cluster, c('opt'))
  index1Reads <- Reduce('append', parLapply(cluster, split(index1Reads, ntile(1:length(index1Reads), opt$correctGolayIndexReads_CPUs)), golayCorrection))
  stopCluster(cluster)
}

samples <- read_tsv(opt$sampleConfigFile, col_types=cols())
if(nrow(samples) == 0) stop('Error - no lines of information was read from the sample configuration file.')


opt$outputDir <- file.path(opt$outputDir, 'demultiplex')

if(! 'subject' %in% names(samples))   samples$subject <- 'subject'
if(! 'replicate' %in% names(samples)) samples$replicate <- 1
if(any(grepl('~|\\|', paste(samples$subject, samples$sample, samples$replicate)))) stop('Error -- tildas (~) are reserved characters and can not be used in the subject, sample, or replicate sample configuration columns.')


# Create unique sample identifiers -- subject~sample~replicate
samples$uniqueSample <- paste0(samples$subject, '~', samples$sample, '~', samples$replicate)
if(any(duplicated(samples$uniqueSample))) stop('Error -- all subject, sample, replicate id combinations are not unique.')


cluster <- makeCluster(opt$demultiplexing_CPUs)
clusterExport(cluster, c('opt', 'samples'))


# Quality trim virus reads and break reads.
invisible(parLapply(cluster,     
                    list(c(opt$adriftReadsFile,  opt$demultiplexing_sequenceChunkSize, 'adriftReads',  file.path(opt$outputDir, 'seqChunks')),
                         c(opt$anchorReadsFile,  opt$demultiplexing_sequenceChunkSize, 'anchorReads',  file.path(opt$outputDir, 'seqChunks'))), 
                    function(x){
                      library(ShortRead)
                      source(file.path(opt$softwareDir, 'lib.R'))
                      qualTrimReads(x[[1]], x[[2]], x[[3]], x[[4]])
                    }))

adriftReads <- Reduce('append', lapply(list.files(file.path(opt$outputDir, 'seqChunks'), pattern = 'adriftReads', full.names = TRUE), function(x) shortRead2DNAstringSet(readFastq(x))))
anchorReads <- Reduce('append', lapply(list.files(file.path(opt$outputDir, 'seqChunks'), pattern = 'anchorReads', full.names = TRUE), function(x) shortRead2DNAstringSet(readFastq(x))))

reads <- syncReads(shortRead2DNAstringSet(readFastq(opt$index1ReadsFile)), anchorReads, adriftReads)
index1Reads <- reads[[1]];  anchorReads  <- reads[[2]];  adriftReads  <- reads[[3]]

invisible(file.remove(list.files(file.path(opt$outputDir, 'seqChunks'), full.names = TRUE)))
rm(reads)
gc()


#------------------------------------------
ooo <- readLines('~/chr22-1840213.reads')
anchorReads <- anchorReads[names(anchorReads) %in% ooo]
opt$demultiplexing_CPUs <- 1
reads <- syncReads(index1Reads, anchorReads, adriftReads)
index1Reads <- reads[[1]];  anchorReads  <- reads[[2]];  adriftReads  <- reads[[3]]
samples <- subset(samples, sample == 'GTSP2168')



# Split the trimmed reads into chunks for parallel processing.
chunkNum <- 1
d <- tibble(i = ntile(1:length(index1Reads), opt$demultiplexing_CPUs), n = 1:length(index1Reads))
invisible(lapply(split(d, d$i), function(x){
  index1Reads <- index1Reads[min(x$n):max(x$n)]
  anchorReads  <- anchorReads[min(x$n):max(x$n)]
  adriftReads  <- adriftReads[min(x$n):max(x$n)]
  save(index1Reads, anchorReads, adriftReads, file = file.path(opt$outputDir, 'seqChunks', chunkNum))
  chunkNum <<- chunkNum + 1
}))


# Clean up and free up memory. 
rm(d, chunkNum, index1Reads, anchorReads, adriftReads)
gc()

#-----------------------------------------



#invisible(parLapply(cluster, list.files(file.path(opt$outputDir, 'seqChunks'), full.names = TRUE), function(f){
invisible(lapply(list.files(file.path(opt$outputDir, 'seqChunks'), full.names = TRUE), function(f){
  library(ShortRead)
  library(dplyr)
  library(stringr)
  source(file.path(opt$softwareDir, 'lib.R'))
  
  load(f)
  
  browser()
  
  # Capture the chunk identifier.
  chunk.n <- unlist(str_match_all(f, '(\\d+)$'))[2]
  
  
  # Loop through samples in sample data file to demultiplex and apply read specific filters.
  invisible(lapply(1:nrow(samples), function(r){
    r <- samples[r,]
    
    # Create barcode demultiplexing vectors.
    v1 <- vcountPattern(r$index1.seq, index1Reads, max.mismatch = opt$demultiplexing_index1ReadFileMaxMismatch) > 0
    
    log.report <- tibble(sample = r$uniqueSample, demultiplexedIndex1Reads = sum(v1))
    
    # Create break read linker barcode demultiplexing vector.
    v2 <- rep(TRUE, length(adriftReads))
    if(opt$demultiplexing_useAdriftReadUniqueLinkers){
      testSeq <- substr(r$adriftRead.linker.seq, r$adriftRead.linkerBarcode.start, r$adriftRead.linkerBarcode.end)
      v2 <- vcountPattern(testSeq, subseq(adriftReads, r$adriftRead.linkerBarcode.start, r$adriftRead.linkerBarcode.end), max.mismatch = opt$demultiplexing_adriftReadLinkerBarcodeMaxMismatch) > 0
      log.report$demultiplexedLinkerReads <- sum(v2)
    } else {
      log.report$demultiplexedLinkerReads <- NA
    }
    
    
    # Test to see if any reads demultiplex to this row of the sample table and then subset reads to this sample.
    i <- base::intersect(which(v1), which(v2))
    if(length(i) == 0){
      log.report$demultiplexedReads <- 0
    } else {
      reads <- syncReads(index1Reads[i], anchorReads[i], adriftReads[i])
      index1Reads <- reads[[1]];  anchorReads  <- reads[[2]]; adriftReads  <- reads[[3]]
      
      if(length(index1Reads) == 0){
          log.report$demultiplexedReads <- 0
      } else {
        writeFasta(anchorReads, file.path(opt$outputDir, 'tmp', paste0(r$uniqueSample, '.', chunk.n, '.anchorReads')))
        writeFasta(adriftReads, file.path(opt$outputDir, 'tmp', paste0(r$uniqueSample, '.', chunk.n, '.adriftReads')))
        log.report$demultiplexedReads <- length(index1Reads)
      }
    }
  
    write.table(log.report, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = file.path(opt$outputDir, 'log', paste0(r$uniqueSample, '.', chunk.n, '.logReport')))
  }))
}))

invisible(unlink(file.path(opt$outputDir, 'seqChunks'), recursive = TRUE))
    

# Colate chunked reads and write out sample read files.
invisible(lapply(unique(samples$uniqueSample), function(x){
  f1 <- list.files(file.path(opt$outputDir, 'tmp'), pattern = paste0(x, '\\.\\d+\\.anchorReads'), full.names = TRUE)
  f2 <- list.files(file.path(opt$outputDir, 'tmp'), pattern = paste0(x, '\\.\\d+\\.adriftReads'), full.names = TRUE)
  if(length(f1) == 0 | length(f2) == 0 | length(f1) != length(f2)) return()
  
  anchorReads <- Reduce('append', lapply(f1, readFasta))
  adriftReads <- Reduce('append', lapply(f2, readFasta))
  
  writeFasta(anchorReads, file.path(opt$outputDir, paste0(x, '.anchorReads.fasta')))
  writeFasta(adriftReads, file.path(opt$outputDir, paste0(x, '.adriftReads.fasta')))
}))

invisible(unlink(file.path(opt$outputDir, 'tmp'), recursive = TRUE))

# Collect all the logs from the different computational nodes and create a single report.
logReport <- bind_rows(lapply(list.files(file.path(opt$outputDir, 'log'), pattern = '*.logReport$', full.names = TRUE), function(f){
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
}))

invisible(unlink(file.path(opt$outputDir, 'log'), recursive = TRUE))
write.table(logReport, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = file.path(opt$outputDir, 'log'))


