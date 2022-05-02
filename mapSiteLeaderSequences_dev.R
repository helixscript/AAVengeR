library(ShortRead)
library(dplyr)
library(stringr)
library(parallel)
options(stringsAsFactors = FALSE)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)

# Min. space between blast hits allowed. Larger gaps will be filled with [x].
opt$mapSiteLeaderSequences_minAllowableGap <- 5

source(file.path(opt$softwareDir, 'lib.R'))

invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

sites <- readRDS(file.path(opt$outputDir, opt$mapSiteLeaderSequences_inputFile))
sites$repLeaderSeqLength <- NULL
sites$repLeaderSeqMap <- NULL

dir.create(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir))
dir.create(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'dbs'))

samples <- loadSamples()

sites$s <- paste(sites$subject, sites$sample)

m <- bind_rows(lapply(split(samples, samples$vectorFastaFile), function(x){
  
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'dbs'), full.names = TRUE)))
  
  system(paste0(opt$command_makeblastdb, ' -in ', file.path(opt$softwareDir, 'data', 'vectors', x$vectorFastaFile[1]), 
                ' -dbtype nucl -out ', file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
  
  #waitForFile(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'dbs', 'd.nin'))
  
  browser()
  
  # Subset sites to match this identification sequence string.
  #s <- subset(sites, s %in% unique(paste(x$subject, x$sample)))
  s <- sites
  
  reads <- DNAStringSet(unique(s$repLeaderSeq))
  names(reads) <- paste0('s', 1:length(reads))
  
  
  
  b <- bind_rows(lapply(mixAndChunkSeqs(reads, opt$mapSiteLeaderSequences_alignmentChunkSize), function(a){
    f <- tmpFile()
 
    writeXStringSet(a,  file.path(opt$outputDir, 'tmp', paste0(f, '.fasta')))
    
    system(paste0(opt$command_blastn, ' -word_size 4 -evalue 50 -outfmt 6 -query ',
                  file.path(opt$outputDir, 'tmp', paste0(f, '.fasta')), ' -db ',
                 #' ~/test ', ' -db ',
                  file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'dbs', 'd'),
                  ' -num_threads 4 -out ', file.path(opt$outputDir, 'tmp', paste0(f, '.blast'))),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    #waitForFile(file.path(opt$outputDir, 'tmp', paste0(f, '.blast')))
    browser()
    
    if(file.info(file.path(opt$outputDir, 'tmp', paste0(f, '.blast')))$size == 0) return(tibble())
    
    b <- read.table(file.path(opt$outputDir, 'tmp', paste0(f, '.blast')), sep = '\t', header = FALSE)
    names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
    b$alignmentLength <- b$qend - b$qstart + 1
    
    dplyr::filter(b, pident >= opt$mapSiteLeaderSequences_minAlignmentPercentID, alignmentLength >= opt$mapSiteLeaderSequences_minAlignmentLength)
  }))
  
  r <- blast2rearangements(b, minAlignmentLength = opt$mapSiteLeaderSequences_minAlignmentLength, 
                           minPercentID = opt$mapSiteLeaderSequences_minAlignmentPercentID, CPUs = opt$mapSiteLeaderSequences_CPUs)
  
  left_join(tibble(qname = names(reads), leaderSeq = as.character(reads)), dplyr::rename(r, repLeaderSeqMap = rearrangement), by = 'qname')
}))

invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

sites <- select(sites, -s) %>% left_join(select(m, leaderSeq, repLeaderSeqMap), by = c('repLeaderSeq' = 'leaderSeq'))
sites$repLeaderSeqLength <- nchar(sites$repLeaderSeq)

# Fill in missing pieces.
sites$repLeaderSeqMap <- unlist(lapply(split(sites, 1:nrow(sites)), function(x){
  if(is.na(x$repLeaderSeqMap)) return(NA)
  o <- unlist(strsplit(as.character(x$repLeaderSeqMap), ';'))
  lastRangeEnd <- as.integer(sub('\\.\\.', '', str_extract(o[length(o)], '..(\\d+)')))
  
  if((x$repLeaderSeqLength - lastRangeEnd) > opt$mapSiteLeaderSequences_minAllowableGap){
    x$repLeaderSeqMap <- paste0(x$repLeaderSeqMap, ';', lastRangeEnd+1, '..', x$repLeaderSeqLength, '[x]')
  }
  
  o <- unlist(strsplit(as.character(x$repLeaderSeqMap), ';'))
  if(length(o) > 1){
    invisible(lapply(1:(length(o)-1), function(x){
      lastRangeEnd <- as.integer(sub('\\.\\.', '', str_extract(o[[x]], '..(\\d+)')))
      nextFirstRangeEnd <- as.integer(sub('\\.\\.', '', str_extract(o[[x+1]], '(\\d+)..')))
    
      if((nextFirstRangeEnd - lastRangeEnd) > opt$mapSiteLeaderSequences_minAllowableGap){
        o[[x]] <<- paste0(o[[x]], ';', lastRangeEnd+1, '..', nextFirstRangeEnd-1, '[x];')
      }
    }))
  }
  
  gsub(';;', ';', paste0(o, collapse = ';'))
}))

saveRDS(sites, file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, opt$mapSiteLeaderSequences_outputFile))
