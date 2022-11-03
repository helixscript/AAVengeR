library(ShortRead)
library(dplyr)
library(stringr)
library(parallel)
options(stringsAsFactors = FALSE)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)

source(file.path(opt$softwareDir, 'lib.R'))

invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

sites <- readRDS(file.path(opt$outputDir, opt$mapSiteLeaderSequences_inputFile))
sites$repLeaderSeqLength <- NULL
sites$repLeaderSeqMap <- NULL

dir.create(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir))
dir.create(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'dbs'))

sites$s <- paste(sites$subject, sites$sample)

cluster <- makeCluster(opt$mapSiteLeaderSequences_CPUs)

m <- bind_rows(lapply(split(sites, sites$vectorFastaFile), function(x){
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'dbs'), full.names = TRUE)))
  
  system(paste0(opt$command_makeblastdb, ' -in ', file.path(opt$softwareDir, 'data', 'vectors', x$vectorFastaFile[1]), 
                ' -dbtype nucl -out ', file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
  
  waitForFile(file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'dbs', 'd.nin'))
  
  s <- data.frame(seq = unique(x$repLeaderSeq))
  s$id <- paste0('s', 1:nrow(s))
  
  s$cut <- cut(nchar(s$seq), c(-Inf, seq(0, max(nchar(s$seq)), by = 10), Inf), labels = FALSE)
  s <- group_by(s, cut) %>% mutate(n = ntile(1:n(), opt$mapSiteLeaderSequences_CPUs)) %>% ungroup()
  
  reads <- DNAStringSet(s$seq)
  names(reads) <- s$id
  
  b <- bind_rows(lapply(split(reads, s$n), function(a){
         f <- tmpFile()
 
         writeXStringSet(a,  file.path(opt$outputDir, 'tmp', paste0(f, '.fasta')))
    
         system(paste0(opt$command_blastn, ' -word_size 4 -evalue 50 -outfmt 6 -query ',
                  file.path(opt$outputDir, 'tmp', paste0(f, '.fasta')), ' -db ',
                  file.path(opt$outputDir, opt$mapSiteLeaderSequences_outputDir, 'dbs', 'd'),
                  ' -num_threads 1 -out ', file.path(opt$outputDir, 'tmp', paste0(f, '.blast'))),
                ignore.stdout = TRUE, ignore.stderr = TRUE)
    
         waitForFile(file.path(opt$outputDir, 'tmp', paste0(f, '.blast')))
    
         if(file.info(file.path(opt$outputDir, 'tmp', paste0(f, '.blast')))$size == 0) return(tibble())
    
          b <- read.table(file.path(opt$outputDir, 'tmp', paste0(f, '.blast')), sep = '\t', header = FALSE)
          names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
          b$alignmentLength <- b$qend - b$qstart + 1
          b$strand <- ifelse(b$sstart > b$send, '-', '+')
    
          dplyr::filter(b, pident >= opt$mapSiteLeaderSequences_minAlignmentPercentID, alignmentLength >= opt$mapSiteLeaderSequences_minAlignmentLength)
       }))
  
  # New
  b$i <- group_by(b, qname) %>% group_indices() 
  o <- tibble(i2 = 1:n_distinct(b$i))
  o$n <- ntile(1:nrow(o), opt$prepReads_CPUs)
  b <- left_join(b, o, by = c('i' = 'i2'))
  r <- rbindlist(parallel::parLapply(cluster, split(b, b$n), blast2rearangements_worker))
  
  left_join(tibble(qname = names(reads), leaderSeq = as.character(reads)), dplyr::rename(r, repLeaderSeqMap = rearrangement), by = 'qname')
}))

stopCluster(cluster)

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
