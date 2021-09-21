library(ShortRead)
library(dplyr)
library(parallel)
options(stringsAsFactors = FALSE)

opt <- yaml::read_yaml('config.yml')

source(file.path(opt$softwareDir, 'lib.R'))

sites <- readRDS(file.path(opt$outputDir, 'buildSites', 'sites.rds'))

dir.create(file.path(opt$outputDir, 'mapSiteLeaderSequences'))

dir.create(file.path(opt$outputDir, 'mapSiteLeaderSequences', 'tmp'))
dir.create(file.path(opt$outputDir, 'mapSiteLeaderSequences', 'dbs'))

samples <- loadSamples()

sites$s <- paste(sites$subject, sites$sample)

m <- bind_rows(lapply(split(samples, samples$anchorRead.identification), function(x){
  invisible(file.remove(list.files(file.path(opt$outputDir, 'mapSiteLeaderSequences', 'dbs'), full.names = TRUE)))
  
  o <- strsplit(unlist(strsplit(x$anchorRead.identification[1], ';')), ',')
  d <- DNAStringSet(unlist(lapply(o, '[', 2)))
  names(d) <- unlist(lapply(o, '[', 1))
  writeXStringSet(d, file.path(opt$outputDir, 'mapSiteLeaderSequences', 'dbs', paste0('d.fasta')))
  
  system(paste0(opt$command_makeblastdb, ' -in ', file.path(opt$outputDir, 'mapSiteLeaderSequences', 'dbs', 'd.fasta'), 
                ' -dbtype nucl -out ', file.path(opt$outputDir, 'mapSiteLeaderSequences', 'dbs', 'd')), ignore.stderr = TRUE)
  
  waitForFile(file.path(opt$outputDir, 'mapSiteLeaderSequences', 'dbs', 'd.nin'))
  
  # Subset sites to match this identification sequence string.
  s <- subset(sites, s %in% unique(paste(x$subject, x$sample)))
  
  reads <- DNAStringSet(unique(s$repLeaderSeq))
  names(reads) <- paste0('s', 1:length(reads))
  
  b <- bind_rows(lapply(mixAndChunkSeqs(reads, opt$mapLeaderSequences_alignmentChunkSize), function(a){
    f <- tmpFile()
 
    writeXStringSet(a,  file.path(opt$outputDir, 'mapSiteLeaderSequences', 'tmp', paste0(f, '.fasta')))
    
    system(paste0(opt$command_blastn, ' -word_size 4 -evalue 50 -outfmt 6 -query ',
                  file.path(opt$outputDir, 'mapSiteLeaderSequences', 'tmp', paste0(f, '.fasta')), ' -db ',
                  file.path(opt$outputDir, 'mapSiteLeaderSequences', 'dbs', 'd'),
                  ' -num_threads 4 -out ', file.path(opt$outputDir, 'mapSiteLeaderSequences', 'tmp', paste0(f, '.blast'))),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    waitForFile(file.path(opt$outputDir, 'mapSiteLeaderSequences', 'tmp', paste0(f, '.blast')))
    
    if(file.info(file.path(opt$outputDir, 'mapSiteLeaderSequences', 'tmp', paste0(f, '.blast')))$size == 0) return(tibble())
    
    b <- read.table(file.path(opt$outputDir, 'mapSiteLeaderSequences', 'tmp', paste0(f, '.blast')), sep = '\t', header = FALSE)
    names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
    b$alignmentLength <- b$qend - b$qstart + 1
    
    b <- dplyr::filter(b, pident >= opt$mapLeaderSequences_minAlignmentPercentID, alignmentLength >= opt$mapLeaderSequences_minAlignmentLength)

    bind_rows(lapply(split(b, b$qname), function(x){
      if('s4' %in% x$qname) browser()
      
      x <- subset(x, evalue == min(x$evalue)) 
      x <- subset(x, alignmentLength == max(x$alignmentLength))
      
      tibble(id = x$qname[1], leaderMaping.id = paste(unique(x$sseqid), collapse = '|'), 
             leaderMapping.qStart = min(x$qstart), leaderMapping.qEnd = max(x$qend),
             leaderMapping.sStart = min(x$sstart), leaderMapping.sEnd = max(x$send))
    }))
  }))
  
  r <- tibble(id = names(reads), leaderSeq = as.character(reads))
  left_join(r, b, by = 'id')
})) %>% tidyr::drop_na()

invisible(file.remove(list.files(file.path(opt$outputDir, 'mapSiteLeaderSequences', 'tmp'), full.names = TRUE)))

sites <- select(sites, -s) %>% left_join(select(m, -id), by = c('repLeaderSeq' = 'leaderSeq'))

saveRDS(sites, file.path(opt$outputDir, 'mapSiteLeaderSequences', 'sites.rds'))

