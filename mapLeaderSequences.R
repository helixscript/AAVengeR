library(dplyr)
library(yaml)
library(parallel)
library(Biostrings)

opt <- read_yaml('config.yml')
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, 'mapLeaderSequences'))
dir.create(file.path(opt$outputDir, 'mapLeaderSequences', 'dbs'))
dir.create(file.path(opt$outputDir, 'mapLeaderSequences', 'tmp'))

frags <- readRDS(file.path(opt$outputDir, 'standardizeFragments', 'frags.rds'))

samples <- readr::read_tsv(opt$sampleConfigFile)
samples$uniqueSample <- paste0(samples$subject, '~', samples$sample, '~', samples$replicate)

samples <- subset(samples, uniqueSample %in% frags$uniqueSample)

frags <- left_join(frags, select(samples, uniqueSample, anchorRead.identification), by = 'uniqueSample')  

cluster <- makeCluster(opt$mapLeaderSequences_CPUs)
clusterExport(cluster, c('opt', 'samples'))

frags <- bind_rows(lapply(split(frags, frags$anchorRead.identification), function(x){
  o <- strsplit(unlist(strsplit(x$anchorRead.identification[1], ';')), ',')
  d <- DNAStringSet(unlist(lapply(o, '[', 2)))
  names(d) <- unlist(lapply(o, '[', 1))
  writeXStringSet(d, file.path(opt$outputDir, 'mapLeaderSequences', 'dbs', 'd.fasta'))
  
  system(paste0(opt$command_makeblastdb, ' -in ', file.path(opt$outputDir, 'mapLeaderSequences', 'dbs', 'd.fasta'), 
                ' -dbtype nucl -out ', file.path(opt$outputDir, 'mapLeaderSequences', 'dbs', 'd')), ignore.stderr = TRUE)
  waitForFile(file.path(opt$outputDir, 'mapLeaderSequences', 'dbs', 'd.nin'))
  
  d <- tibble(id = paste0('s', 1:length(unique(x$repLeaderSeq))), seq = unique(x$repLeaderSeq))
  
  d$s <- dplyr::ntile(1:nrow(d), opt$mapLeaderSequences_CPUs)
  
  invisible(parLapply(cluster, split(d, d$s), function(a){
  #invisible(lapply(split(d, d$s), function(a){  
    library(Biostrings)
    o <- DNAStringSet(a$seq)
    names(o) <- a$id
    writeXStringSet(o, file.path(opt$outputDir, 'mapLeaderSequences', 'tmp', paste0('s', a$s[1])))
    
    system(paste0(opt$command_blastn, ' -word_size 5 -evalue 50 -outfmt 6 -query ',
                  file.path(opt$outputDir, 'mapLeaderSequences', 'tmp', paste0('s', a$s[1])), ' -db ',
                  file.path(opt$outputDir, 'mapLeaderSequences', 'dbs', 'd'),
                  ' -num_threads 5 -out ', file.path(opt$outputDir, 'mapLeaderSequences', 'tmp', paste0('s', a$s[1], '.blast'))),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
  }))
  
  f <- tibble(file = list.files(file.path(opt$outputDir, 'mapLeaderSequences', 'tmp'), pattern = '.blast$', full.names = TRUE))
  f$fileSize <- sapply(f$file, function(x) file.info(x)$size)
  f <- subset(f, fileSize > 0)

  b <- bind_rows(lapply(f$file, function(x){  
         o <- read.table(x, sep = '\t', header = FALSE)
         names(o) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
         o$alignmentLength <- o$qend - o$qstart + 1
         o
       })) %>% dplyr::filter(pident >= opt$mapLeaderSequences_minAlignmentPercentID,
                             alignmentLength >= opt$mapLeaderSequences_minAlignmentLength)
  
  b2 <- bind_rows(lapply(split(b, b$qname), function(x){
          x <- subset(x, evalue == min(x$evalue)) 
          x <- subset(x, alignmentLength == max(x$alignmentLength))
          tibble(id = x$qname[1], leaderMaping.id = paste(unique(x$sseqid), collapse = '|'), 
                 leaderMapping.qStart = min(x$qstart), leaderMapping.qEnd = max(x$qend),
                 leaderMapping.sStart = min(x$sstart), leaderMapping.sEnd = max(x$send))
        }))
  
  d <- left_join(d, b2, by = 'id') %>% dplyr::select(-id, -s)
  left_join(x, d, by = c('repLeaderSeq' = 'seq')) %>% dplyr::select(-anchorRead.identification)
}))

unlink(file.path(opt$outputDir, 'mapLeaderSequences', 'tmp'), recursive = TRUE)
unlink(file.path(opt$outputDir, 'mapLeaderSequences', 'dbs'), recursive = TRUE)

saveRDS(frags, file.path(opt$outputDir, 'mapLeaderSequences', 'frags.rds'))

