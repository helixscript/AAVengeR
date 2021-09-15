library(ShortRead)
library(dplyr)
library(parallel)
options(stringsAsFactors = FALSE)

opt <- yaml::read_yaml('config.yml')

source(file.path(opt$softwareDir, 'lib.R'))

opt$inputFastaDir <- file.path(opt$outputDir, 'uniqueFasta')
opt$outputDir <- file.path(opt$outputDir, 'vectorFilter')
dir.create(opt$outputDir)
dir.create(file.path(opt$outputDir, 'tmp'))
dir.create(file.path(opt$outputDir, 'dbs'))

samples <- readr::read_tsv(opt$sampleConfigFile, col_types = readr::cols())
if(! 'subject' %in% names(samples))   samples$subject <- 'subject'
if(! 'replicate' %in% names(samples)) samples$replicate <- 1
samples$uniqueSample <- paste0(samples$subject, '~', samples$sample, '~', samples$replicate)

cluster <- makeCluster(opt$vectorFilter_CPUs)
clusterExport(cluster, c('opt', 'samples'))

d <- tibble(file = list.files(opt$inputFastaDir, full.names = TRUE, pattern = 'anchor'))
d$uniqueSample <- unlist(lapply(d$file, function(x){ 
  o <- unlist(strsplit(x, '/'))
  unlist(strsplit(o[length(o)], '\\.'))[1]
})) 
d <- left_join(d, select(samples, uniqueSample, anchorRead.seqFilter.file), by = 'uniqueSample')

m <- bind_rows(lapply(split(d, d$anchorRead.seqFilter.file), function(x){
  invisible(file.remove(list.files(file.path(opt$outputDir, 'dbs'), full.names = TRUE)))
  
  system(paste0(opt$command_makeblastdb, ' -in ', x$anchorRead.seqFilter.file[1], ' -dbtype nucl -out ', file.path(opt$outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
  
  waitForFile(file.path(opt$outputDir, 'dbs', 'd.nin'))
  
  reads <- Reduce('append', lapply(x$file, readFasta))
  
  bind_rows(parLapply(cluster, split(reads, ntile(1:length(reads), opt$vectorFilter_CPUs)), function(a){
  #bind_rows(lapply(split(reads, ntile(1:length(reads), opt$vectorFilter_CPUs)), function(a){  
    library(Biostrings)
    library(dplyr)
    source(file.path(opt$softwareDir, 'lib.R'))
    
    f <- tmpFile()
    a.ids <- as.character(a@id)
    a <- a@sread
    names(a) <- a.ids 
    
    a <- subseq(a, (width(a) - (opt$vectorFilter_anchorReadSampleLength -1)), width(a))
    
    writeXStringSet(a,  file.path(opt$outputDir, 'tmp', paste0(f, '.fasta')))
    
    system(paste0(opt$command_blastn, ' -word_size 4 -evalue 50 -outfmt 6 -query ',
                  file.path(opt$outputDir, 'tmp', paste0(f, '.fasta')), ' -db ',
                  file.path(opt$outputDir, 'dbs', 'd'),
                  ' -num_threads 4 -out ', file.path(opt$outputDir, 'tmp', paste0(f, '.blast'))),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    waitForFile(file.path(opt$outputDir, 'tmp', paste0(f, '.blast')))
    
    if(file.info(file.path(opt$outputDir, 'tmp', paste0(f, '.blast')))$size == 0) return(tibble())
    
    b <- read.table(file.path(opt$outputDir, 'tmp', paste0(f, '.blast')), sep = '\t', header = FALSE)
    names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
    b$percentQueryAligned <- ((b$qend - b$qstart + 1) / opt$vectorFilter_anchorReadSampleLength) * 100
    
    dplyr::filter(b, pident >= opt$vectorFilter_percentCoverage, percentQueryAligned >= opt$vectorFilter_percentCoverage)
  }))
})) %>% tidyr::drop_na()

stopCluster(cluster)
invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

saveRDS(m, file.path(opt$outputDir, 'mappings.rds'))
