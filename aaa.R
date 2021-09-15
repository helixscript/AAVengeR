library(ShortRead)
library(dplyr)
library(parallel)
library(readr)
options(stringsAsFactors = FALSE)

opt <- yaml::read_yaml('config.yml')

source(file.path(opt$softwareDir, 'lib.R'))

opt$inputFastaDir <- file.path(opt$outputDir, 'uniqueFasta')
opt$outputDir <- file.path(opt$outputDir, 'leaderSeqMappings')
dir.create(opt$outputDir)
dir.create(file.path(opt$outputDir, 'tmp'))
dir.create(file.path(opt$outputDir, 'dbs'))

samples <- read_tsv(opt$sampleConfigFile, col_types=cols())
samples$uniqueSample <- paste0(samples$subject, '~', samples$sample, '~', samples$replicate)

cluster <- makeCluster(opt$alignReads_CPUs)
clusterExport(cluster, c('opt', 'samples'))

anchorReads <- shortRead2DNAstringSet(Reduce('append', lapply(f[grepl('anchorReads', f)], readFasta)))

CPUs <- 20

m <- bind_rows(lapply(list.files(opt$inputFastaDir, full.names = TRUE, pattern = 'anchor'), function(x){
       library(Biostrings)
       library(dplyr)
       r <- shortRead2DNAstringSet(readFasta(x))
       o <- unlist(strsplit(x, '/'))
       u <- unlist(strsplit(o[length(o)], '\\.'))[1]
  
       o <- strsplit(unlist(strsplit(subset(samples, uniqueSample == u)$anchorRead.identification, ';')), ',')
       d <- DNAStringSet(unlist(lapply(o, '[', 2)))
       names(d) <- unlist(lapply(o, '[', 1))
       writeXStringSet(d, file.path(opt$outputDir, 'dbs', paste0(u, '.fasta')))
       writeXStringSet(r, file.path(opt$outputDir, 'tmp', paste0(u, '.fasta')))
  
       system(paste0(opt$command_makeblastdb, ' -in ', file.path(opt$outputDir, 'dbs', paste0(u, '.fasta')), 
                    ' -dbtype nucl -out ', file.path(opt$outputDir, 'dbs', u)), ignore.stderr = TRUE)
  
       waitForFile(file.path(opt$outputDir, 'dbs', paste0(u, '.nin')))
  
       system(paste0(opt$command_blastn, ' -word_size 4 -evalue 50 -outfmt 6 -query ',
                     file.path(opt$outputDir, 'tmp', paste0(u, '.fasta')), ' -db ',
                     file.path(opt$outputDir, 'dbs', u),
                     ' -num_threads 5 -out ', file.path(opt$outputDir, 'tmp', paste0(u, '.blast'))),
             ignore.stdout = TRUE, ignore.stderr = TRUE)
  
      waitForFile(file.path(opt$outputDir, 'tmp', paste0(u, '.blast')))
  
      if(file.info(file.path(opt$outputDir, 'tmp', paste0(u, '.blast')))$size == 0) return(tibble())
  
      b <- read.table(file.path(opt$outputDir, 'tmp', paste0(u, '.blast')), sep = '\t', header = FALSE)
      names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
      b$alignmentLength <- b$qend - b$qstart + 1
  
      b <- dplyr::filter(b, pident >= opt$mapLeaderSequences_minAlignmentPercentID, alignmentLength >= opt$mapLeaderSequences_minAlignmentLength)
  
      bind_rows(lapply(split(b, b$qname), function(x){
       x <- subset(x, evalue == min(x$evalue)) 
       x <- subset(x, alignmentLength == max(x$alignmentLength))
    
       tibble(id = x$qname[1], leaderMaping.id = paste(unique(x$sseqid), collapse = '|'), 
              leaderMapping.qStart = min(x$qstart), leaderMapping.qEnd = max(x$qend),
              leaderMapping.sStart = min(x$sstart), leaderMapping.sEnd = max(x$send))
     }))
}))
  