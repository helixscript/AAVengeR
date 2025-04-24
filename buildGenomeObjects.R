#!/usr/bin/Rscript

library(yaml)
library(dplyr)
library(Biostrings)
library(parallel)
library(GenomicRanges)
library(tidyr)
library(parallel)

tmpFile <- paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')

# Read in the configuration file and perform basic sanity checks.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) stop('Expected at least one command line argument')
opts <- yaml::read_yaml(args)

logFile <- paste0(opts$refGenomeName, '_buildGenomeObjects.log')
write(date(), logFile, append = FALSE)

dir.create(opts$outputDir, showWarnings = FALSE)
dir.create(file.path(opts$outputDir, 'genomeAnnotations'), showWarnings = FALSE)
dir.create(file.path(opts$outputDir, 'referenceGenomes'), showWarnings = FALSE)
dir.create(file.path(opts$outputDir, 'referenceGenomes', 'blat'), showWarnings = FALSE)
dir.create(file.path(opts$outputDir, 'referenceGenomes', 'bwa2'), showWarnings = FALSE)

# Download genome data from UCSC.
write(paste0(date(), ' Downloading genome data.'), logFile, append = TRUE)
system(paste0('wget -q -O ', opts$outputDir, '/', tmpFile, '.gz ', opts$refGenomeSeq_URL))
system(paste0('gunzip ', opts$outputDir, '/', tmpFile, '.gz'))

# Read in genome NTs and limit to standard chromosomes.
o <- readDNAStringSet(paste0(opts$outputDir, '/', tmpFile))
invisible(file.remove(paste0(opts$outputDir, '/', tmpFile)))

if(opts$forceStdChromosomes) o <- o[names(o) %in% c(paste0('chr', 1:100), paste0('chr', as.roman(1:100)), 'chrY', 'chrM')]

# # Run RepeatMasker.
# dir.create(paste0(opts$outputDir, '/repeatMasker'), showWarnings = FALSE)
# 
# n <- 1
# invisible(lapply(split(o, names(o)), function(x){
#   write(paste0(date(), ' Starting RepeatMasker for chromosome ', n, '/', n_distinct(names(o)), ' length: ', width(x)), logFile, append = TRUE)
#   n <<- n+1
#   tmpFile <- paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')
#   dir.create(paste0(opts$outputDir, '/repeatMasker/', tmpFile), showWarnings = FALSE)
#   writeXStringSet(x, paste0(opts$outputDir, '/repeatMasker/', tmpFile, '/seq.fasta'))
#   system(paste0(opts$repeatMaskerPath, ' -s -pa ', opts$CPUs, ' -e rmblast -species "', 
#                 opts$repeatMaskerSpecies, '" -dir ', paste0(opts$outputDir, '/repeatMasker/', tmpFile), ' ', 
#                 paste0(opts$outputDir, '/repeatMasker/', tmpFile, '/seq.fasta')))
#   }))
# 
# write(paste0(date(), ' Reading RepeatMasker outputs.'), logFile, append = TRUE)
# 
# cluster <- makeCluster(opts$CPUs)
# clusterExport(cluster, c('tempfile'))
# 
# r <- bind_rows(parLapply(cluster, system(paste0('find ', opts$outputDir, ' -name seq.fasta.out'), intern = TRUE), function(x){
#        o <- readLines(x)
#        f <- tempfile(tmpdir = '.')
#        write(o[4:length(o)], file = f)
#        r <- read.table(f, fill = TRUE, row.names = NULL, header = FALSE, quote = '', comment.char = '', stringsAsFactors = FALSE)
#        r$V1 <- as.integer(r$V1)
#        r <- r[! is.na(r$V1),]
#        file.remove(f)
#        r
#       }))
# 
# stopCluster(cluster)
# 
# unlink(paste0(opts$outputDir, '/repeatMasker'), recursive = TRUE)
# 
# names(r) <-  c('SW_score', 'percent_div', 'percent_del', 'percent_ins', 'query_seq', 'query_start', 'query_end', 'query_after', 'strand', 'repeat_name', 'repeat_class', 'repeat_start', 'repeat_end', 'repeat_after',  'ID',  'alt')
# 
# write(paste0(date(), ' Writing RepeatMasker result.'), logFile, append = TRUE)
# 
# readr::write_tsv(r, file.path(opts$outputDir, 'genomeAnnotations', paste0(opts$refGenomeName, '.repeatTable.gz')))

writeXStringSet(o, file.path(opts$outputDir, 'genome.fasta'))

# system(paste0(opts$bwa2Path, ' index -p ', opts$refGenomeName, ' ', file.path(opts$outputDir, 'genome.fasta')))
# system(paste0('mv ', paste0(list.files(pattern = paste0('^', opts$refGenomeName, '\\.'), full.names = TRUE), collapse = ' '), ' ', file.path(opts$outputDir, 'referenceGenomes', 'bwa2'), '/'))

system(paste(opts$faToTwoBitPath, file.path(opts$outputDir, 'genome.fasta'), file.path(opts$outputDir, 'referenceGenomes', 'blat', paste0(opts$refGenomeName, '.2bit'))))

invisible(file.remove(file.path(opts$outputDir, 'genome.fasta')))

# Download annotation data from UCSC.

humanGeneFilter <- function(d){
  message('Calling humanGeneFilter()')
  d$name <- sub('\\.\\d+$', '', d$name)
  message('Total gene names in source: ', n_distinct(toupper(d$name)))
  d <- subset(d, toupper(name) %in% toupper(readLines('data/geneLists/humanGeneIDs')))
  message('Gene names remaning after filter: ', n_distinct(toupper(d$name)))
  d
}

createRefSeqObjects <- function(file, label, humanGeneFilter = FALSE){
  d <- read.table(file, sep = '\t', header = FALSE, quote = '')
  
  if(length(d) == 15) d <- tibble::add_column(d, bin = NA, .before = 'V1')
  
  names(d) <- c('bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart',
                'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2',
                'cdsStartStat', 'cdsEndStat', 'exonFrames')
  
  if(humanGeneFilter) d <- humanGeneFilter(d)
  
  g <- makeGRangesFromDataFrame(d, seqnames.field = 'chrom', start.field = 'txStart',
                                end.field = 'txEnd', strand.field = 'strand',
                                keep.extra.columns = TRUE)
  
  d$exonStarts2 <- strsplit(d$exonStarts, ',')
  d$exonEnds2 <- strsplit(d$exonEnds, ',')
  
  e <- group_by(d, 1:nrow(d)) %>%
    unnest(cols = c(exonStarts2, exonEnds2)) %>%
    mutate(name = paste('exon', 1:n()),
           exonStarts = as.integer(exonStarts2),
           exonEnds = as.integer(exonEnds2)) %>%
    ungroup() %>%
    makeGRangesFromDataFrame(seqnames.field = 'chrom', start.field = 'exonStarts2',
                             end.field = 'exonEnds2', strand.field = 'strand',
                             keep.extra.columns = TRUE)
  
  saveRDS(g, file = file.path(opts$outputDir, 'genomeAnnotations', paste0(label, '.TUs.rds')))
  saveRDS(e, file = file.path(opts$outputDir, 'genomeAnnotations', paste0(label, '.exons.rds')))
}

write(paste0(date(), ' Building refSeq data.'), logFile, append = TRUE)

if('refSeqGeneAnnotations_URL' %in% names(opts)){
  system(paste0('wget -q -O ', opts$outputDir, '/', tmpFile, '.gz ', opts$refSeqGeneAnnotations_URL))
  
  if(grepl('\\.gtf', opts$refSeqGeneAnnotations_URL, ignore.case = TRUE)){
    system(paste(opts$gtfToGenePredPath, '-genePredExt', paste0(opts$outputDir, '/', tmpFile, '.gz'), paste0(opts$outputDir, '/', tmpFile, '2.gz')))
    system(paste('mv', paste0(opts$outputDir, '/', tmpFile, '2.gz'), paste0(opts$outputDir, '/', tmpFile, '.gz')))
  }
  
  createRefSeqObjects(paste0(opts$outputDir, '/', tmpFile, '.gz'), opts$refGenomeName)
  invisible(file.remove(paste0(opts$outputDir, '/', tmpFile, '.gz')))
}

if('xenoRefSeqGeneAnnotations_URL' %in% names(opts)){
  system(paste0('wget -q -O ', opts$outputDir, '/', tmpFile, '.gz ', opts$xenoRefSeqGeneAnnotations_URL))
  
  if(grepl('\\.gtf', opts$xenoRefSeqGeneAnnotations_URL, ignore.case = TRUE)){
    system(paste(opts$gtfToGenePredPath, '-genePredExt', paste0(opts$outputDir, '/', tmpFile, '.gz'), paste0(opts$outputDir, '/', tmpFile, '2.gz')))
    system(paste('mv', paste0(opts$outputDir, '/', tmpFile, '2.gz'), paste0(opts$outputDir, '/', tmpFile, '.gz')))
  }
  
  createRefSeqObjects(paste0(opts$outputDir, '/', tmpFile, '.gz'), paste0(opts$refGenomeName, '.humanXenoRef'), humanGeneFilter = TRUE)
  invisible(file.remove(paste0(opts$outputDir, '/', tmpFile, '.gz')))
}

dir.create(file.path(opts$outputDir, 'data'), showWarnings = FALSE)

system(paste('mv', file.path(opts$outputDir, 'genomeAnnotations'), file.path(opts$outputDir, 'data')))
system(paste('mv', file.path(opts$outputDir, 'referenceGenomes'), file.path(opts$outputDir, 'data')))

write(paste0(date(), ' Building final tar ball.'), logFile, append = TRUE)

# system(paste('tar cvf', file.path(opts$outputDir, paste0(opts$refGenomeName, '.tar')), ' -C ', opts$outputDir, ' data'))
# 
# unlink(paste0(opts$outputDir, '/data'), recursive = TRUE)

write(paste0(date(), ' Done.'), logFile, append = TRUE)
