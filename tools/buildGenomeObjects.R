#!/usr/bin/Rscript

for (p in c('yaml', 'dplyr', 'Biostrings', 'parallel', 'GenomicRanges', 'tidyr', 'parallel')) suppressPackageStartupMessages(library(p, character.only = TRUE))

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
### dir.create(file.path(opts$outputDir, 'referenceGenomes', 'bwa2'), showWarnings = FALSE)

# Download genome data from UCSC.
write(paste0(date(), ' Downloading genome data.'), logFile, append = TRUE)
system(paste0('wget -q -O ', opts$outputDir, '/', tmpFile, '.gz ', opts$refGenomeSeq_URL))
system(paste0('gunzip ', opts$outputDir, '/', tmpFile, '.gz'))

# Read in genome NTs and limit to standard chromosomes.
o <- readDNAStringSet(paste0(opts$outputDir, '/', tmpFile))
invisible(file.remove(paste0(opts$outputDir, '/', tmpFile)))

if(opts$forceStdChromosomes){
  if(opts$refGenomeName =='felCat9'){
    allowedChromosomes <- c("chrA1", "chrA2", "chrA3", "chrB1", "chrB2", "chrB3", "chrB4", "chrC1", "chrC2", "chrD1", "chrD2", "chrD3", "chrD4", "chrE1", "chrE2", "chrE3", "chrF1", "chrF2", "chrM", "chrX")
  } else {
    allowedChromosomes <- c(paste0('chr', as.roman(1:100)), 'chrY', 'chrM')
  }
  
  o <- o[names(o) %in% allowedChromosomes]
}

# Run RepeatMasker.
dir.create(paste0(opts$outputDir, '/repeatMasker'), showWarnings = FALSE)

n <- 1
invisible(lapply(split(o, names(o)), function(x){
  write(paste0(date(), ' Starting RepeatMasker for chromosome ', unique(names(o)), ' - ',   n, '/', n_distinct(names(o)), ' length: ', width(x)), logFile, append = TRUE)
  n <<- n+1
  tmpFile <- paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')
  dir.create(paste0(opts$outputDir, '/repeatMasker/', tmpFile), showWarnings = FALSE)
  writeXStringSet(x, paste0(opts$outputDir, '/repeatMasker/', tmpFile, '/seq.fasta'))
  comm <- paste0(opts$repeatMaskerPath, ' -s -pa ', opts$CPUs, ' -e rmblast -species "',
                 opts$repeatMaskerSpecies, '" -dir ', paste0(opts$outputDir, '/repeatMasker/', tmpFile), ' ',
                 paste0(opts$outputDir, '/repeatMasker/', tmpFile, '/seq.fasta'))
  message(comm)
  system(comm)
}))

write(paste0(date(), ' Reading RepeatMasker outputs.'), logFile, append = TRUE)

cluster <- makeCluster(opts$CPUs)
clusterExport(cluster, c('tempfile'))

r <- bind_rows(parLapply(cluster, system(paste0('find ', opts$outputDir, ' -name seq.fasta.out'), intern = TRUE), function(x){
       o <- readLines(x)
       f <- tempfile(tmpdir = '.')
       write(o[4:length(o)], file = f)
       r <- read.table(f, fill = TRUE, row.names = NULL, header = FALSE, quote = '', comment.char = '', stringsAsFactors = FALSE)
       r$V1 <- as.integer(r$V1)
       r <- r[! is.na(r$V1),]
       file.remove(f)
       r
      }))

stopCluster(cluster)

unlink(paste0(opts$outputDir, '/repeatMasker'), recursive = TRUE)

repeatMaskerColNames <- c('SW_score', 'percent_div', 'percent_del', 'percent_ins', 'query_seq', 'query_start', 'query_end', 'query_after', 'strand', 'repeat_name', 'repeat_class', 'repeat_start', 'repeat_end', 'repeat_after',  'ID',  'alt')
if(length(r) == 15) repeatMaskerColNames <- c('SW_score', 'percent_div', 'percent_del', 'percent_ins', 'query_seq', 'query_start', 'query_end', 'query_after', 'strand', 'repeat_name', 'repeat_class', 'repeat_start', 'repeat_end', 'repeat_after',  'ID')
names(r) <-  repeatMaskerColNames

write(paste0(date(), ' Writing RepeatMasker result.'), logFile, append = TRUE)

readr::write_tsv(r, file.path(opts$outputDir, 'genomeAnnotations', paste0(opts$refGenomeName, '.repeatTable.gz')))

writeXStringSet(o, file.path(opts$outputDir, 'genome.fasta'))

### system(paste0(opts$bwa2Path, ' index -p ', opts$refGenomeName, ' ', file.path(opts$outputDir, 'genome.fasta')))
### system(paste0('mv ', paste0(list.files(pattern = paste0('^', opts$refGenomeName, '\\.'), full.names = TRUE), collapse = ' '), ' ', file.path(opts$outputDir, 'referenceGenomes', 'bwa2'), '/'))

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
  library(tidyverse)
  d <- read.table(file, sep = '\t', header = FALSE, quote = '')
  
  # Expected
  # 'bin',  'name',      'chrom',   'strand',  'txStart',      'txEnd',         'cdsStart',     'cdsEnd',   'exonCount',   'exonStarts',                                                                                                                                                    'exonEnds',                                                                                                           'score',   'name2',    'cdsStartStat', 'cdsEndStat',       'exonFrames'
  # 1145    NM_000477.7     chr4      +         73404286        73421482        73404327        73420298        15           73404286,73405115,73406628,73408593,73409354,73410311,73411995,73413419,73415034,73416255,73417530,73418087,73419506,73420253,73421091,    73404406,73405173,73406761,73408805,73409487,73410409,73412125,73413634,73415167,73416353,73417669,73418311,73419639,73420321,73421482,      0         ALB            cmpl           cmpl          0,1,2,0,2,0,2,0,2,0,2,0,2,0,-1,

  if(any(d$V6 == '+') & any(d$V6 == '-')){
    # BED 12+ format.
    # Convert to expect format above.
    names(d) <- c("chrom","txStart","txEnd","name","score","strand",
                  "cdsStart","cdsEnd","itemRgb","exonCount","blockSizes","blockStarts",
                  "name2","cdsStartStat","cdsEndStat","exonFrames","V17","transcript_id","gene_name","V20")
    
    d <- d %>%
      mutate(
        across(c(txStart, txEnd, cdsStart, cdsEnd, exonCount), as.integer),
        blockSizes  = str_remove(blockSizes,  ",+$"),     # drop trailing comma
        blockStarts = str_remove(blockStarts, ",+$"),
        bs = str_split(blockSizes,  ","),
        bt = str_split(blockStarts, ",")
      ) %>%
      pmap_dfr(function(chrom, txStart, txEnd, name, score, strand,
                        cdsStart, cdsEnd, itemRgb, exonCount, blockSizes, blockStarts,
                        name2, cdsStartStat, cdsEndStat, exonFrames, ...){
        bs <- as.integer(unlist(strsplit(blockSizes,  ",")))
        bt <- as.integer(unlist(strsplit(blockStarts, ",")))
        n  <- min(exonCount, length(bs), length(bt))
        exStarts <- paste(txStart + bt[seq_len(n)], collapse = ",")
        exEnds   <- paste(txStart + bt[seq_len(n)] + bs[seq_len(n)], collapse = ",")
        tibble(
          bin = 0L, name, chrom, strand,
          txStart, txEnd, cdsStart, cdsEnd,
          exonCount,
          exonStarts = paste0(exStarts, ","),
          exonEnds   = paste0(exEnds,   ","),
          score, name2, cdsStartStat, cdsEndStat, exonFrames
        )
      })
    
  } else {
    if(length(d) == 15) d <- tibble::add_column(d, bin = NA, .before = 'V1')
  
    names(d) <- c('bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart',
                  'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2',
                  'cdsStartStat', 'cdsEndStat', 'exonFrames')
  }
  
  if(humanGeneFilter) d <- humanGeneFilter(d)
  
  g <- makeGRangesFromDataFrame(d, seqnames.field = 'chrom', start.field = 'txStart',
                                end.field = 'txEnd', strand.field = 'strand',
                                keep.extra.columns = TRUE)
  
  d$exonStarts <- sub(',\\s*$', '', d$exonStarts)
  d$exonEnds <- sub(',\\s*$', '', d$exonEnds)

  e <- d %>%
    mutate(row_id = row_number()) %>%     
    group_by(row_id) %>%
    tidyr::separate_rows(exonStarts, exonEnds, sep = ",") %>%
    mutate(name       = paste("exon", row_number()),
           exonStarts = as.integer(exonStarts),
           exonEnds   = as.integer(exonEnds)) %>%
    ungroup() %>%
    select(-row_id) %>%
    makeGRangesFromDataFrame(seqnames.field = 'chrom', start.field = 'exonStarts',
                             end.field = 'exonEnds', strand.field = 'strand',
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

system(paste('tar cvf', file.path(opts$outputDir, paste0(opts$refGenomeName, '.tar')), ' -C ', opts$outputDir, ' data'))

unlink(paste0(opts$outputDir, '/data'), recursive = TRUE)

write(paste0(date(), ' Done.'), logFile, append = TRUE)
