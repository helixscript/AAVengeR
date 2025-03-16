#!/usr/bin/Rscript

# AAVengeR/annotateRepeats.R
# John K. Everett, Ph.D.
# 
# This scripts accepts integration sites buildSites or any module following buildSites.
# For each genome, a repeat table created by RepeatMasker should be available in AAVengeR's
# data directory. This script will report both the type and class of repeats intersecting 
# with integrations.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(GenomicRanges))

# Read in the configuration file and perform basic sanity checks.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) stop('Expected at least one command line argument')
source(file.path(yaml::read_yaml(args[1])$softwareDir, 'lib.R'))
opt <- startModule(args)

createOuputDir()
dir.create(file.path(opt$outputDir, opt$annotateRepeats_outputDir), showWarnings = FALSE)

# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, opt$annotateRepeats_outputDir, 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)
write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)

quitOnErorr <- function(msg){
  updateLog(msg)
  message(msg)
  message(paste0('See log for more details: ', opt$defaultLogFile))
  q(save = 'no', status = 1, runLast = FALSE) 
}

updateLog('Starting annotateRepeats.')

if(! file.exists(file.path(opt$outputDir, opt$annotateRepeats_inputFile))){
  updateLog('Error - input file does not exist.')
  q(save = 'no', status = 1, runLast = FALSE)
}

sites <- readRDS(file.path(opt$outputDir, opt$annotateRepeats_inputFile))

if(nrow(sites) == 0){
  updateLog('Error - sites file was read and contains zero rows of data.')
  q(save = 'no', status = 1, runLast = FALSE)
}

if(! opt$annotateRepeats_addAfter %in% names(sites)){
  updateLog(paste0('Error - ', opt$annotateRepeats_addAfter, ' is not a column in your input data frame.'))
  q(save = 'no', status = 1, runLast = FALSE) 
}

invisible(lapply(unique(sites$refGenome), function(x){
  if(! file.exists(file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x, '.repeatTable.gz')))){
    updateLog(paste0('Error - ', file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x, '.repeatTable.gz')), ' does not exists.'))
    stop('Error - ', file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x, '.repeatTable.gz')), ' does not exists.')
  }
}))

sites <- bind_rows(lapply(unique(sites$refGenome), function(x){
  updateLog(paste0('Loading annotation table: ', file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x, '.repeatTable.gz'))))
  r <- readr::read_tsv(file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x, '.repeatTable.gz')))
  
  r$strand <- sub('C', '-', r$strand)
  r <- subset(r, strand %in% c('+', '-'))
  g <- makeGRangesFromDataFrame(select(r, query_seq, query_start, query_end, strand, repeat_name, repeat_class), 
                                seqnames.field = 'query_seq', start.field = 'query_start', end.field = 'query_end', 
                                strand.field = 'strand', keep.extra.columns = TRUE)
  
  s <- subset(sites, refGenome == x)
  
  a <- strsplit(s$posid, '[\\+\\-]')
  s$chromosome <- unlist(lapply(a, '[', 1))
  s$strand <- stringr::str_extract(s$posid, '[\\+\\-]')
  s$position <- unlist(lapply(a, '[', 2))
  
  s$position <- sub('\\.\\S+$', '', s$position)
  
  s2 <- makeGRangesFromDataFrame(select(s, chromosome, strand, position, posid), 
                                start.field = 'position', end.field = 'position', keep.extra.columns = TRUE)
  
  o <- GenomicRanges::findOverlaps(s2, g, ignore.strand = TRUE)
  k <- tibble(posid = s2[queryHits(o)]$posid, 
              repeat_name = g[subjectHits(o)]$repeat_name,
              repeat_class = g[subjectHits(o)]$repeat_class)
  
  r <- bind_rows(lapply(split(k, k$posid), function(x){
         x$repeat_name  <- paste0(sort(unique(x$repeat_name)), collapse = ',')
         x$repeat_class <- paste0(sort(unique(x$repeat_class)), collapse = ',')
         x[1,]
       }))
  
  if(nrow(r) > 0){
   return(left_join(s, select(r, posid, repeat_name, repeat_class), by = 'posid', relationship = 'many-to-many'))
  } else {
    return(s)
  }
})) 

sites <- dplyr::relocate(sites, repeat_name, repeat_class, .after = opt$annotateRepeats_addAfter) %>% 
         dplyr::select(-chromosome, -strand, -position)

sites <- arrange(sites, desc(sonicLengths))

if(opt$databaseConfigGroup != 'none'){
  suppressPackageStartupMessages(library(RMariaDB))
  uploadSitesToDB(sites)
}

saveRDS(sites, file = file.path(opt$outputDir, opt$annotateRepeats_outputDir, 'sites.rds'))
openxlsx::write.xlsx(sites, file = file.path(opt$outputDir, opt$annotateRepeats_outputDir, 'sites.xlsx'))
readr::write_tsv(sites, file.path(opt$outputDir, opt$annotateRepeats_outputDir, 'sites.tsv.gz'))


if(file.exists(file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitClusters.rds'))){
  updateLog(paste0('Updating multi-hit clusters found at: ', file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitClusters.rds')))
  
  o <- readRDS(file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitClusters.rds'))
  
  o <- left_join(o, distinct(select(sites, sample, refGenome)), by = 'sample', relationship = 'many-to-many')
  
  o <- bind_rows(lapply(unique(o$refGenome), function(x){
    r <- readr::read_tsv(file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x, '.repeatTable.gz')))
    r$strand <- sub('C', '-', r$strand)
    r <- subset(r, strand %in% c('+', '-'))
    g <- makeGRangesFromDataFrame(select(r, query_seq, query_start, query_end, strand, repeat_name, repeat_class), 
                                  seqnames.field = 'query_seq', start.field = 'query_start', end.field = 'query_end', 
                                  strand.field = 'strand', keep.extra.columns = TRUE)
    
    s <- subset(o, refGenome == x)
    
    bind_rows(lapply(1:nrow(s), function(x){
      x <- s[x,]
      nodes <- unlist(x$posids)
    
      a <- strsplit(nodes, '[\\+\\-]')
      d <- tibble(chromosome = unlist(lapply(a, '[', 1)),
                  strand = stringr::str_extract(nodes, '[\\+\\-]'),
                  position = unlist(lapply(a, '[', 2)),
                  posid = nodes)
    
      s2 <- makeGRangesFromDataFrame(select(d, chromosome, strand, position, posid), 
                                   start.field = 'position', end.field = 'position', keep.extra.columns = TRUE)
      
      s2$position <- sub('\\.\\S+$', '', s2$position)
    
      o <- GenomicRanges::findOverlaps(s2, g, ignore.strand = TRUE)
      k <- tibble(posid = s2[queryHits(o)]$posid, 
                  repeat_name = g[subjectHits(o)]$repeat_name,
                  repeat_class = g[subjectHits(o)]$repeat_class)
    
      r <- bind_rows(lapply(split(k, k$posid), function(x){
             x$repeat_name  <- paste0(sort(unique(x$repeat_name)), collapse = ',')
             x$repeat_class <- paste0(sort(unique(x$repeat_class)), collapse = ',')
             x[1,]
           }))
    
    if(nrow(r) > 0){
      topRepeatClass <- names(sort(table(r$repeat_class), decreasing = TRUE)[1])
      x$topRepeatClass <- topRepeatClass
      x$percentNodesInTopClass <- sprintf("%.1f%%", n_distinct(subset(r, repeat_class = topRepeatClass)$posid) / n_distinct(r$posid)*100)
    } else {
      x$topRepeatClass <- NA
      x$percentNodesInTopClass <- NA
    }
      
    x
  }))
    
 }))
  
 saveRDS(o, file = file.path(opt$outputDir, opt$annotateRepeats_outputDir, 'multiHitClusters.rds'))
 openxlsx::write.xlsx(o, file.path(opt$outputDir, opt$annotateRepeats_outputDir, 'multiHitClusters.xlsx'))
 readr::write_tsv(o, file.path(opt$outputDir, opt$annotateRepeats_outputDir, 'multiHitClusters.tsv.gz'))
}

updateLog('annotateRepeats completed.')

q(save = 'no', status = 0, runLast = FALSE) 
