library(dplyr)
library(lubridate)
library(readr)
library(GenomicRanges)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$annotateRepeats_outputDir))

sites <- readRDS(file.path(opt$outputDir, opt$annotateRepeats_inputFile))

samples <- loadSamples()

# Assuming replicates collapsed.
sites$uniqueSample <- paste0(sites$trial, '~', sites$subject, '~', sites$sample)
samples$uniqueSample <- paste0(samples$trial, '~', samples$subject, '~', samples$sample)

# Add refGenome to sites.
sites$refGenome.id <- NULL
sites <- distinct(left_join(sites, select(samples, uniqueSample, refGenome.id), by = 'uniqueSample'))

invisible(lapply(unique(sites$refGenome.id), function(x){
  if(! file.exists(file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x, '.repeatTable.gz')))){
    stop('Error - ', file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x, '.repeatTable.gz')), ' does not exists.')
  }
}))

sites <- bind_rows(lapply(unique(sites$refGenome.id), function(x){
  r <- readr::read_tsv(file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x, '.repeatTable.gz')))
  r$strand <- sub('C', '-', r$strand)
  r <- subset(r, strand %in% c('+', '-'))
  g <- makeGRangesFromDataFrame(select(r, query_seq, query_start, query_end, strand, repeat_name, repeat_class), 
                                seqnames.field = 'query_seq', start.field = 'query_start', end.field = 'query_end', 
                                strand.field = 'strand', keep.extra.columns = TRUE)
  
  s <- subset(sites, refGenome.id == x)
  
  a <- strsplit(s$posid, '[\\+\\-]')
  s$chromosome <- unlist(lapply(a, '[', 1))
  s$strand <- stringr::str_extract(s$posid, '[\\+\\-]')
  s$position <- unlist(lapply(a, '[', 2))
  
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
   return(left_join(s, select(r, posid, repeat_name, repeat_class), by = 'posid'))
  } else {
    return(s)
  }
})) 

saveRDS(select(sites, -chromosome, -strand, -position, -refGenome.id, -uniqueSample), file = file.path(opt$outputDir, opt$annotateRepeats_outputDir, opt$annotateRepeats_outputFile))



if(file.exists(file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitClusters.rds'))){
  o <- readRDS(file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitClusters.rds'))
  
  o$refGenome.id <- NULL
  samples$s <- paste(samples$trial, samples$subject, samples$sample)
  o$s <- paste(o$trial, o$subject, o$sample)
  o <- left_join(o, distinct(select(samples, s, refGenome.id)), by = 's')
  
  
  
  o <- bind_rows(lapply(unique(o$refGenome.id), function(x){
    r <- readr::read_tsv(file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x, '.repeatTable.gz')))
    r$strand <- sub('C', '-', r$strand)
    r <- subset(r, strand %in% c('+', '-'))
    g <- makeGRangesFromDataFrame(select(r, query_seq, query_start, query_end, strand, repeat_name, repeat_class), 
                                  seqnames.field = 'query_seq', start.field = 'query_start', end.field = 'query_end', 
                                  strand.field = 'strand', keep.extra.columns = TRUE)
    
    s <- subset(o, refGenome.id == x)
    
    bind_rows(lapply(1:nrow(s), function(x){
      x <- s[x,]
      nodes <- unlist(x$clusters)
    
      a <- strsplit(nodes, '[\\+\\-]')
      d <- tibble(chromosome = unlist(lapply(a, '[', 1)),
                  strand = stringr::str_extract(nodes, '[\\+\\-]'),
                  position = unlist(lapply(a, '[', 2)),
                  posid = nodes)
    
      s2 <- makeGRangesFromDataFrame(select(d, chromosome, strand, position, posid), 
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
      topRepeatClass <- names(sort(table(r$repeat_class), decreasing = TRUE)[1])
      x$topRepeatClass <- topRepeatClass
      x$percentNodesInTopClass <- sprintf("%.1f%%", n_distinct(subset(r, repeat_class = topRepeatClass)$posid) / n_distinct(r$posid)*100)
    } else {
      x$topRepeatClass <- NA
      x$percentNodesInTopClass <- NA
    }
      
    return(dplyr::select(x, -s, -refGenome.id))
  })) 
 }))
  
  saveRDS(o, file = file.path(opt$outputDir, opt$annotateRepeats_outputDir, 'multiHitClusters.rds'))
}

q(save = 'no', status = 0, runLast = FALSE) 
