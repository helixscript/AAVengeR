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

if(! opt$annotateRepeats_addAfter %in% names(sites)){
  write(c(paste(now(), paste0('   Error - ', opt$annotateRepeats_addAfter, ' is not a column in your input data frame.'))), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}

invisible(lapply(unique(sites$refGenome), function(x){
  if(! file.exists(file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x, '.repeatTable.gz')))){
    stop('Error - ', file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x, '.repeatTable.gz')), ' does not exists.')
  }
}))

sites <- bind_rows(lapply(unique(sites$refGenome), function(x){
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

sites <- dplyr::relocate(sites, repeat_name, repeat_class, .after = opt$annotateRepeats_addAfter) %>% 
         dplyr::select(-chromosome, -strand, -position)

saveRDS(sites, file = file.path(opt$outputDir, opt$annotateRepeats_outputDir, 'sites.rds'))
openxlsx::write.xlsx(sites, file = file.path(opt$outputDir, opt$annotateRepeats_outputDir, 'sites.xlsx'))
readr::write_csv(sites, file.path(opt$outputDir, opt$annotateRepeats_outputDir, 'sites.csv'))


if(file.exists(file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitClusters.rds'))){
  o <- readRDS(file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitClusters.rds'))
  
  o <- left_join(o, distinct(select(sites, sample, refGenome)), by = 'sample')
  
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
      
    x
  }))
    
 }))
  
 o2 <- bind_cols(o[,1:3], o[,14:16], o[,5:13], o[,4]) %>% arrange(desc(reads))
 saveRDS(o2, file = file.path(opt$outputDir, opt$annotateRepeats_outputDir, 'multiHitClusters.rds'))
 readr::write_csv(o2, file.path(opt$outputDir, opt$annotateRepeats_outputDir, 'multiHitClusters.csv'))
 openxlsx::write.xlsx(o2, file.path(opt$outputDir, opt$annotateRepeats_outputDir, 'multiHitClusters.xlsx'))
}

q(save = 'no', status = 0, runLast = FALSE) 
