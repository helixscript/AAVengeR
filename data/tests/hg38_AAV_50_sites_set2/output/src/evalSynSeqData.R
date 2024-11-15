#!/usr/bin/Rscript

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2) stop('usage:  evalSynSeqData.R <analysis dir path> <output file>')

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GenomicRanges))

posidToGrange <- function(x, expand = 0){
  oo <- tibble(x = x)
  oo <- tidyr::separate(oo, x, c('seqname', 'start'), sep = '[\\+\\-]', remove = FALSE)
  oo$strand <- stringr::str_extract(oo$x, '[\\+\\-]')
  oo$end <- oo$start
  makeGRangesFromDataFrame(oo) + expand
}

r <- tibble()

if(file.exists(file.path(args[1], 'output', 'core', 'sites.rds'))){
    o  <- readRDS(file.path(args[1], 'output', 'core', 'sites.rds'))
    t <- readr::read_tsv(file.path(args[1], 'truth.tsv'), show_col_types = FALSE)
        
    tg <- posidToGrange(t$posid, expand = 3)
    og <- posidToGrange(o$posid, expand = 3)
        
    i <- findOverlaps(og, tg)
        
    # Find percentage of unique called sites.
    p <- n_distinct(subjectHits(i)) / n_distinct(t$posid)
        
    missingSites <- t[which(! c(1:nrow(t)) %in% subjectHits(i)),]$posid
        
    z <- o[which(! c(1:nrow(o)) %in% queryHits(i)),]$posid
    unexpectedSites <- 0
    if(length(z) > 0) unexpectedSites <- n_distinct(z)
        
    if(file.exists(file.path(args[1], 'output', 'core', 'multiHitClusters.rds'))){
          message('Found multiHit result.')
          m <- readRDS(file.path(args[1], 'output', 'core', 'multiHitClusters.rds'))
         
          missingSites <- posidToGrange(unique(sub('\\.\\d+$', '', missingSites)), expand = 3)
          multiHitSites <- posidToGrange(unique(sub('\\.\\d+$', '', unlist(m$posids))), expand = 3) 
          
          i2 <- findOverlaps(missingSites, multiHitSites)
          p2 <- (n_distinct(subjectHits(i)) + n_distinct(queryHits(i2))) / n_distinct(t$posid)
    } else {
          
          p2 <- NA
    }
        
    e <- unlist(lapply(unique(subjectHits(i)), function(k){
               r <- o[queryHits(i[subjectHits(i) == k]),]$repLeaderSeq
               e <- t[k,]$leaderSeq
               stringdist::stringdist(r, e)
          }))
        
    
    r <- tibble(exp = gsub('_', '', unlist(stringr::str_extract_all(args[1], '_[^\\d]+_')))[1],
                  nSites = as.integer(gsub('_', '', (stringr::str_extract_all(args[1], '_\\d+_')))), 
                  set = stringr::str_extract(args[1], '\\d+$'), 
                  percentUniqueRecovery = sprintf("%.1f%%", p*100), 
                  percentTotalRecovery = ifelse(is.na(p2), sprintf("%.1f%%", p*100), sprintf("%.1f%%", p2*100)), 
                  unexpectedSites = unexpectedSites,
                  leaderSeqDist.mean = round(mean(e), 2),
                  leaderSeqDist.sd = round(sd(e), 2)) %>% arrange(as.integer(nSites), exp, as.integer(set))
} else {
  message('Error -- could not locate analysis folder.')
}

readr::write_tsv(r, args[2])

q()


