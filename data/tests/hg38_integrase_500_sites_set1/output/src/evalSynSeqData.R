library(dplyr)
library(GenomicRanges)

args <- commandArgs(trailingOnly=TRUE)

posidToGrange <- function(x, expand = 0){
  oo <- tibble(x = x)
  oo <- tidyr::separate(oo, x, c('seqname', 'start'), sep = '[\\+\\-]', remove = FALSE)
  oo$strand <- stringr::str_extract(oo$x, '[\\+\\-]')
  oo$end <- oo$start
  makeGRangesFromDataFrame(oo) + expand
}

r <- bind_rows(lapply(d, function(x){
      message(x)
      if(file.exists(file.path(x, 'output', 'core', 'sites.rds'))){

        o  <- readRDS(file.path(x, 'output', 'core', 'sites.rds'))
        t <- readr::read_tsv(file.path(x, 'truth.tsv'), show_col_types = FALSE)
        
        
        tg <- posidToGrange(t$posid, expand = 3)
        og <- posidToGrange(o$posid, expand = 3)
        
        i <- findOverlaps(og, tg)
        
        # Find percentage of unique called sites.
        p <- n_distinct(subjectHits(i)) / n_distinct(t$posid)
        
        missingSites <- t[which(! c(1:nrow(t)) %in% subjectHits(i)),]$posid
        
        z <- o[which(! c(1:nrow(o)) %in% queryHits(i)),]$posid
        unexpectedSites <- 0
        if(length(z) > 0) unexpectedSites <- n_distinct(z)
        
        if(file.exists(file.path(x, 'output', 'core', 'multiHitClusters.rds'))){
          m <- readRDS(file.path(x, 'output', 'core', 'multiHitClusters.rds'))
         
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
        
        return(tibble(exp = gsub('_', '', unlist(stringr::str_extract_all(x, '_[^\\d]+_')))[1],
                      nSites = as.integer(gsub('_', '', (stringr::str_extract_all(x, '_\\d+_')))), 
                      set = stringr::str_extract(x, '\\d+$'), 
                      percentUniqueRecovery = sprintf("%.1f%%", p*100), 
                      percentTotalRecovery = ifelse(is.na(p2), sprintf("%.1f%%", p*100), sprintf("%.1f%%", p2*100)), 
                      unexpectedSites = unexpectedSites,
                      leaderSeqDist.mean = round(mean(e), 2),
                      leaderSeqDist.sd = round(sd(e), 2)))
     } else {
       return(tibble())
     }
   })) %>% filter(nSites >= 100) %>% arrange(as.integer(nSites), exp, as.integer(set))

