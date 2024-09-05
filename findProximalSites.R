library(dplyr)
library(GenomicRanges)
expand <- 3

args <- commandArgs(trailingOnly=TRUE)

if(length(args) != 2) stop('Error - two command line arguments were expectes; input and output files.')
if(! file.exists(args[1])) stop('Error - input file could not be found.')

sites <- readRDS(args[1])

oo <- tidyr::separate(tibble(sample = sites$sample, posid = sites$posid), posid, c('seqnames', 'start'), sep = '[\\+\\-]', remove = FALSE)
oo$strand <- stringr::str_extract(oo$posid, '[\\+\\-]')
oo$start <- as.integer(sub('\\.\\d+$', '', oo$start))
oo$end <- oo$start

r <- bind_rows(lapply(split(oo, oo$sample), function(x){
       clusterNum <- 0
       g <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE) + expand

       o <- GenomicRanges::reduce(g, with.revmap = TRUE)

       bind_rows(lapply(o$revmap, function(i){
         if(length(i) > 1){
           clusterNum <<- clusterNum + 1
           m <- subset(sites, posid %in% g[i]$posid & sample == x$sample[1])
           t <- tibble(sample = x$sample[1], 
                       cluster = clusterNum, 
                       sites = m$posid, 
                       reads = m$reads, 
                       sonicLengths = m$sonicLengths,
                       repLeaderSeq = m$repLeaderSeq)
         } else {
            return(tibble())
         }
       }))
     }))

readr::write_tsv(r, args[2])

q()
