library(GenomicRanges)
library(dplyr)
source('../lib.R')

o <- makeGRangesFromDataFrame(data.frame(
       seqnames = 'chr1', 
       strand = '+', 
       start  = c(98,99,100,101,102,200,201,202,203,204),
       end    = c(150,151,149,150,150,250,251,250,250,249),
       reads  = c(2,4,99,10,2,1,4,80,10,2)), keep.extra.columns = TRUE)


o1 <- refine_breakpoints(o, counts.col = 'reads')
o2 <- standardize_sites(o1, counts.col = 'reads')

expected <- makeGRangesFromDataFrame(data.frame(
              seqnames = 'chr1', 
              strand = '+', 
              start = c(100, 202),
              end = c(149, 250)),  keep.extra.columns = TRUE)

message(ifelse(all(GenomicRanges::reduce(o2) == expected), 'pass', 'fail'))