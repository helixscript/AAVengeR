library(GenomicRanges)
source('../stdPos.lib.R')

o <- makeGRangesFromDataFrame(data.frame(
          seqnames = 'chr1', 
          strand = '+', 
          start  = c(98,99,100,101,102,200,201,202,203,204),
          end    = c(150,151,149,150,150,250,251,250,250,249),
          reads  = c(2,4,99,10,2,1,4,80,10,2)), keep.extra.columns = TRUE)

o <- standardize_sites(o, counts.col = 'reads', sata.gap = 5)
o <- refine_breakpoints(o, counts.col = 'reads', sata.gap = 3)
o