library(GenomicRanges)
library(ggplot2)
library(dplyr)
source('../stdPos.lib.R')

# Simulate standardization of dual detection sites.
d <- makeGRangesFromDataFrame(data.frame(seqnames = 'chr1', 
                strand = c('+','+','+','+','+',  '-','-','-','-','-'),
                start  = c(98,99,100,101,102,    90,91,97,98,100),
                end    = c(150,151,165,166,171,  102,103,104,105,106),
                reads  = c(2,4,99,10,2,              1,4,80,10,2)), keep.extra.columns = TRUE)

o1 <- standardize_sites(d, counts.col = 'reads', sata.gap = 5)
o2 <- refine_breakpoints(o1, counts.col = 'reads', sata.gap = 3)

# Simulate two postions next to one another with the same orientation.

d <- makeGRangesFromDataFrame(data.frame(seqnames = 'chr1', 
                                         strand = c('+','+','+','+','+', '+','+','+','+','+'),
                                         start  = c(98,99,100,101,102,    108,109,110,111,112),
                                         end    = c(150,151,160,161,171,  160,161,170,171,181),
                                         reads  = c(2,4,99,10,2,              1,4,80,10,2)), keep.extra.columns = TRUE)

o1 <- standardize_sites(d, counts.col = 'reads', sata.gap = 5)
o2 <- refine_breakpoints(o1, counts.col = 'reads', sata.gap = 3)
