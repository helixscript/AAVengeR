library(dplyr)
library(ggplot2)
library(scales)
library(ShortRead)

# R1 251   R2 
path <- '/data/sequencingData/220222_M03249_0243_000000000-JHTY6/Undetermined_S0_L001_R2_001.fastq.gz'

r <- readFastq(path)
r <- trimTailw(r, 2, '?', 5)
r <- r[width(r) >= 200]

s <- as.character(subseq(r@sread, 1, 200))


readSamplePlot <- function(reads, n){
  ds <- sample(unique(reads), n)
  dp <- lapply(strsplit(sort(ds), ''), function(x){ tibble(base = x, n = 1:length(x)) })
  dp <- bind_rows(mapply(function(x, n){ x$read <- n; x}, dp, 1:length(dp), SIMPLIFY = FALSE))
  dp$base <- factor(dp$base, levels = c('A', 'T', 'C', 'G', 'N'))
  
  ggplot(dp, aes(n, read, fill = base)) + theme_bw() + geom_tile() +
    scale_fill_manual(values =  c('red', 'green', 'blue', 'gold', 'gray50')) +
    scale_x_continuous(limits = c(0, width(reads[1])), expand = c(0, 0)) +
    scale_y_continuous(label=comma, limits = c(0, n), expand = c(0, 0)) +
    labs(x = 'Position', y = 'Reads') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}

readSamplePlot(s, 10000)

