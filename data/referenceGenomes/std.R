library(Biostrings)

g <- readDNAStringSet('macFas5.2bit.org.fasta')
r <- c('chrX', 'chrY', paste0('chr', 1:50))
g <- g[names(g) %in% r]
writeXStringSet(g, 'macFas5.fasta')