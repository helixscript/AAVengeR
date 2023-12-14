library(Biostrings)

fastaFile <- '../data/referenceGenomes/bwa2/chlSab2.fa'
outputFile <- '../data/referenceGenomes/bwa2/chlSab2.std.fa'

g <- readDNAStringSet(fastaFile)

writeXStringSet(g[names(g) %in% paste0('chr', c(1:100, 'X', 'Y'))], outputFile)