library(Biostrings)

# Create SNP within a provided sequence to create a faux MSA.
# Select a range such that the SNPs are not introduced into the beginning 
# or end ofthe sequence since we want the resulting HMM to strongly match 
# these regions.

seq <- 'TCAGCGGGGGTCTTTCA'
rangeStart <- 6
rangeEnd <- 14
num <- 100
pool <- c('A', 'T', 'C', 'G')

set.seed(1)
r <- sapply(1:num, function(x){
     n <- sample(rangeStart:rangeEnd, 1)
     o <- unlist(strsplit(seq, ''))
     o[n] <- sample(pool[pool != o[n]], 1)
     paste0(o, collapse = '')
})

d <- DNAStringSet(r)
names(d) <- paste0('s', 1:length(d))
writeXStringSet(d, file = 'simulateSeqErrorMSA.fasta')