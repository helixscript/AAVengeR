library(dplyr)
library(stringr)
library(readr)
library(Biostrings)

# This script reads in a vector sequence and a second FASTA formatted file
# and aligns the vector against the FASTA file. Significant alignments to the 
# second file are masked with lowercase letters and written out to a second 
# vector file.
#
# Input vectors should already be masked with repeatMasker, eg.
# RepeatMasker -s -xsmall -species "Mus musculus" -dir tmp vector.ff

vectorIN <- 'data/vectors/Peranteau-AAV8-CAG-GFP.fasta'  
vectorOUT <- 'data/vectors/Peranteau-AAV8-CAG-GFP.filtered.fasta'
blastDB <- '/home/ubuntu/projects/Peranteau_AAV/guideRegions.fasta'
tmpDir <- '~/tmp'
minPercentID <- 93
minAlnLength <- 15

invisible(file.remove(list.files(tmpDir, full.names = TRUE)))

v <- readLines(vectorIN)
v.id <- v[1]
v <- paste0(v[2:length(v)], collapse = '') # Assume Fasta format. Can not use read fasta functions or will lose lower case filtering.

system(paste0('makeblastdb -dbtype nucl -in ', blastDB, ' -out ', file.path(tmpDir, 'db')))
system(paste0('blastn -dust no -soft_masking false -word_size 7 -evalue 100 -outfmt 6 -query ', 
              vectorIN, ' -db ', file.path(tmpDir, 'db'), ' -out ', file.path(tmpDir, 'blastOut')))

if(file.info(file.path(tmpDir, 'blastOut'))$size > 0){
  b <- read.table(file.path(tmpDir, 'blastOut'), sep = '\t', header = FALSE)
  names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
  b <- filter(b, length >= minAlnLength & pident >= minPercentID & gapopen <= 1)
  
  if(nrow(b) > 0){
    invisible(lapply(split(b, 1:nrow(b)), function(x){
      message(tolower(substr(v, x$qstart, x$qend)))
      substr(v, x$qstart, x$qend) <<- tolower(substr(v, x$qstart, x$qend))
    }))
  }
}

write(c(v.id, v), file = vectorOUT)
