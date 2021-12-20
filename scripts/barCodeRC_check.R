library(ShortRead)
library(dplyr)
library(Biostrings)
config <- yaml::read_yaml('/home/everett/projects/CSL/201218_M02973_0294_000000000-J96T7/config.yml')
I1 <- as.character(readFastq(config$demultiplex_index1ReadsFile)@sread)
samples <- readr::read_tsv(config$sampleConfigFile)
                           
d <- select(samples, subject, sample, replicate, index1.seq)

d <- bind_rows(lapply(split(d, 1:nrow(d)), function(x){
       x$barcodePercent <- sum(I1 %in% x$index1.seq)/length(I1) * 100
       x$barcodePercentRC <- sum(I1 %in% as.character(reverseComplement(DNAString(x$index1.seq))))/length(I1) * 100
       x
     }))

d