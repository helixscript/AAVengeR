library(dplyr)
library(ShortRead)

I1_path <- '~/projects/Ype_deJong_AAV/data/Undetermined_S0_I1_001.fastq.gz'
R1_path <- '~/projects/Ype_deJong_AAV/data/Undetermined_S0_R1_001.fastq.gz'
R2_path <- '~/projects/Ype_deJong_AAV/data/Undetermined_S0_R2_001.fastq.gz'

I1_out <- '~/projects/Ype_deJong_AAV/data/I1_controlsOnly.fastq.gz'
R1_out <- '~/projects/Ype_deJong_AAV/data/R1_controlsOnly.fastq.gz'
R2_out <- '~/projects/Ype_deJong_AAV/data/R2_controlsOnly.fastq.gz'

# 50mers from sacCer3 following the 6 AAV2 ITR remnants found in the Bushman positive controls.
seqs <- c('TATACGCGTTTGCACAATATACTCTATCTTATCCGTATCTATACGAGCCGG',
          'TGAGCTCAAGCACTTGCCTACTACTAAGTATGACGTCGTAATTGACCAGAA',
          'CTCAAGATAAAGCCATCATTCTGGGTGCCGAGGGTAACTTCCACGGGAGAA',
          'TTGGACAGTAGGTATGATTTTTTCAAGTTTTGGGAACCGGCAAAGAAATCG',
          'ATACATTCCGAGGGCGCCCGCACAAGGCCTATTATTAGAGGGACCTGTGTT',
          'GCTACACGTCCGATAACACGGACTCAATGACGTCCGGAGAAATCTCGGAGC')

I1 <- readFastq(I1_path)
R1 <- readFastq(R1_path)
R2 <- readFastq(R2_path)

ids <- unique(unlist(sapply(seqs, function(x){
         i <- vcountPattern(x, R2@sread, max.mismatch = 2)
         sub('\\s.+$', '', as.character(R2[i == 1]@id))
       })))

i <- sub('\\s.+$', '', as.character(R2@id))

writeFastq(I1[i %in% ids], I1_out)
writeFastq(R1[i %in% ids], R1_out)
writeFastq(R2[i %in% ids], R2_out)
