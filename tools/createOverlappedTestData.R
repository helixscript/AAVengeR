library(dplyr)
library(Biostrings)
dbPath    <- '/home/ubuntu/AAVengeR/data/referenceGenomes/blat/GCA_009914755.4.2bit'
outputDir <- '/home/ubuntu/projects/AAVengeR_synOverlappedReads/data'

remnants <- c('TCTGCGCGCTCGCTCGCTCACTGAG',
              'TCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAA',
              'TCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGA',
              'TCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTT',
              'TCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTG')
               
linker <- 'GAACGAGCACTAGTAAGCCCAAAAAAAAAAAACTCCGCTTAAGGGACT'
                               
g <- rtracklayer::import.2bit(dbPath)

targets <- c('chr18+58922000', 'chr11+69662500', 'chr8+22963550', 'chr17+61185000', 'chr10+63671820', 'chr19+6829550')

r <- bind_rows(lapply(targets, function(x){
       o <- unlist(strsplit(x, '[\\+\\-]'))
       s <- stringr::str_extract(x, '[\\+\\-]')
       seq <- subseq(g[names(g) == o[1]], as.integer(o[2]), as.integer(o[2]) + 500)

       bind_rows(lapply(1:3, function(n){
         A   <- subseq(seq, 10+n, 10+n + 150)
         AL1 <- subseq(seq, 9+n,   9+n + 150); AR1 <- subseq(seq, 11+n, 11+n + 150)
         AL2 <- subseq(seq, 8+n,   8+n + 150); AR2 <- subseq(seq, 12+n, 12+n + 150)
    
         B1 <- reverseComplement(subseq(seq, width(seq) - 150 - 0,  width(seq) - 0))
         B2 <- reverseComplement(subseq(seq, width(seq) - 150 - 10, width(seq) - 10))
         B3 <- reverseComplement(subseq(seq, width(seq) - 150 - 20, width(seq) - 20))
         B4 <- reverseComplement(subseq(seq, width(seq) - 150 - 30, width(seq) - 30))
         B5 <- reverseComplement(subseq(seq, width(seq) - 150 - 40, width(seq) - 40))
    
         R2 <- c(rep(paste0(remnants[n], as.character(AL2)), 1),
                 rep(paste0(remnants[n], as.character(AL1)), 3),
                 rep(paste0(remnants[n], as.character(A)),   5),
                 rep(paste0(remnants[n], as.character(AR1)), 3),
                 rep(paste0(remnants[n], as.character(AR2)), 1))
    
         R1 <- c(rep(paste0(linker, as.character(B1)), 1),
                 rep(paste0(linker, as.character(B2)), 3),
                 rep(paste0(linker, as.character(B3)), 5),
                 rep(paste0(linker, as.character(B4)), 3),
                 rep(paste0(linker, as.character(B5)), 1))
    
         tibble(readID = paste0(o[1], s, as.integer(o[2])+10+n, '_step', n ,'_read', 1:length(R1)), 
                R1 = R1, 
                R2 = R2)
       }))
}))


R1 <- ShortRead::ShortReadQ(sread = DNAStringSet(r$R1), 
                            id = BStringSet(r$readID), 
                            quality = BStringSet(rep(paste0(rep('?', 199), collapse = ''), nrow(r))))

R2 <- ShortRead::ShortReadQ(sread = DNAStringSet(r$R2), 
                            id = BStringSet(r$readID), 
                            quality = BStringSet(sapply(nchar(r$R2), function(x) paste0(rep('?', x), collapse = ''))))

I1 <- ShortRead::ShortReadQ(sread = DNAStringSet(rep('CATAGTTTAGCC', nrow(r))),
                            id = BStringSet(r$readID), 
                            quality = BStringSet(rep('????????????', nrow(r))))

writeFastq(R1, file.path(outputDir, 'testData_R1.fastq.gz'), compress = TRUE)
writeFastq(R2, file.path(outputDir, 'testData_R2.fastq.gz'), compress = TRUE)
writeFastq(I1, file.path(outputDir, 'testData_I1.fastq.gz'), compress = TRUE)

f <- tibble(trial = 'test', subject= 'test', sample = 'test', replicate = 1, 
            adriftReadLinkerSeq = 'GAACGAGCACTAGTAAGCCCNNNNNNNNNNNNCTCCGCTTAAGGGACT',
            index1Seq = 'CATAGTTTAGCC', refGenome = 'GCA_009914755.4', 
            vectorFastaFile = 'pAAV_GFP.fasta', anchorReadStartSeq ='TCTGCGCGCTCGCTC',
            flags = 'AAV')

readr::write_tsv(f, file.path(outputDir, 'testData_sampleData.tsv'))


#---

f <- f[grepl('chr10\\+6367', f$readID),]
f$step <- stringr::str_extract(f$readID, 'step\\d+')
f <- arrange(f, step)

Expect 1 - 39
>chr10+63671832_step2_read1
TCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAA CAGTGGCAGGTGGTGTGGGTTATAGGCGGCGGCGGCGGCGG
CTGCTGCTGCTGCTGAGCTGCTACCAGGCAGCCGGGGTCCCTCGCTGCCTCCACTGGCCCCTGGTCCGGTCACCCCAGCA
CTGGGGCCCCCCAGGGTAAAGTTACGGATT

>chr10+63671832_step2_read2
TCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAA AGTGGCAGGTGGTGTGGGTTATAGGCGGCGGCGGCGGCGGC
TGCTGCTGCTGCTGAGCTGCTACCAGGCAGCCGGGGTCCCTCGCTGCCTCCACTGGCCCCTGGTCCGGTCACCCCAGCAC
TGGGGCCCCCCAGGGTAAAGTTACGGATTG

# Already had wrong leader seq
> subset(a, readID == 'chr10+63671832_step2_read1')
uniqueSample                     readID       refGenome strand qSize qStart qEnd tName     tSize   tStart     tEnd queryPercentID tAlignmentWidth queryWidth alignmentPercentID percentQueryCoverage
1: test~test~test~1 chr10+63671832_step2_read1 GCA_009914755.4      +   151      0  151 chr10 134758134 63671828 63671978            100             151        152                100             100.6623
vectorFastaFile seqRunID flags                               leaderSeq
1:  pAAV_GFP.fasta      xxx   AAV TCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAA

> subset(a, readID == 'chr10+63671832_step2_read2')
uniqueSample                     readID       refGenome strand qSize qStart qEnd tName     tSize   tStart     tEnd queryPercentID tAlignmentWidth queryWidth alignmentPercentID percentQueryCoverage
1: test~test~test~1 chr10+63671832_step2_read2 GCA_009914755.4      +   145      0  145 chr10 134758134 63671835 63671979            100             145        146                100             100.6897
vectorFastaFile seqRunID flags                                     leaderSeq
1:  pAAV_GFP.fasta      xxx   AAV TCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAAGTGGC

