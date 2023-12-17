library(dplyr)
library(Biostrings)
dbPath    <- '/home/ubuntu/AAVengeR/data/referenceGenomes/blat/GCA_009914755.4.2bit'
outputDir <- '/home/ubuntu/projects/AAVengeR_synOverlappedReads/data'

nSites <- 1000

ltrBit <- 'GAAAATCTCTAGCA'
linker <- 'GAACGAGCACTAGTAAGCCCAAAAAAAAAAAACTCCGCTTAAGGGACT'
                               
g <- rtracklayer::import.2bit(dbPath)
g <- g[! names(g) %in% c('chrM', 'chrY')]

set.seed(1)
targets <- bind_rows(lapply(sample(names(g), nSites, replace = TRUE), function(x){
            tibble(chromosome = x,
                   strand = sample(c('+', '-'), 1),
                   position = sample((1+10000):(width(g[names(g) == x])-10000), 1))
           }))

r <- bind_rows(lapply(split(targets, 1:nrow(targets)), function(x){
       seq <- subseq(g[names(g) == x$chromosome], x$position, x$position+1000)

       bind_rows(lapply(c(0, 10, 20, 30, 40), function(f){
           if(x$strand == '+'){
             A <- subseq(seq,  100, 150 + 100)
             AL1 <- subseq(seq, 99, 150 + 99);  AR1 <- subseq(seq, 101, 150 + 101)
             AL2 <- subseq(seq, 98, 150 + 98);  AR2 <- subseq(seq, 102, 150 + 102)
             
             B <- subseq(seq, 495 + f - 150, 495 + f)
             BL1 <- subseq(seq, 495 + f - 151, 494 + f);  BR1 <- subseq(seq, 496 + f - 150, 496 + f)
             BL2 <- subseq(seq, 495 + f - 152, 493 + f);  BR2 <- subseq(seq, 497 + f - 150, 497 + f)
             
             R2 <- c(rep(paste0(ltrBit, as.character(AL2)), 1),
                     rep(paste0(ltrBit, as.character(AL1)), 3),
                     rep(paste0(ltrBit, as.character(A)),   5),
                     rep(paste0(ltrBit, as.character(AR1)), 3),
                     rep(paste0(ltrBit, as.character(AR2)), 1))
             
             R1 <- c(rep(paste0(linker, as.character(BL2)), 1),
                     rep(paste0(linker, as.character(BL1)), 3),
                     rep(paste0(linker, as.character(B)),   5),
                     rep(paste0(linker, as.character(BR1)), 3),
                     rep(paste0(linker, as.character(BR2)), 1))
             
             R1 <- as.character(reverseComplement(DNAStringSet(R1)))
             
           } else {
             A <- subseq(seq,  100 - f, 150 + 100 - f)
             AL1 <- subseq(seq, 99 - f, 150 + 99 - f);  AR1 <- subseq(seq, 101 - f, 150 + 101 - f)
             AL2 <- subseq(seq, 98 - f, 150 + 98 - f);  AR2 <- subseq(seq, 102 - f, 150 + 102 - f)
             
             B <- subseq(seq, 495 - 150, 495)
             BL1 <- subseq(seq, 495 - 151, 494);  BR1 <- subseq(seq, 496 - 150, 496)
             BL2 <- subseq(seq, 495 - 152, 493);  BR2 <- subseq(seq, 497 - 150, 497)
             
             R1 <- c(rep(paste0(linker, as.character(AL2)), 1),
                     rep(paste0(linker, as.character(AL1)), 3),
                     rep(paste0(linker, as.character(A)),   5),
                     rep(paste0(linker, as.character(AR1)), 3),
                     rep(paste0(linker, as.character(AR2)), 1))
             
             R2 <- c(rep(paste0(ltrBit, as.character(BL2)), 1),
                     rep(paste0(ltrBit, as.character(BL1)), 3),
                     rep(paste0(ltrBit, as.character(B)),   5),
                     rep(paste0(ltrBit, as.character(BR1)), 3),
                     rep(paste0(ltrBit, as.character(BR2)), 1))
             
             R2 <- as.character(reverseComplement(DNAStringSet(R2)))
           }
    
    
           tibble(readID = paste0(x$chromosome, x$strand, x$position, '_step_', f, '_read_', 1:length(R1)), 
                  R1 = R1, R2 = R2)
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

