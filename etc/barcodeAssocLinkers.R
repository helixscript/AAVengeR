library(ShortRead)
library(dplyr)

nCodes <- 50
adriftLinkerSeq_sampleLength <- 20

I1_path    <- 'Undetermined_S0_I1_001.fastq.gz'
R1_path    <- 'Undetermined_S0_R1_001.fastq.gz'
sampleData <- 'sampleData.tsv'

s <- readr::read_tsv(sampleData)
s$indexLinker <- paste0(s$index1Seq, substr(s$adriftReadLinkerSeq, 1, adriftLinkerSeq_sampleLength)) 

I1 <- readFastq(I1_path)@sread
R1 <- readFastq(R1_path)@sread

tab <- table(as.character(I1))
tab <- tab[! grepl('GGGGG', names(tab))]  # Exclude poly-G (MiniSeq no signal for phi-X)
tab <- data.frame(sort(tab, decreasing = TRUE)[1:nCodes])

r <- bind_rows(lapply(split(tab, 1:nrow(tab)), function(x){
       i <- which(as.character(I1) == x$Var1)
       k <- sort(table(subseq(R1[i], 1, 20)), decreasing = TRUE)[1]
       x$topAssocLinker <- names(k)
       x$percentTopAssocLinker <- sprintf("%.1f%%", (k/length(i))*100)
       
       x$barcodeInSampleData <- x$Var1 %in% s$index1Seq
       
       # Test if the linker sequence associated with this barcode in the metadata is the same as the top linker found.
       testSeq <- subset(s, index1Seq == x$Var1)$indexLinker
       if(length(testSeq) == 1){
         x$barcodeLinkerPairInSampleData <- paste0(x$Var1, x$topAssocLinker) == testSeq
       } else {
         x$barcodeLinkerPairInSampleData <- NA
       }
       
       names(x) <- c('barcode', 'nObs', 'mostAssocLinker', 'percentAssocLinkers', 'barcodeInSampleData', 'barcodeLinkerPairInSampleData')
       x
     }))

openxlsx::write.xlsx(r, 'barcodeAssocLinkers.xlsx')