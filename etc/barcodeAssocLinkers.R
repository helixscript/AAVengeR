library(ShortRead)
library(dplyr)

nCodes <- 50

I1_path <- 'Undetermined_S0_I1_001.fastq.gz'
R1_path <- 'Undetermined_S0_R1_001.fastq.gz'

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
       names(x) <- c('barcode', 'nObs', 'mostAssocLinker', 'percentAssocLinkers')
       x
      }))

openxlsx::write.xlsx(r, 'barcodeAssocLinkers.xlsx')