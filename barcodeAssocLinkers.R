library(ShortRead)
library(dplyr)

# Read in configuration file.
configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)

if(! dir.exists(file.path(opt$outputDir, opt$barcodeAssocLinkers_outputDir))) dir.create(file.path(opt$outputDir, opt$barcodeAssocLinkers_outputDir))

s <- readr::read_tsv(opt$barcodeAssocLinkers_sampleDataFile)

s$indexLinker <- paste0(s$index1Seq, substr(s$adriftReadLinkerSeq, 1, opt$barcodeAssocLinkers_adriftReadUniqueLinkerLength)) 

I1 <- readFastq(opt$barcodeAssocLinkers_index1ReadsFile)@sread
R1 <- readFastq(opt$barcodeAssocLinkers_adriftReadsFile)@sread

I1 <- as.character(I1)

tab <- table(I1)
if(opt$barcodeAssocLinkers_excludePolyG_index1BarCodes) tab <- tab[! grepl('GGGGG', names(tab))]  # Exclude poly-G (MiniSeq no signal for phi-X)
tab <- data.frame(sort(tab, decreasing = TRUE)[1:opt$barcodeAssocLinkers_nCodes])

r <- bind_rows(lapply(split(tab, 1:nrow(tab)), function(x){
       i <- which(I1 == x$I1)
       k <- sort(table(subseq(R1[i], 1, opt$barcodeAssocLinkers_adriftReadUniqueLinkerLength)), decreasing = TRUE)[1]
       x$topAssocLinker <- names(k)
       x$percentTopAssocLinker <- sprintf("%.1f%%", (k/length(i))*100)
      
       x$barcodeInSampleData <- x$I1 %in% s$index1Seq
       
       # Test if the linker sequence associated with this barcode in the metadata is the same as the top linker found.
       testSeq <- subset(s, index1Seq == x$I1)$indexLinker
       if(length(testSeq) == 1){
         x$barcodeLinkerPairInSampleData <- paste0(x$I1, x$topAssocLinker) == testSeq
       } else {
         x$barcodeLinkerPairInSampleData <- NA
       }
       
       names(x) <- c('barcode', 'nObs', 'mostAssocLinker', 'percentAssocLinkers', 'barcodeInSampleData', 'barcodeLinkerPairInSampleData')
       x
     }))

openxlsx::write.xlsx(r, file.path(opt$outputDir, opt$barcodeAssocLinkers_outputDir, 'barcodeTable.xlsx'))
readr::write_tsv(r, file.path(opt$outputDir, opt$barcodeAssocLinkers_outputDir, 'barcodeTable.tsv'))