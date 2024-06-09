# AAVengeR/barcodeAssocLinkers.R
# John K. Everett, Ph.D.
# 
# This script is a tool for working through barcoding and linker issues when
# the demultiplexing module does not work as expected. This module provides 
# useful summaries for technicians that can be used to check sample barcode/linker assignments.

suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(dplyr))

# Parse the config file from the command line.
configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - the configuration file does not exists.')

# Read config file.
opt <- yaml::read_yaml(configFile)

# Test for key items needed to run sanity tests.
if(! 'softwareDir' %in% names(opt)) stop('Error - the softwareDir parameter was not found in the configuration file.')
if(! dir.exists(opt$softwareDir)) stop(paste0('Error - the softwareDir directory (', opt$softwareDir, ') does not exist.'))

# Run config sanity tests.
source(file.path(opt$softwareDir, 'lib.R'))
optionsSanityCheck()

if(! dir.exists(file.path(opt$outputDir, opt$barcodeAssocLinkers_outputDir))) dir.create(file.path(opt$outputDir, opt$barcodeAssocLinkers_outputDir))

s <- readr::read_tsv(opt$barcodeAssocLinkers_sampleDataFile)

s$indexLinker <- paste0(s$index1Seq, substr(s$adriftReadLinkerSeq, 1, opt$barcodeAssocLinkers_adriftReadUniqueLinkerLength)) 

I1 <- readFastq(opt$barcodeAssocLinkers_index1ReadsFile)@sread
R1 <- readFastq(opt$barcodeAssocLinkers_adriftReadsFile)@sread

#  Add as an option.
### I1 <- reverseComplement(I1)

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

R1 <- as.character(subseq(R1, 1, opt$barcodeAssocLinkers_adriftReadUniqueLinkerLength))
tab <- sort(table(R1), decreasing = TRUE)
tab <- data.frame(sort(tab, decreasing = TRUE)[1:opt$barcodeAssocLinkers_nCodes])
  
r2 <- bind_rows(lapply(split(tab, 1:nrow(tab)), function(x){
        i <- which(R1 == x$R1)
        k <- data.frame(sort(table(I1[i]), decreasing = TRUE)[1:5]) %>% 
             mutate(uniqueLinker = x$R1, uniqueLinkerFreq = x$Freq, .before = 'Var1')
        names(k) <- c('uniqueLinker', 'uniqueLinkerFreq', 'barcode', 'barcodeFreq')
        k
      }))

openxlsx::write.xlsx(r, file.path(opt$outputDir, opt$barcodeAssocLinkers_outputDir, 'barcode2LinkerTable.xlsx'))
readr::write_tsv(r, file.path(opt$outputDir, opt$barcodeAssocLinkers_outputDir, 'barcode2LinkerTable.tsv'))

openxlsx::write.xlsx(r2, file.path(opt$outputDir, opt$barcodeAssocLinkers_outputDir, 'linker2barcodeTable.xlsx'))
readr::write_tsv(r2, file.path(opt$outputDir, opt$barcodeAssocLinkers_outputDir, 'linker2barcodeTable.tsv'))

