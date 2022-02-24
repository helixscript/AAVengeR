library(dplyr)
library(stringr)
library(readr)
options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)

outputFile <- '/data/project/Spark/previousSampleData.tsv'

x <- read_csv('/data/project/Spark/previousSampleData.csv')
d <- read_csv('/data/project/Spark/previousSampleDetails.csv')
d$subject <- gsub(' ', '_', d$subject)
x <- left_join(x, select(d, sample, subject), by = 'sample')

o <- tibble(trial = 'Spark',
            subject = x$subject,
            sample = x$sample,
            replicate = 1,
            adriftRead.linkerBarcode.start = 1,
            adriftRead.linkerBarcode.end = str_locate(x$linker, 'NNNNNNNNNNNN')[,1]-1,
            adriftRead.linker.seq = x$linker,
            index1.seq = x$barcode,
            refGenome.id = x$refGenome,
            adriftRead.linkerRandomID.start =  str_locate(x$linker, 'NNNNNNNNNNNN')[,1],
            adriftRead.linkerRandomID.end = str_locate(x$linker, 'NNNNNNNNNNNN')[,2],
            vectorFastaFile = 'Spark_2020.fasta',
            anchorRead.startSeq = 'TCCCTCTCTGCGCGC',
            flags = 'AAV')
            
write.table(o, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

