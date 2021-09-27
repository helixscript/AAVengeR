library(dplyr)
library(stringr)
library(RMySQL)
options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)

sampleSheet <- '/data/sequencing/Illumina-archive/210625_M03249_0193_000000000-DCLJC/SampleSheet.csv'
outputFile  <- '/home/everett/projects/Samelson-Jones_AAV/210625_M03249_0193_000000000-DCLJC/sampleConfig.tsv'

refGenome.id <- 'mm9'
anchorRead.seqFilter.file <- NA
ITRseqs <- paste0('Flip,TCTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCAAAGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAG;',
                  'Flop,TCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAG')

f <- readLines(sampleSheet)
s <- read.csv(textConnection(f[(which(grepl('\\[metaData\\]', f))+1):length(f)]), header = TRUE)


invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
gtsp <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)

samples <- sub('\\-\\d+$', '', s$alias)
subjects <- gtsp[match(samples, gtsp$SpecimenAccNum),]$Patient
i <- which(is.na(subjects))
subjects[i] <- sub('\\-\\d+$', '', samples[i])

r <- tibble(subject = subjects,
            sample = samples,
            replicate = str_extract(s$alias, '(\\d+)$'),
            adriftRead.linkerBarcode.start = 1,
            adriftRead.linkerBarcode.end = str_locate(s$linkerSequence, 'NNNNNNNNNNNN')[,1]-1,
            adriftRead.linker.seq = s$linkerSequence,
            index1.seq = s$bcSeq,
            anchorRead.identification = ITRseqs,
            anchorRead.seqFilter.file = anchorRead.seqFilter.file,
            refGenome.id = refGenome.id,
            adriftRead.linkerRandomID.start =  str_locate(s$linkerSequence, 'NNNNNNNNNNNN')[,1],
            adriftRead.linkerRandomID.end = str_locate(s$linkerSequence, 'NNNNNNNNNNNN')[,2])

write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

