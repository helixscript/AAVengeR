library(dplyr)
library(stringr)
library(RMySQL)
options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)

sampleSheet <- '/data/sequencing/Illumina-archive/201221_M02973_0295_000000000-JDY29/SampleSheet.csv'
outputFile  <- '/home/everett/projects/AAVengeR2/scripts/201221_M02973_0295_000000000-JDY29_sampleConfig.tsv'
refGenome.id <- 'hg38'
trialID <- 'CSL'

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

r <- tibble(trial = trialID,
            subject = subjects,
            sample = samples,
            replicate = str_extract(s$alias, '(\\d+)$'),
            adriftRead.linkerBarcode.start = 1,
            adriftRead.linkerBarcode.end = str_locate(s$linkerSequence, 'NNNNNNNNNNNN')[,1]-1,
            adriftRead.linker.seq = s$linkerSequence,
            index1.seq = s$bcSeq,
            refGenome.id = refGenome.id,
            adriftRead.linkerRandomID.start =  str_locate(s$linkerSequence, 'NNNNNNNNNNNN')[,1],
            adriftRead.linkerRandomID.end = str_locate(s$linkerSequence, 'NNNNNNNNNNNN')[,2],
            vectorFastaFile = s$vectorSeq,
            flags = 'AAV')
            
write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

