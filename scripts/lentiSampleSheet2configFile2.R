library(dplyr)
library(stringr)
library(RMySQL)
options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)

sampleSheet <- '/data/sequencing/Illumina-archive/211221_M03249_0230_000000000-K4WDM/SampleSheet.csv'
outputFile  <- '/home/everett/projects/AAV_controls/211221_M03249_0230_000000000-K4WDM/sampleConfig.tsv'
refGenome.id <- 'sacCer3'
trialID <- 'AAV_pos_control'

f <- readLines(sampleSheet)
s <- read.csv(textConnection(f[(which(grepl('\\[metaData\\]', f))+1):length(f)]), header = TRUE)


invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
gtsp <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)

samples <- str_extract(s$SampleName, 'GTSP\\d+')
samples[is.na(samples)] <- s$SampleName[is.na(samples)]

samples <- gsub('\\s', '_', s$alias)

subjects <- ifelse(grepl('Dilution', samples), 'Dilution', samples)

subjects[grepl('PositiveControl', subjects)] <- 'PositiveControl'
subjects[grepl('NoTemplateControl', subjects)] <- 'NoTemplateControl'
subjects[grepl('UninfectedControl', subjects)] <- 'UninfectedControl'

replicates <- sub('\\-', '', str_extract(s$alias, '\\-\\d$'))
samples <- sub('\\-\\d$', '', samples)


r <- tibble(trial = trialID,
            subject = subjects,
            sample = samples,
            replicate = replicates,
            adriftRead.linkerBarcode.start = 1,
            adriftRead.linkerBarcode.end = str_locate(s$linkerSequence, 'NNNNNNNNNNNN')[,1]-1,
            adriftRead.linker.seq = s$linkerSequence,
            index1.seq = s$bcSeq,
            refGenome.id = refGenome.id,
            adriftRead.linkerRandomID.start =  str_locate(s$linkerSequence, 'NNNNNNNNNNNN')[,1],
            adriftRead.linkerRandomID.end = str_locate(s$linkerSequence, 'NNNNNNNNNNNN')[,2],
            vectorFastaFile = 'AAV_pos_control.fasta',
            flags = 'AAV')
            
write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

