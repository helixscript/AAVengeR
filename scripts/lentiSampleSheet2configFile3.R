library(dplyr)
library(stringr)
library(RMySQL)
options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)

sampleSheet <- '/data/sequencingData/220318_MN01490_0067_A000H3M32V_cf/SampleSheet.csv'
outputFile  <- '/data/sequencingData/220318_MN01490_0067_A000H3M32V_cf/sampleData.tsv'
refGenome.id <- 'hg38'
trialID <- 'Gill_cellFree'

f <- readLines(sampleSheet)
s <- read.csv(textConnection(f[(which(grepl('\\[metaData\\]', f))+1):length(f)]), header = TRUE)
#s <- read.csv(sampleSheet, header = TRUE)

invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
gtsp <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)

samples <- str_extract(s$alias, 'GTSP\\d+')
samples[is.na(samples)] <- s$alias[is.na(samples)]

subjects <- gtsp[match(samples, gtsp$SpecimenAccNum),]$Patient
subjects[is.na(subjects)] <- s$alias[is.na(subjects)]
subjects[grepl('PositiveControl', subjects)] <- 'PositiveControl'
subjects[grepl('NoTemplateControl', subjects)] <- 'NoTemplateControl'
subjects[grepl('UninfectedControl', subjects)] <- 'UninfectedControl'
subjects[grepl('Extract', subjects)] <- 'ExtractionControl'

replicates <- sub('\\-', '', str_extract(s$alias, '\\-\\d+'))

i <- which(is.na(subjects))
subjects[i] <- sub('\\-\\d+$', '', samples[i])

r <- tibble(trial = trialID,
            subject = subjects,
            sample = samples,
            replicate = replicates,
            adriftRead.linkerBarcode.start = 1,
            adriftRead.linkerBarcode.end = str_locate(s$linker, 'NNNNNNNNNNNN')[,1]-1,
            adriftRead.linker.seq = s$linker,
            index1.seq = s$bcSeq,
            refGenome.id = refGenome.id,
            adriftRead.linkerRandomID.start =  str_locate(s$linker, 'NNNNNNNNNNNN')[,1],
            adriftRead.linkerRandomID.end = str_locate(s$linker, 'NNNNNNNNNNNN')[,2],
            vectorFastaFile = 'huCART19BbzC2137.fa',
            leaderSeqHMM = 'lenti_GAAAATCTCTAGCA.hmm',
            flags = 'HIV_u5')
            
write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

