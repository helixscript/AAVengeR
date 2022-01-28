library(dplyr)
library(stringr)
library(RMySQL)
options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)

sampleSheet <- '/data/sequencingData/220124_MN01490_0059_A000H3M2VG/SampleSheet.csv'
outputFile  <- '/data/project/Persaud_HIV/x'
refGenome.id <- 'hg38'
trialID <- 'Persaud_HIV'

f <- readLines(sampleSheet)
s <- read.csv(textConnection(f[(which(grepl('\\[metaData\\]', f))+1):length(f)]), header = TRUE)


invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
gtsp <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)

samples <- str_extract(s$SampleName, 'GTSP\\d+')
samples[is.na(samples)] <- s$SampleName[is.na(samples)]

subjects <- gtsp[match(samples, gtsp$SpecimenAccNum),]$Patient
subjects[is.na(subjects)] <- s$SampleName[is.na(subjects)]
subjects[grepl('PositiveControl', subjects)] <- 'PositiveControl'
subjects[grepl('NoTemplateControl', subjects)] <- 'NoTemplateControl'
subjects[grepl('UninfectedControl', subjects)] <- 'UninfectedControl'

replicates <- sub('\\-', '', str_extract(s$SampleName, '\\-\\d+'))

i <- which(is.na(subjects))
subjects[i] <- sub('\\-\\d+$', '', samples[i])

r <- tibble(trial = trialID,
            subject = subjects,
            sample = samples,
            replicate = replicates,
            adriftRead.linkerBarcode.start = 1,
            adriftRead.linkerBarcode.end = str_locate(s$linker, 'NNNNNNNNNNNN')[,1]-1,
            adriftRead.linker.seq = s$linker,
            index1.seq = s$barcode,
            refGenome.id = refGenome.id,
            adriftRead.linkerRandomID.start =  str_locate(s$linker, 'NNNNNNNNNNNN')[,1],
            adriftRead.linkerRandomID.end = str_locate(s$linker, 'NNNNNNNNNNNN')[,2],
            vectorFastaFile = 'HIV-1_reference.fasta',
            leaderSeqHMM = ifelse(s$uniqueRegion == 'U5', 'u5_100.hmm', 'u3_100_rc.hmm'),
            flags = ifelse(s$uniqueRegion == 'U5', 'HIV_u5', 'HIV_u3'))
            
write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

