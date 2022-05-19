library(dplyr)
library(stringr)
library(RMySQL)

sampleSheet <- '~/zeta/220405_MN01490_0071_A000H3VTVJ/SampleSheet.csv'
outputFile  <- '~/zeta/sampleData.tsv'
refGenome.id <- 'hg38'
trialID <- 'cd4Zeta'

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
            index1.seq = s$bcSeq,
            refGenome.id = refGenome.id,
            adriftRead.linker.seq = s$linkerSequence,
            vectorFastaFile = s$vectorSeq,
            leaderSeqHMM = 'cd4Zeta.hmm',
            flags = 'HIV_u5')
            
write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)


#adriftRead.linkerBarcode.start = 1,
#adriftRead.linkerBarcode.end = str_locate(s$linkerSequence, 'NNNNNNNNNNNN')[,1]-1,
#adriftRead.linkerRandomID.start =  str_locate(s$linkerSequence, 'NNNNNNNNNNNN')[,1],
#adriftRead.linkerRandomID.end = str_locate(s$linkerSequence, 'NNNNNNNNNNNN')[,2],


