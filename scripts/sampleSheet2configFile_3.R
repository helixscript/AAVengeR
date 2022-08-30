library(dplyr)
library(stringr)
library(RMySQL)

sampleSheet <- '/home/ubuntu/projects/ViiV_HIV/data/220510_MN01490_0079_A000H3W5LH/SampleSheet.csv'
outputFile  <- '/home/ubuntu/projects/ViiV_HIV/data/220510_MN01490_0079_A000H3W5LH/sampleData.tsv'
refGenome.id <- 'hg38'
trialID <- 'ViiV'

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
            leaderSeqHMM = 'vector_HIV_89.6.hmm', 
            vectorFastaFile = 'vector_HIV_89.6.fa',
            flags = 'IN_u5')
            
write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)





