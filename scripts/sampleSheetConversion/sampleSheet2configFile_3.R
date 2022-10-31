library(dplyr)
library(stringr)
library(RMySQL)

sampleSheet <- '/home/ubuntu/projects/CD4zeta/data/SampleSheet.csv'
outputFile  <- '/home/ubuntu/projects/CD4zeta/samplData.tsv'
refGenome.id <- 'hg38'
trialID <- 'CD4zeta'

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
            leaderSeqHMM = 'cd4Zeta.hmm', 
            vectorFastaFile = 'vector_ZETA.fa',
            flags = 'IN_u5')
            
write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)





