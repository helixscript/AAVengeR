library(dplyr)
library(stringr)
library(RMySQL)

sampleSheet <- '/home/ubuntu/projects/Encoded/20220817/data/SampleSheet.csv'
outputFile  <- '/home/ubuntu/projects/Encoded/20220817/sampleData.tsv'
refGenome.id <- 'macFas6_alu_before_MALAT1_site'
trialID <- 'Encoded'

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
            vectorFastaFile = 'Encoded.ALUs_masked.fasta',
            flags = 'AAV')
            
write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)





