library(dplyr)
library(RMySQL)

inputFile <- '/home/ubuntu/projects/AAVengeR_230306_inward/data/metaData'
outputFile <- '/home/ubuntu/projects/AAVengeR_230306_inward/data/sampleData.tsv'

m <- readr::read_csv(inputFile)

invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
gtsp <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)

samples <- m$sample
reps <- m$replicate
flags <- 'AAV'

subjects <- gtsp[match(samples, gtsp$SpecimenAccNum),]$Patient
subjects <- ifelse(is.na(subjects), samples, subjects)

trial <- gtsp[match(samples, gtsp$SpecimenAccNum),]$Trial
trial <- ifelse(is.na(trial), 'Control', trial)


r <- tibble(trial = trial,
            subject = subjects,
            sample = samples,
            replicate = reps,
            index1Seq = m$index1Seq,
            refGenome = ifelse(grepl('DogHemo', trial, ignore.case = TRUE), 'canFam3', 'hg38'),
            adriftReadLinkerSeq = m$adriftReadLinkerSeq,
            vectorFastaFile = m$vectorFastaFile,
            flags = flags)

write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
