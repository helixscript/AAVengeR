library(dplyr)
library(RMySQL)

inputFile <- '/home/ubuntu/projects/Jones_CD8/230511_MN01490_0132_A000H5JTML/data/metaData'
outputFile <- '/home/ubuntu/projects/Jones_CD8/230511_MN01490_0132_A000H5JTML/data/sampleData.tsv'

m <- readr::read_csv(inputFile)

invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
gtsp <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)

samples <- sub('\\-\\d+.+$', '', m$SampleName)
reps <- as.integer(unlist(lapply(stringr::str_extract_all(m$SampleName, '\\d+'), '[', 2)))


subjects <- gtsp[match(samples, gtsp$SpecimenAccNum),]$Patient
subjects <- ifelse(is.na(subjects), samples, subjects)

trial <- gtsp[match(samples, gtsp$SpecimenAccNum),]$Trial
trial <- ifelse(is.na(trial), 'Control', trial)

hmms <- ifelse(grepl('u5', m$SampleName, ignore.case = TRUE), 'HIV1_1-100_U5.hmm', 'HIV1_1-100_U3_RC.hmm') 
flags <- ifelse(grepl('u5', m$SampleName, ignore.case = TRUE), 'IN_u5', 'IN_u3') 

m$vectorSeq <- 'HXB2.fasta'

r <- tibble(trial = trial,
            subject = subjects,
            sample = samples,
            replicate = reps,
            index1Seq = m$barcode,
            refGenome = ifelse(grepl('Positive', subject, ignore.case = TRUE), 'sacCer3', 'hg38'),
            adriftReadLinkerSeq = m$linker,
            vectorFastaFile = m$vectorSeq,
            leaderSeqHMM = hmms,
            flags = flags)

write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
