library(dplyr)
library(RMySQL)

inputFile <- '/home/ubuntu/projects/Weinberger/data/metaData'
outputFile <- '/home/ubuntu/projects/Weinberger/data/sampleData.tsv'

trial2HMM <- list()
trial2HMM[['UPENN_CART19_CLL']] = 'Bushman_CART19.hmm'
trial2HMM[['UPENN_CART19_ALL']] = 'Bushman_CART19.hmm'
trial2HMM[['SCID1_Paris_Cavazzana']] = 'Bushman_SCID1.hmm'
trial2HMM[['Zeta_HIV']] = 'Bushman_CD4zeta.hmm'
trial2HMM[['Control']] = 'Bushman_CART19.hmm'
trial2HMM[['Weinberger_HIV']] = 'NA'

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

hmms <- ifelse(m$uniqueRegion == 'U5', 'HIV1_1-100_U5.hmm', 'HIV1_1-100_U3_RC.hmm') 
flags <- ifelse(m$uniqueRegion == 'U5', 'IN_u5', 'IN_u3') 

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
