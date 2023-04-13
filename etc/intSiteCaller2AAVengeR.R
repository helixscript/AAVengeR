library(dplyr)
library(RMySQL)

inputFile <- '/home/ubuntu/projects/CD4_zeta/190910_M03249_0009_000000000-CLHJ9/data/metaData'
outputFile <- '/home/ubuntu/projects/CD4_zeta/190910_M03249_0009_000000000-CLHJ9/data/sampleData.tsv'

trial2HMM <- list()
trial2HMM[['UPENN_CART19_CLL']] = 'Bushman_CART19.hmm'
trial2HMM[['CART19_CLL']] = 'Bushman_CART19.hmm'
trial2HMM[['UPENN_CART19_ALL']] = 'Bushman_CART19.hmm'
trial2HMM[['SCID1_Paris_Cavazzana']] = 'Bushman_SCID1.hmm'
trial2HMM[['Zeta_HIV']] = 'Bushman_CD4zeta.hmm'
trial2HMM[['Control']] = 'Bushman_CART19.hmm'
trial2HMM[['WAS_Paris_Cavazzana']] = 'Bushman_CART19.hmm'
trial2HMM[['WAS_UK_Thrasher']] = 'Bushman_CART19.hmm'

m <- readr::read_csv(inputFile)

invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
gtsp <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)

samples <- sub('\\-\\d+$', '', m$alias)
reps <- as.integer(stringr::str_extract(m$alias, '\\d+$'))
flags <- 'IN_u5'

subjects <- gtsp[match(samples, gtsp$SpecimenAccNum),]$Patient
subjects <- ifelse(is.na(subjects), samples, subjects)

trial <- gtsp[match(samples, gtsp$SpecimenAccNum),]$Trial
trial <- ifelse(is.na(trial), 'Control', trial)

hmms <- sapply(trial, function(x) trial2HMM[[x]]) %>% unlist() %>% unname()

m$vectorSeq <- 'Bushman_CD4zeta.fasta'

r <- tibble(trial = trial,
            subject = subjects,
            sample = samples,
            replicate = reps,
            index1Seq = m$bcSeq,
            refGenome = ifelse(grepl('Positive', subject, ignore.case = TRUE), 'sacCer3', 'hg38'),
            adriftReadLinkerSeq = m$linkerSequence,
            vectorFastaFile = m$vectorSeq,
            leaderSeqHMM = hmms,
            flags = flags)

write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
