library(dplyr)
library(stringr)
library(RMySQL)

sampleSheet <- '/home/ubuntu/projects/Lichterfeld/data/SampleSheet.csv'
outputFile  <- '/home/ubuntu/projects/Lichterfeld/data/sampleData.tsv'
refGenome.id <- 'hg38'
trialID <- 'Lichterfeld'

f <- readLines(sampleSheet)
s <- read.csv(textConnection(f[(which(grepl('\\[metaData\\]', f))+1):length(f)]), header = TRUE)

if('Sampleme' %in% names(s)) s$SampleName <- s$Sampleme


invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
gtsp <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)

# s$SampleName <- sub('20221121\\-', '20221121_', s$SampleName)

samples <- sub('\\-\\d+$', '', s$SampleName)
reps <- as.integer(sub('\\-', '', stringr::str_extract(samples, '\\-\\d+')))
flags <- ifelse(grepl('u5', samples, ignore.case = TRUE), 'IN_u5', 'IN_u3')

samples <- stringr::str_extract(samples, 'GTSP\\d+')
samples <- ifelse(is.na(samples), sub('\\-\\S+$', '', s$SampleName), samples)

subjects <- gtsp[match(samples, gtsp$SpecimenAccNum),]$Patient
subjects <- ifelse(is.na(subjects), samples, subjects)

r <- tibble(trial = trialID,
            subject = subjects,
            sample = samples,
            replicate = reps,
            index1.seq = s$barcode,
            refGenome.id = refGenome.id,
            adriftRead.linker.seq = s$linker,
            leaderSeqHMM = ifelse(flags == 'IN_u5', 'HIV1_1-100_U5.hmm', 'HIV1_1-100_U3_RC.hmm'), 
            vectorFastaFile = 'HXB2.fasta',
            flags = flags)
            
write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
