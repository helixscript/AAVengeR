library(dplyr)
library(stringr)
library(RMySQL)

sampleSheet <- '/home/ubuntu/projects/Persaud/201123_M03249_0121_000000000-D6KV5/data/SampleSheet.csv'
outputFile  <- '/home/ubuntu/projects/Persaud/201123_M03249_0121_000000000-D6KV5/data/sampleData.tsv'
refGenome.id <- 'hg38'
trialID <- 'Persaud'

f <- readLines(sampleSheet)
s <- read.csv(textConnection(f[(which(grepl('\\[metaData\\]', f))+1):length(f)]), header = TRUE)


invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
gtsp <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)

samples <- sub('\\-\\d+$', '', s$SampleName)
reps <- as.integer(sub('\\-', '', stringr::str_extract(samples, '\\-\\d+')))
flags <- ifelse(grepl('u5', samples, ignore.case = TRUE), 'IN_u5', 'IN_u3')
samples <- stringr::str_extract(samples, 'GTSP\\d+')
subjects <- gtsp[match(samples, gtsp$SpecimenAccNum),]$Patient

r <- tibble(trial = trialID,
            subject = subjects,
            sample = samples,
            replicate = reps,
            index1.seq = s$barcode,
            refGenome.id = refGenome.id,
            adriftRead.linker.seq = s$linker,
            leaderSeqHMM = ifelse(flags == 'IN_u5', 'u5_100.hmm', 'u3_100_rc.hmm'), 
            vectorFastaFile = 'HIV-1_reference.fasta',
            flags = flags)
            
write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)





