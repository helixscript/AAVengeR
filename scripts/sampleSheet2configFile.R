library(dplyr)
library(stringr)
library(RMySQL)
options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)

sampleSheet <- '/data/sequencing/Illumina-archive/210316_MN01490_0011_A000H3FJL2/SampleSheet.csv'
outputFile  <- '/home/everett/projects/AAVengeRvsIntSiteCaller/210316_MN01490_0011_A000H3FJL2/sampleConfig.tsv'
refGenome.id <- 'hg38'

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

r <- tibble(subject = subjects,
            sample = samples,
            replicate = str_extract(s$alias, '(\\d+)$'),
            adriftRead.linkerBarcode.start = 1,
            adriftRead.linkerBarcode.end = str_locate(s$linkerSequence, 'NNNNNNNNNNNN')[,1]-1,
            adriftRead.linker.seq = s$linkerSequence,
            index1.seq = s$bcSeq,
            leaderSeqHMM = '/home/everett/projects/AAVengeR2/data/hmms/Bushman_SCID1.hmm',
            refGenome.id = refGenome.id,
            adriftRead.linkerRandomID.start =  str_locate(s$linkerSequence, 'NNNNNNNNNNNN')[,1],
            adriftRead.linkerRandomID.end = str_locate(s$linkerSequence, 'NNNNNNNNNNNN')[,2],
            vectorFastaFile = '/home/everett/data/BushmanGeneTherapy/vectorSequences/vector_SCID.fa',
            flags = 'HIV_u5')
            
write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

