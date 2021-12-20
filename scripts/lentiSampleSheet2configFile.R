library(dplyr)
library(stringr)
library(RMySQL)
options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)

sampleSheet <- '/data/sequencing/Illumina-archive/211208_M03249_0229_000000000-G9RT3/SampleSheet.csv'
outputFile  <- '/home/everett/projects/AAVengeR2/scripts/211208_M03249_0229_000000000-G9RT3_sampleConfig.tsv'
refGenome.id <- 'hg38'
trialID <- 'Persaud_HIV'

f <- readLines(sampleSheet)
s <- read.csv(textConnection(f[(which(grepl('\\[metaData\\]', f))+1):length(f)]), header = TRUE)


invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
gtsp <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)

samples <- str_extract(s$SampleName, 'GTSP\\d+')
samples[is.na(samples)] <- s$SampleName[is.na(samples)]

subjects <- gtsp[match(samples, gtsp$SpecimenAccNum),]$Patient
subjects[is.na(subjects)] <- s$SampleName[is.na(subjects)]
subjects[grepl('PositiveControl', subjects)] <- 'PositiveControl'
subjects[grepl('NoTemplateControl', subjects)] <- 'NoTemplateControl'
subjects[grepl('UninfectedControl', subjects)] <- 'UninfectedControl'

replicates <- sub('\\-', '', str_extract(s$SampleName, '\\-\\d+'))

i <- which(is.na(subjects))
subjects[i] <- sub('\\-\\d+$', '', samples[i])

r <- tibble(trial = trialID,
            subject = subjects,
            sample = samples,
            replicate = replicates,
            adriftRead.linkerBarcode.start = 1,
            adriftRead.linkerBarcode.end = str_locate(s$linker, 'NNNNNNNNNNNN')[,1]-1,
            adriftRead.linker.seq = s$linker,
            index1.seq = s$barcode,
            refGenome.id = refGenome.id,
            adriftRead.linkerRandomID.start =  str_locate(s$linker, 'NNNNNNNNNNNN')[,1],
            adriftRead.linkerRandomID.end = str_locate(s$linker, 'NNNNNNNNNNNN')[,2],
            vectorFastaFile = ifelse(s$uniqueRegion == 'U5', 'u5_100_v1.fasta', 'u3_100_v1.fasta'),
            flags = ifelse(s$uniqueRegion == 'U5', 'HIV_u5', 'HIV_u3'))
            
write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

