library(dplyr)
library(stringr)

sampleSheet <- '/home/ubuntu/projects/AAV_dev/20220819/data/SampleSheet.csv'
outputFile  <- '/home/ubuntu/projects/AAV_dev/20220819/sampleData.tsv'
trialID <- 'AAV_dev'

f <- readLines(sampleSheet)
s <- read.csv(textConnection(f[(which(grepl('\\[metaData\\]', f))+1):length(f)]), header = TRUE)

subjects <- ifelse(grepl('Control', s$alias, ignore.case = TRUE), 'Positive_control', 'Negative_control')
samples <- sub('\\-\\d+$', '', s$alias)
replicates <- str_extract(s$alias, '(\\d+)$')

r <- tibble(trial = trialID,
            subject = subjects,
            sample = samples,
            replicate = replicates,
            index1.seq = s$bcSeq,
            refGenome.id = s$refGenome,
            adriftRead.linker.seq = s$linkerSequence,
            vectorFastaFile = s$vectorSeq,
            flags = 'AAV')
            
write.table(r, file = outputFile, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)





