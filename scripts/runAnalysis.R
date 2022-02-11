library(ShortRead)
library(dplyr)
library(readr)
library(Biostrings)

n <- 500
BushmanBarCodes <- readr::read_tsv('/home/ubuntu/software/AAVengeR/data/misc/BushmanGroup_barCodes')

sampleConfig <- read_tsv('/data/project/Encoded/220208_M03249_0241_000000000-K4YRL/sampleConfig.tsv')
#sampleConfig$index1.seq <- as.character(reverseComplement(DNAStringSet(sampleConfig$index1.seq)))

I1 <- '/data/sequencingData/220208_M03249_0241_000000000-K4YRL/Undetermined_S0_I1_001.fastq.gz'
R1 <- '/data/sequencingData/220208_M03249_0241_000000000-K4YRL/Undetermined_S0_R1_001.fastq.gz'

I1 <- readFastq(I1)
R1 <- readFastq(R1)

I1.ids <- as.character(I1@id)
R1.ids <- as.character(R1@id)

I1 <- I1@sread
R1 <- R1@sread

names(I1) <- I1.ids
names(R1) <- R1.ids

names(I1) <- sub('\\s+.+$', '', names(I1))
names(R1) <- sub('\\s+.+$', '', names(R1))

# Identify the to n bar code sequences.
i <- sort(table(as.character(I1)), decreasing = TRUE)[1:n]

I1.s <- I1[as.character(I1) %in% names(i)]
R1.s <- R1[names(R1) %in% names(I1.s)]

barCodes <- unique(as.character(I1.s))

I1.s.char <- as.character(I1.s)

r <- bind_rows(lapply(barCodes, function(x){
  o <- I1.s[I1.s.char == x]
  a <- R1[names(R1) %in% names(o)]
  b <- sort(table(as.character(subseq(a, 1, 20))), decreasing = TRUE)[1:3]
  seqs <- paste0(names(b), collapse = ',')
  p <- paste0(sprintf("%.1f%%", (b/length(o))*100), collapse = ', ')
  tibble(barcode = x, reads = length(o), inSampleData = ifelse(x %in% sampleConfig$index1.seq, 'yes', 'no'),
         BushmanCode = ifelse(x %in% BushmanBarCodes$Sequence, 'yes', 'no'),
         BushmanCode_RC = ifelse(x %in% BushmanBarCodes$RevComp, 'yes', 'no'), topLinkerPercents = p, topLinkers = seqs )
})) %>% arrange(desc(reads))

openxlsx::write.xlsx(r, 'runAnalysis.xlsx')

