# John K. Everett, PhD
# Jan. 2024

library(dplyr)
library(Biostrings)
library(parallel)

# Selecting either U3 or U5 will change the output labeling, LTR bit, 
# and HMM in resulting sampleData file.

direction <- 'U3' # Select U3 or U5

# Data paths.
dbPath    <- '../../../data/referenceGenomes/blat/hg38.2bit'
outputDir <- './'


# Configuration.
nCPUs           <- 5
readWidth       <- 150
fragWidthStart  <- 500
nFragsPerSite   <- 5
nReadsPerFrag   <- 30

set.seed(sum(utf8ToInt(direction)) + 1)

# Target file selection determines the number of sites and their specifications.

if(direction == 'U3'){
  filePrefix <- 'U3_'
  hmmName <- 'HXB2_U3_RC.hmm'
  sampleFlag <- 'IN_u3'
  ltrBit <- 'AATTAGCCCTTCCA'
  targets <- readr::read_tsv('U3_targets.tsv')
} else {
  filePrefix <- 'U5_'
  hmmName <- 'HXB2_U5.hmm'
  sampleFlag <- 'IN_u5'
  ltrBit <- 'GAAAATCTCTAGCA'
  targets <- readr::read_tsv('U5_targets.tsv')
}
                               
# Read in genome, exclude mitochondria and Y chromosomes.
g <- rtracklayer::import.2bit(dbPath)
g <- g[! names(g) %in% c('chrM', 'chrY')]

cluster <- makeCluster(nCPUs)
clusterSetRNGStream(cluster, 1)
clusterExport(cluster, c('ltrBit', 'g', 'readWidth'))

# Create a series of values around a mean beneath a tight Gaussian curve. 
createPositions <- function(pos, n, d = 3){
  
  # Make half of the positions exactly the requested length.
  h <- rep(pos, floor(n/2))
  
  o <- round(rnorm((n*10), pos, sd = 0.75))  # Create 10x more than needed
  o <- o[abs(o - pos) <= d]                  # Select values within a distance of d from pos
  c(o[1:(n - length(h))], h)                 # return a selection of positions
}

# Function to create random DNA sequences.
randomDNAseq <- function(x, n = 12){
  paste(paste0(stringi::stri_rand_strings(n, 1, '[ATCG]')), collapse = '')
}


# Create fragment boundaries for selected targets.
r <- bind_rows(lapply(split(targets, 1:nrow(targets)), function(x){
  
       if(x$strand == '+'){
         d <- tibble(target = paste0(x$chromosome, x$strand, x$position),
                     startPos = x$position,
                     fragID = sapply(0:(nFragsPerSite-1), function(n){ paste0(x$chromosome, x$strand, x$position, '_fragment', n) }),
                     endPos = sapply(0:(nFragsPerSite-1), function(n){ x$position + fragWidthStart + n*10}))
       } else {
         d <- tibble(target = paste0(x$chromosome, x$strand, x$position),
                     startPos = sapply(0:(nFragsPerSite-1), function(n){ x$position - (fragWidthStart + n*10)}),
                     fragID = sapply(0:(nFragsPerSite-1), function(n){ paste0(x$chromosome, x$strand, x$position, '_fragment', n) }),
                     endPos = x$position)
       } 
       
       bind_rows(lapply(split(d, 1:nrow(d)), function(k){
         tibble(target = k$target,
                chromosome = x$chromosome,
                strand = x$strand,
                fragID = k$fragID,
                readStartPos = as.character(sapply(k$startPos, createPositions, n = nReadsPerFrag, d = 3)),
                readEndPos = as.character(sapply(k$endPos, createPositions, n = nReadsPerFrag, d = 3)))
       })) %>% mutate(readID = paste0(target, '_', fragID, '_read', 1:n()))
     }))

r$readStartPos <- as.integer(r$readStartPos)
r$readEndPos <- as.integer(r$readEndPos)

# Add UMI sequences to the proto-reads.
r <- bind_rows(lapply(split(r, paste(r$target, r$fragID)), function(x){
       x$UMI = randomDNAseq()
       x
     }))

# Extract genomic sequences and build read sequences.
r <- bind_rows(parLapply(cluster, split(r, ntile(1:nrow(r), nCPUs)), function(k){
#r <- bind_rows(lapply(split(r, ntile(1:nrow(r), nCPUs)), function(k){
       library(dplyr)
       library(Biostrings)
  
       bind_rows(lapply(split(k, 1:nrow(k)), function(x){
         
         if((x$readStartPos + readWidth) > width(g[names(g) == x$chromosome]) |
            x$readEndPos > width(g[names(g) == x$chromosome])) return(tibble())
         
         a <- subseq(g[names(g) == x$chromosome], x$readStartPos, x$readStartPos + readWidth)
         b <- subseq(g[names(g) == x$chromosome], x$readEndPos - readWidth, x$readEndPos)
  
         if(x$strand == '+'){
            x$R1 <- paste0(x$UMI, 'CTCCGCTTAAGGGACT', as.character(reverseComplement(b)))
            x$R2 <- paste0(ltrBit, as.character(a))
         } else {
            x$R1 <- paste0(x$UMI, 'CTCCGCTTAAGGGACT', as.character(a))
            x$R2 <- paste0(ltrBit, as.character(reverseComplement(b)))
         }
       
         x
       }))
    }))

# Split the reads across replicates while not breaking fragments across replicates.
o <- split(r, r$fragID)
r <- bind_rows(mapply(function(x, n){
       x$replicate <- n
       x
     }, o, ntile(1:length(o), 4), SIMPLIFY = FALSE))

r$sample <- 'sample1'

r <- bind_rows(lapply(split(r, r$sample), function(x){
       x$uniqueLinker = randomDNAseq(n = 20)
  
       bind_rows(lapply(split(x, x$replicate), function(x2){
         x2$barcode <-  randomDNAseq(n = 12)
         x2
       }))
     }))

r$R1 <-paste0(r$uniqueLinker, r$R1)

r$trial <- 'validation'
r$subject <- 'validationSubject'
r$adriftReadLinkerSeq <- paste0(r$uniqueLinker, 'NNNNNNNNNNNNCTCCGCTTAAGGGACT')
r$refGenome <- 'hg38'
r$vectorFastaFile <- 'HXB2.fasta'
r$leaderSeqHMM <- hmmName
r$flags <- sampleFlag

R1 <- ShortRead::ShortReadQ(sread = DNAStringSet(r$R1), 
                            id = BStringSet(r$readID), 
                            quality = BStringSet(sapply(nchar(r$R1), function(x) paste0(rep('?', x), collapse = ''))))

R2 <- ShortRead::ShortReadQ(sread = DNAStringSet(r$R2), 
                            id = BStringSet(r$readID), 
                            quality = BStringSet(sapply(nchar(r$R2), function(x) paste0(rep('?', x), collapse = ''))))

I1 <- ShortRead::ShortReadQ(sread = DNAStringSet(r$barcode), 
                            id = BStringSet(r$readID), 
                            quality = BStringSet(sapply(nchar(r$barcode), function(x) paste0(rep('?', x), collapse = ''))))

ShortRead::writeFastq(R1, file.path(outputDir, paste0(filePrefix, 'syn_R1.fastq.gz')), compress = TRUE)
ShortRead::writeFastq(R2, file.path(outputDir, paste0(filePrefix, 'syn_R2.fastq.gz')), compress = TRUE)
ShortRead::writeFastq(I1, file.path(outputDir, paste0(filePrefix, 'syn_I1.fastq.gz')), compress = TRUE)

sampleData <- distinct(select(r, trial, subject, sample, replicate, adriftReadLinkerSeq, barcode, refGenome, vectorFastaFile, leaderSeqHMM, flags))
names(sampleData) <- c('trial',	'subject', 'sample',	'replicate', 'adriftReadLinkerSeq', 'index1Seq', 'refGenome', 'vectorFastaFile', 'leaderSeqHMM', 'flags')

readr::write_tsv(sampleData, file.path(outputDir, paste0(filePrefix, 'sampleData.tsv')))

stopCluster(cluster)