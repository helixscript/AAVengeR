library(dplyr)
library(Biostrings)
library(parallel)

set.seed(1)

dbPath    <- '/home/ubuntu/AAVengeR/data/referenceGenomes/blat/GCA_009914755.4.2bit'
outputDir <- '/home/ubuntu/AAVengeR/data/validation'

nCPUs  <- 15
nSites <- 1000
readWidth <- 150
fragWidthStart <- 500
nFragsPerSite  <- 5
nReadsPerFrag  <- 25


ltrBit <- 'GAAAATCTCTAGCA'
                               
# Read in genome, exclude mitochondria and Y chromosomes.
g <- rtracklayer::import.2bit(dbPath)
g <- g[! names(g) %in% c('chrM', 'chrY')]

cluster <- makeCluster(nCPUs)
clusterExport(cluster, c('ltrBit', 'g', 'readWidth'))

# Select positions positions in the genome that are not near any Ns.
targets <- bind_rows(lapply(sample(names(g), nSites, replace = TRUE), function(x){
             # Retrieve 10 1000 NT wide fragments.
             pos <- sample(10000:(width(g[names(g) == x])-10000), 10)
             frags <- Reduce('append', lapply(pos, function(p) subseq(g[names(g) == x], start = p, width = 10000)))
             names(frags) <- paste('pos', pos)  
             
             # Remove fragments with any Ns and select the first one.
             frags <- frags[! sapply(as.character(frags), function(x) stringr::str_detect(x, 'N'))]
             frag  <- frags[1]
             
             # Define the position as the center of the 10,000 NT frag without Ns.
             tibble(chromosome = x, 
                    strand = sample(c('+', '-'), 1), 
                    position = as.integer(stringr::str_extract(names(frag), '\\d+')) + 5000)
}))
                 

# Create a series of values around a mean beneath a tight Gaussian curve. 
createPositions <- function(pos, n, d = 3){
  o <- round(rnorm((n*10), pos, sd = 0.75)) # Create 10x more than needed
  o <- o[abs(o - pos) <= d]                 # Select values within a distance of d from pos
  o[1:n]                                    # return a selection of positions
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
                     fragID = sapply(0:(nFragsPerSite-1), function(n){ paste0('fragment', n) }),
                     endPos = sapply(0:(nFragsPerSite-1), function(n){ x$position + fragWidthStart + n*10}))
       } else {
         d <- tibble(target = paste0(x$chromosome, x$strand, x$position),
                     startPos = sapply(0:(nFragsPerSite-1), function(n){ x$position - (fragWidthStart + n*10)}),
                     fragID = sapply(0:(nFragsPerSite-1), function(n){ paste0('fragment', n) }),
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
       library(dplyr)
       library(Biostrings)
  
       bind_rows(lapply(split(k, 1:nrow(k)), function(x){
  
         a <- subseq(g[names(g) == x$chromosome], x$readStartPos, x$readStartPos + readWidth)
         b <- subseq(g[names(g) == x$chromosome], x$readEndPos - readWidth, x$readEndPos)
  
         if(x$strand == '+'){
            x$R1 <- paste0('xxx', x$UMI, 'CTCCGCTTAAGGGACT', as.character(reverseComplement(b)))
            x$R2 <- paste0(ltrBit, as.character(a))
         } else {
            x$R1 <- paste0('xxx', x$UMI, 'CTCCGCTTAAGGGACT', as.character(a))
            x$R2 <- paste0(ltrBit, as.character(reverseComplement(b)))
         }
       
         x
       }))
    }))

# Distribute targets across five samples.
o <- tibble(target = unique(r$target))
o$nSample <- ntile(o$target, 5) 
r <- left_join(r, o, by = 'target')
r$sample <- paste0('sample', r$nSample)

# Within each sample, randomize reads across four replicates. 
r <- bind_rows(lapply(split(r, r$sample), function(x){
       x <- x[sample(1:nrow(x)),]
       x$replicate <- ntile(1:nrow(x), 4)
       x
     }))

# Build sample data table.
barcodes <- tibble(sample = unique(r$sample),
                   uniqueLinker = sapply(1:n_distinct(r$sample), randomDNAseq, n = 20))

barcodes <- bind_rows(lapply(split(barcodes, barcodes$uniqueLinker), function(x){
  tibble(sample = x$sample, replicate = 1:4, uniqueLinker = x$uniqueLinker, index1Seq = sapply(1:4, randomDNAseq))
}))

# Go back to the reads table and add unique linkers and barcodes. 

r <- bind_rows(lapply(split(r, paste(r$sample, r$replicate)), function(x){
       b <- subset(barcodes, sample == x$sample[1] & replicate == x$replicate[1])
       x$R1 <- sub('xxx', b$uniqueLinker, x$R1)
       x$index1Seq <- b$index1Seq
       x
     }))

barcodes$trial <- 'validation'
barcodes$subject <- 'validationSubject'
barcodes$adriftReadLinkerSeq <- paste0(barcodes$uniqueLinker, 'NNNNNNNNNNNNCTCCGCTTAAGGGACT')
barcodes$refGenome <- 'GCA_009914755.4'
barcodes$vectorFastaFile <- 'HXB2.fasta'
barcodes$leaderSeqHMM <- 'HXB2_U5.hmm'
barcodes$flags <- 'IN_u5'

R1 <- ShortRead::ShortReadQ(sread = DNAStringSet(r$R1), 
                            id = BStringSet(r$readID), 
                            quality = BStringSet(sapply(nchar(r$R1), function(x) paste0(rep('?', x), collapse = ''))))

R2 <- ShortRead::ShortReadQ(sread = DNAStringSet(r$R2), 
                            id = BStringSet(r$readID), 
                            quality = BStringSet(sapply(nchar(r$R2), function(x) paste0(rep('?', x), collapse = ''))))

I1 <- ShortRead::ShortReadQ(sread = DNAStringSet(r$index1Seq), 
                            id = BStringSet(r$readID), 
                            quality = BStringSet(sapply(nchar(r$index1Seq), function(x) paste0(rep('?', x), collapse = ''))))

ShortRead::writeFastq(R1, file.path(outputDir, 'fullTestData_R1.fastq.gz'), compress = TRUE)
ShortRead::writeFastq(R2, file.path(outputDir, 'fullTestData_R2.fastq.gz'), compress = TRUE)
ShortRead::writeFastq(I1, file.path(outputDir, 'fullTestData_I1.fastq.gz'), compress = TRUE)

readr::write_tsv(barcodes, file.path(outputDir, 'sampleData.tsv'))
