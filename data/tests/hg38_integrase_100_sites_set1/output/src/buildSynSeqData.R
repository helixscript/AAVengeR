#!/usr/bin/Rscript

library(dplyr)
library(Biostrings)
library(parallel)
library(yaml)

args <- commandArgs(trailingOnly=TRUE)

if(length(args) == 0) stop('Error - configuration file not provided as first argument.')
if(! file.exists(args)) stop('Error - the cofiguration file could not be found.')

opt <- read_yaml(args)

if(opt$mode == 'AAV') opt$addIntegraseCorrection <- FALSE

set.seed(opt$seed)

remnant0   <- 'TCTGCGCGCTCGCTCGCTCA' # To be used for integrase mode.
remnant    <- 'TCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTG' 
linker     <- 'GAACGAGCACTAGTAAGCCCNNNNNNNNNNNNCTCCGCTTAAGGGACT' 

# Create output directories.
if(! dir.exists(opt$outputDir)) dir.create(opt$outputDir)

# Create a collection of remnant pieces to draw from to simulate rearrangements.
remnantChunks <- unique(unlist(sapply(1:1000, function(x){
  p1 <- sample(1:(nchar(remnant) - 15), 1)
  p2 <-  sample(c((p1+15):(p1+45)), 1)
  o <- substr(remnant, p1, p2)
  o[nchar(o) >= 15 & nchar(o) <= 30]
})))

# Read in genome, already filtered for main chromosomes.
g <- rtracklayer::import.2bit(opt$dbPath)
g <- g[! names(g) %in% c('chrM', 'chrY')]

# Create fragments as dual detection pairs with some wiggle added to their start positions.
n <- 1
f <- bind_rows(lapply(sample(names(g), opt$nSitesPool, replace = TRUE), function(x){
       message(paste0(n, ' / ', opt$nSitesPool, ': ', sprintf("%.2f%%", (n/opt$nSitesPool)*100)))
  
       # Retrieve 10 1000 NT wide fragments.
       pos <- sample(10000:(width(g[names(g) == x])-10000), 10)
       frags <- Reduce('append', lapply(pos, function(p) subseq(g[names(g) == x], start = p, width = 10000)))
       names(frags) <- paste('pos', pos)  
  
       # Remove fragments with any Ns and select the first one.
       frags <- frags[! sapply(as.character(frags), function(x) stringr::str_detect(x, 'N'))]
       if(length(frags) == 0) stop('Error -- could not find a fragment without Ns.')
       frag  <- frags[1]
       
       pos1 <- as.integer(stringr::str_extract(names(frag), '\\d+')) + 500 - 1

       posFragStartPositions <- round(rnorm(opt$nFragsperSite * opt$nReadsPerFrag, mean = 500, sd = 0.5))
   
       posFrags <- tibble(chr = x,
                          id = paste0(x, '+', pos1, ':sitePair', n, ':', unlist(lapply(1:opt$nFragsperSite, function(x) rep(paste0('frag', x), opt$nReadsPerFrag))), paste0(':read', 1:opt$nReadsPerFrag)),
                          position = pos1, 
                          strand = '+', 
                          seq = sapply(posFragStartPositions, function(s) substr(as.character(frag), s, 1000)))
       pos2 <- pos1 - 3
       
       negFragEndPositions <- round(rnorm(opt$nFragsperSite * opt$nReadsPerFrag, mean = 497, sd = 0.5))
       negFrags <- tibble(chr = x,
                          id = paste0(x, '-', pos2, ':sitePair', n, ':',  unlist(lapply(1:opt$nFragsperSite, function(x) rep(paste0('frag', x), opt$nReadsPerFrag))), paste0(':read', 1:opt$nReadsPerFrag)),
                          position = pos2, 
                          strand = '-', 
                          seq = sapply(negFragEndPositions, function(s) as.character(reverseComplement(DNAStringSet(substr(as.character(frag), 1, s))))))

       if(opt$buildSitePairs){
         o <- bind_rows(posFrags, negFrags)
       } else {
         if(sample(c(TRUE, FALSE), 1)){
           o <- posFrags
         } else {
           o <- negFrags
         }
       }
       
       # Update positions for the actions of buildSites when integrase mode is used.
       if(opt$addIntegraseCorrection){
         o$position <- ifelse(o$strand == '+', o$position+2, o$position-2)
         o <- bind_rows(lapply(split(o, 1:nrow(o)), function(x){ x$id <- sub('\\d+:', paste0(x$position, ':'), x$id); x}))
       }
       
       o$posid <- paste0(o$chr, o$strand, o$position)
       
       n <<- n+1
       o
     }))


if(opt$buildSitePairs){
  f$pair <- stringr::str_extract(f$id, 'sitePair\\d+')
  pairs <- sample(unique(f$pair), opt$nSites)
  f <- subset(f, pair %in% pairs)
} else {
  f <- subset(f, posid %in% sample(unique(f$posid), opt$nSites))
}


# Break points are refined within replicates -- posid fragments can not be split across replicates.

# Create fragment breaks and remnant sequences.
f2 <- bind_rows(lapply(split(f, f$posid), function(x){
        if(opt$mode == 'AAV'){
          x$r <- substr(remnant, 1, sample(15:25, 1))
          n <- sample(0:3, 1)
          if(n != 0) x$r <- paste0(x$r[1], sample(remnantChunks, n), collapse = '')
      } else {
        x$r <- remnant0
      }
      
      x$endPositions <- unlist(lapply((1:opt$nFragsperSite * opt$fragStepSize) + opt$fragStartSize, function(x) round(rnorm(opt$nReadsPerFrag, mean = x, sd = 0.5))))
      x$seq <- substr(x$seq, 1, x$endPositions)
      x
     }))

# Add adapters, linkers, and over-reading to read sequences.
f3 <- bind_rows(lapply(split(f2, f2$id), function(x){
        x$R1 <- paste0(linker, as.character(reverseComplement(DNAString(x$seq))))
        x$R2 <- paste0(x$r, x$seq)
        select(x, - seq, -r, -endPositions)
        
        # Simulate 10-15 NTs of over-reading for fragment 1 reads.
        if(grepl('_frag1_', x$id)){
          x$R2 <- paste0(x$R2, as.character(reverseComplement(DNAString(substr(linker, nchar(linker) - sample(10:15, 1), nchar(linker))))))
          x$R1 <- paste0(x$R1, as.character(reverseComplement(DNAString(substr(x$r, nchar(x$r) - sample(10:15, 1), nchar(x$r))))))
        }
        x
      }))

randomDNAseq <- function(x, n = 12){
  paste(paste0(stringi::stri_rand_strings(n, 1, '[ATCG]')), collapse = '')
}

f3 <- tidyr::separate(f3, id, c('posid', 'sitePair', 'fragNum', 'readNum'), sep = ':', remove = FALSE)

f3 <- bind_rows(lapply(split(f3, paste(f3$sitePair, f3$fragNum)), function(x){
        seq <- randomDNAseq()
        x$R1 <- gsub('[N]+', seq, x$R1)
        x
      }))

o <- list()
o[['sample1']] <- c('GGCTAAACTATG', 'TCAACCCGTGAA')
o[['sample2']] <- c('ATCAGAGCCCAT', 'GTGTGCTAACGT')
o[['sample3']] <- c('CTTGCGGCAATC', 'TACCTAGTGAGA')

f3 <- bind_rows(lapply(split(f3, f3$posid), function(x){
        # Assign random sample to this site.
        x$sample <- sample(names(o), 1)
        
        # Assign each fragment to one of two isolates.
        x$replicate <- (as.integer(stringr::str_extract(x$fragNum, '\\d+')) %% 2) + 1
        
        # Assign sample / replicate specific bar code.
        x$barcode <- o[[x$sample[1]]][x$replicate]
        
        x
      }))

f3$trial <- 'validation'
f3$subject <- 'validationSubject'
f3$refGenome <- 'hg38'
f3$leaderSeqHMM <- 'validation.hmm'
f3$vectorFastaFile <- 'validationVector.fasta'
f3$adriftReadLinkerSeq <- 'GAACGAGCACTAGTAAGCCCNNNNNNNNNNNNCTCCGCTTAAGGGACT'

R1 <- ShortRead::ShortReadQ(sread = DNAStringSet(f3$R1), 
                            id = BStringSet(f3$id), 
                            quality = BStringSet(sapply(nchar(f3$R1), function(x) paste0(rep('?', x), collapse = ''))))

R2 <- ShortRead::ShortReadQ(sread = DNAStringSet(f3$R2), 
                            id = BStringSet(f3$id), 
                            quality = BStringSet(sapply(nchar(f3$R2), function(x) paste0(rep('?', x), collapse = ''))))

I1 <- ShortRead::ShortReadQ(sread = DNAStringSet(f3$barcode), 
                            id = BStringSet(f3$id), 
                            quality = BStringSet(sapply(nchar(f3$barcode), function(x) paste0(rep('?', x), collapse = ''))))

if(file.exists(file.path(opt$outputDir, 'syn_I1.fastq.gz'))) invisible(file.remove(file.path(opt$outputDir, 'syn_I1.fastq.gz')))
if(file.exists(file.path(opt$outputDir, 'syn_R1.fastq.gz'))) invisible(file.remove(file.path(opt$outputDir, 'syn_R1.fastq.gz')))
if(file.exists(file.path(opt$outputDir, 'syn_R2.fastq.gz'))) invisible(file.remove(file.path(opt$outputDir, 'syn_R2.fastq.gz')))

ShortRead::writeFastq(R1, file.path(opt$outputDir, 'syn_R1.fastq.gz'), compress = TRUE, mode = 'w')
ShortRead::writeFastq(R2, file.path(opt$outputDir, 'syn_R2.fastq.gz'), compress = TRUE, mode = 'w')
ShortRead::writeFastq(I1, file.path(opt$outputDir, 'syn_I1.fastq.gz'), compress = TRUE, mode = 'w')

sampleData <- distinct(select(f3, trial, subject, sample, replicate, adriftReadLinkerSeq, barcode, refGenome, vectorFastaFile, leaderSeqHMM))

if(opt$mode == 'AAV'){
  sampleData <- select(sampleData, -leaderSeqHMM)
  sampleData$anchorReadStartSeq <- substr(remnant, 1, 10)
  names(sampleData) <- c('trial',	'subject', 'sample',	'replicate', 'adriftReadLinkerSeq', 'index1Seq', 'refGenome', 'vectorFastaFile', 'anchorReadStartSeq')
  sampleData$flags <- 'AAV'
} else {
  names(sampleData) <- c('trial',	'subject', 'sample',	'replicate', 'adriftReadLinkerSeq', 'index1Seq', 'refGenome', 'vectorFastaFile', 'leaderSeqHMM')
  sampleData$flags <- 'IN_u5'
}

readr::write_tsv(sampleData, file.path(opt$outputDir,'sampleData.tsv'), append = FALSE)

group_by(f3, trial, subject, sample, posid) %>% 
summarise(nReads = n_distinct(id), 
          nFrags = n_distinct(fragNum), 
          nUMIs = n_distinct(substr(R1, 21, 32)),
          leaderSeq = r[1]) %>% 
ungroup() %>%
readr::write_tsv(file.path(opt$outputDir,'truth.tsv'), append = FALSE)

file.copy('config.yml', file.path(opt$outputDir,'config.yml'))

p <- file.path(opt$outputDir, 'syn_R2.fastq.gz')
p <- gsub("/", "\\\\/", p)
system(paste0("sed -i 's/Undetermined_S0_R2_001.fastq.gz/", p, "/' ", file.path(opt$outputDir,'config.yml')))

p <- file.path(opt$outputDir, 'syn_R1.fastq.gz')
p <- gsub("/", "\\\\/", p)
system(paste0("sed -i 's/Undetermined_S0_R1_001.fastq.gz/", p, "/' ", file.path(opt$outputDir,'config.yml')))

p <- file.path(opt$outputDir, 'syn_I1.fastq.gz')
p <- gsub("/", "\\\\/", p)
system(paste0("sed -i 's/Undetermined_S0_I1_001.fastq.gz/", p, "/' ", file.path(opt$outputDir,'config.yml')))

p <- file.path(opt$outputDir, 'sampleData.tsv')
p <- gsub("/", "\\\\/", p)
system(paste0("sed -i 's/sampleData.tsv/", p, "/' ", file.path(opt$outputDir,'config.yml')))

p <- paste0(getwd(), '/', opt$outputDir, '/output')
p <- gsub("/", "\\\\/", p)
system(paste0("sed -i -E '0,/outputDir:\\s\\S+/s//outputDir: ",  p, "/' ", file.path(opt$outputDir,'config.yml')))

p <- getwd()
p <- gsub("/", "\\\\/", p)
system(paste0("sed -i -E 's/softwareDir:\\s\\S+/softwareDir: ",  p, "/' ", file.path(opt$outputDir,'config.yml')))

p <- ifelse(grepl('integrase', args), 'integrase', 'AAV')
system(paste0("sed -i -E 's/mode:\\s\\S+/mode: ",  p, "/' ", file.path(opt$outputDir,'config.yml')))

p <- '~\\/.my.cnf'
system(paste0("sed -i -E 's/databaseConfigFile:\\s\\S+/databaseConfigFile: ",  p, "/' ", file.path(opt$outputDir,'config.yml')))
 
parallel::detectCores()
system(paste0("sed -i -E 's/core_CPUs:\\s\\S+/core_CPUs: ", parallel::detectCores(), "/' ", file.path(opt$outputDir,'config.yml')))

system(paste0('./aavenger.R ', file.path(opt$outputDir,'config.yml')))


q()
