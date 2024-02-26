# AAVengeR/anchorReadRearrangements.R
# John K. Everett, Ph.D.
# 
# This scripts accepts the a demultiplexed read table from the demultiplex module
# and identifies instances of vector rearrangements. This analysis focus on anchor
# reads and it is often applied to sequencing experiments where anchor reads lengths
# are prioritized over adrift read lengths.

suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))

# Read the configuration file and set additional parameters.
configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')

opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

createOuputDir()
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir))
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs'))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs'))
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'))
invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), full.names = TRUE)))

setMissingOptions()
setOptimalParameters()
set.seed(1)


# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)
write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)


quitOnErorr <- function(msg){
  updateLog(msg)
  message(msg)
  message(paste0('See log for more details: ', opt$defaultLogFile))
  q(save = 'no', status = 1, runLast = FALSE) 
}

updateLog('Reading sample data.')

cluster <- makeCluster(opt$anchorReadRearrangements_CPUs)
clusterSetRNGStream(cluster, 1)

if(! file.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_inputFile))) quitOnErorr('Error - the input data file does not exist.')

# Read in the reads table.
updateLog('Reading in demultiplexed reads.')
reads <- readRDS(file.path(opt$outputDir, opt$anchorReadRearrangements_inputFile))

reads <- reads[nchar(reads$anchorReadSeq) >= opt$anchorReadRearrangements_minAnchorReadLength,]

# Hot fix for SampleSheet errors during demultiplexing.
#
# o <- readr::read_csv('correctedMetaData.csv', col_names = TRUE)
# o$uniqueSample <- paste0(o$trial, '~', o$subject, '~', o$sample, '~', o$replicate)
# reads <- rbindlist(lapply(split(reads, reads$uniqueSample), function(x){
#            x$vectorFastaFile <- subset(o, uniqueSample == x$uniqueSample[1])$vectorFastaFile 
#            x
#          }))
# reads$vectorFastaFile <- sub('BushmanAAVcontrols.fasta', 'BushmanAAVcontrolsPlasmidLargestRemnant.fasta', reads$vectorFastaFile)

if(nrow(reads) == 0) quitOnErorr('Error - the input table contained no rows.')

reads <- select(reads, uniqueSample, readID, anchorReadSeq, adriftReadRandomID, adriftReadTrimSeq, vectorFastaFile)

reads$i <- paste(reads$uniqueSample, reads$anchorReadSeq, reads$adriftReadRandomID, reads$vectorFastaFile)
reads <- reads[! duplicated(reads$i)]
reads$i <- NULL

if(! opt$anchorReadRearrangements_vectorLeaderSeqFilter == 'none'){
  updateLog('Limiting reads to select anchor read start sequences.')
  
  reads$sample <- sub('~\\d+$', '', reads$uniqueSample)
  
  leaderSeqReport <- bind_rows(lapply(split(reads, paste0(reads$vectorFastaFile,reads$sample)), function(x){
    tab <- data.frame(sort(table(substr(x$anchorReadSeq, 1, 15)), decreasing = TRUE)[1:3])
    tab$percent <- sprintf("%.2f%%", (tab$Freq/nrow(x))*100)
    tab$vector <- x$vectorFastaFile[1]
    tab$sample <- x$sample[1]
    names(tab) <- c('first15nts', 'nReads', 'percent', 'vector', 'sample')
    select(tab, vector, sample, nReads, first15nts, percent)
  }))
  
  
  nReadsPreFilter <- nrow(reads)
  reads <- rbindlist(lapply(unlist(strsplit(opt$anchorReadRearrangements_vectorLeaderSeqFilter, '\\|')), function(x){
           o <- unlist(strsplit(x, ';'))
           updateLog(paste0('Limit reads for vector:  ', o[1]))
           
           x <- subset(reads, reads$vectorFastaFile == o[1])
           if(nrow(x) == 0) return(data.table())
           
           rbindlist(lapply(unlist(strsplit(o[2], ',')), function(k){
                updateLog(paste0('  Limiting reads for this vector to those starting with ', k))
                r <- x[vcountPattern(k, 
                                     DNAStringSet(substr(x$anchorReadSeq, 1, nchar(k))), 
                                     max.mismatch = opt$anchorReadRearrangements_vectorLeaderSeqFilterMaxMismatch) == 1]
                updateLog(paste0('    ', ppNum(nrow(r)), ' reads matched the leader sequence.'))
                r
           }))
        }))
  
  if(nrow(reads) == 0) quitOnErorr('Error - not reads matched any of the provided leader sequences.')
  
  updateLog(paste0(sprintf("%.2f%%", (nrow(reads) / nReadsPreFilter)*100), ' reads retained.'))
}


# Trim anchor read over-reading with cutadapt using the RC of the common linker in adrift reads.
# The trim sequence is defined in the reads table created by demultiplex.R.

updateLog('Trimming anchor read over-reading.')


if(tolower(opt$anchorReadRearrangements_additionalAnchorReadOverReadingSeqs) != 'none'){
  opt$anchorReadRearrangements_additionalAnchorReadOverReadingSeqs <- paste0(' -a ', paste0(unlist(strsplit(opt$anchorReadRearrangements_additionalAnchorReadOverReadingSeqs, ',')), collapse = ' -a '))
} else{
  opt$anchorReadRearrangements_additionalAnchorReadOverReadingSeqs = ''
}

clusterExport(cluster, 'opt')

reads <- data.table::rbindlist(parLapply(cluster, split(reads, dplyr::ntile(1:nrow(reads), opt$anchorReadRearrangements_CPUs)), function(x){
  source(file.path(opt$softwareDir, 'lib.R'))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(Biostrings))
  
  # Split reads in this CPU chunk by their adapter sequences.
  data.table::rbindlist(lapply(split(x, x$adriftReadTrimSeq), function(y){
    
    # Write out reads in this chunk to a tmp file.
    f <- file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp',  tmpFile())
    o <- DNAStringSet(y$anchorReadSeq)
    names(o) <- y$readID
    Biostrings::writeXStringSet(o, f)
    
    # Use cut adapt to trim anchor read over-reading.
    system(paste0('cutadapt -e ', opt$anchorReadRearrangements_cutAdaptErrorRate, ' -a ', y$adriftReadTrimSeq[1], ' ', opt$anchorReadRearrangements_additionalAnchorReadOverReadingSeqs, ' ', f, ' > ', paste0(f, '.cutAdapt')), ignore.stderr = TRUE)
    
    # Read in the cutadapt trimmed sequences.
    t <- readDNAStringSet(paste0(f, '.cutAdapt'))
    invisible(file.remove(f, paste0(f, '.cutAdapt')))
    
    # Remove reads that were trimmed too short.
    t <- t[width(t) >= opt$anchorReadRearrangements_minAnchorReadLength]
    if(length(t) == 0) return(data.table())
    
    # Limit incoming reads to those found in the cutAdapt output.
    o <- subset(y, readID %in% names(t))
    trimmed <- data.table(readID = names(t), anchorReadSeq2 = as.character(t))
    
    # Attach the trimmed sequences to the input sequence labeled as anchorReadSeq2.
    data.table(left_join(o, trimmed, by = 'readID'))
  }))
})) %>% select(uniqueSample, readID, adriftReadRandomID, anchorReadSeq, anchorReadSeq2, vectorFastaFile)


# Report trimming stats.
trimReport <- mutate(reads, sample = sub('~\\d+$', '', uniqueSample)) %>% 
  group_by(sample) %>% 
  summarise(reads = n_distinct(readID), percentTrimmed = sprintf("%.2f%%", (sum(anchorReadSeq != anchorReadSeq2)/n())*100)) %>% 
  ungroup() %>%
  readr::write_tsv(file = file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'commonLinkerTrimReport.tsv'), append = FALSE, col_names = TRUE)


# Switch original sequences for trimmed sequences.
reads$anchorReadSeq <- NULL
reads <- dplyr::rename(reads, anchorReadSeq = anchorReadSeq2)


# Identify UMIs that were sequenced several times and use these as a guide
# to correct other UMIs that are different by an edit distance of 1.

tab <- table(reads$adriftReadRandomID)
k <- sort(tab[tab >= opt$anchorReadRearrangements_abundSeqMinCount], decreasing = TRUE)

# We need to limit the number of abundant UMIs otherwise the string distance matrix may grow too large.
if(length(k) > 100) k <- k[1:100]

if(length(k) > 0){
  # stringdistmatrix() will use all CPUs by default.
  m <- stringdist::stringdistmatrix(reads$adriftReadRandomID, names(k))

  singleEditDist <- function(x){
    i <- which(x == 1)
  
    if(length(i) == 1){
      names(k)[i]
    } else {
      return(NA)
    }
  }

  reads$altUMI <- apply(m, 1, singleEditDist)
  reads$adriftReadRandomID <- ifelse(is.na(reads$altUMI), reads$adriftReadRandomID, reads$altUMI)
  reads$altUMI <- NULL
}

g <- group_by(reads, adriftReadRandomID, uniqueSample) %>% summarise(n = n_distinct(readID), .groups = 'drop') 

k <- tibble(reshape2::dcast(g, adriftReadRandomID ~ uniqueSample, value.var = 'n'))
rm(g)

k2 <- k[, 2:length(k)]
k2[is.na(k2)] <- 0
k2[k2 > 0] <- 1
k2$sum <- apply(k2, 1, sum)

o <- bind_cols(k[,1], k2)

rm(k, k1)

# k1 contains UMIs that do not span multiple samples.
o1 <- o[o$sum == 1,]
reads1 <-subset(reads, reads$adriftReadRandomID %in% o1$adriftReadRandomID)

reads2 <- tibble()

# o2 contains UMIs that span multiple samples and requires more scrutiny.
o2 <- o[o$sum > 1,]

if(nrow(o2) > 0){
   o <- dplyr::filter(k, adriftReadRandomID %in% o2$adriftReadRandomID) 

   o$n <- ntile(1:nrow(o2), opt$anchorReadRearrangements_CPUs)

   o2 <- bind_rows(parLapply(cluster, split(o, o$n), function(a){
           library(dplyr)
           a$n <- NULL

           bind_rows(lapply(split(a, 1:nrow(a)), function(d){ 
             x <- d[2:length(d)]
             a <- x[x > 0]
             a <- a[rev(order(a))]
  
             # A sample needs to have 3x the number of reads compared to the next highest sample to be the winner.
             if(a[1] >= a[2]*3){
                  x$keep <- names(x)[which(x == a[1])]
                } else {
                  x$keep <- NA
             }
         
             bind_cols(d[,1], x)
           }))
         }))

  if(any(! is.na(o2$keep))){
    o2 <- o2[! is.na(o2$keep),]
           
    reads$f <- paste(reads$uniqueSample, reads$adriftReadRandomID)
    toKeep  <- paste(o2$keep, o2$adriftReadRandomID)
    reads2  <- subset(reads, f %in% toKeep)
    reads2$f <- NULL
    reads$f <- NULL
  } else {
     reads2 <- tibble()
  }
}

reads <- data.table(bind_rows(reads1, reads2))
rm(o1, o2, reads1, reads2)

d <- reads$adriftReadRandomID[duplicated(reads$adriftReadRandomID)]
a <- reads[reads$adriftReadRandomID %in% d]
b <- reads[! reads$adriftReadRandomID %in% d]
rm(d, reads)


# Create splitting vectors for parLapply that will not break appart UMI groupings.
a$i <- group_by(a, adriftReadRandomID) %>% group_indices
v <-   tibble(i = 1:max(a$i), n = ntile(i, opt$anchorReadRearrangements_CPUs))
a <- left_join(a, v, by = 'i')


# Determine consensus UMI sequences where possible.
a1 <- rbindlist(parLapply(cluster, split(a, a$n), function(k){
        suppressPackageStartupMessages(library(data.table))
        suppressPackageStartupMessages(library(Biostrings))
  
         rbindlist(lapply(split(k, k$adriftReadRandomID), function(x){
           # Limit reads to those near the median read length.
           x <- x[abs(nchar(x$anchorReadSeq) - median(nchar(x$anchorReadSeq))) <= 3,]
           
           # Only return a sequence if two or more reads persist. 
           # Data is in this loop because two or more reads had the same UMI.
           if(nrow(x) >= 2){
             return(tibble(uniqueSample = x$uniqueSample[1], adriftReadRandomID = x$adriftReadRandomID[1], vector = x$vectorFastaFile[1], seq = as.character(consensusString(DNAStringSet(x$anchorReadSeq)))))
           } else {
             return(data.table())
            }
         }))
       }))

# Ensure that consensus sequences are reasonable and did not create several ambiguous codes.
i <- (stringr::str_count(toupper(a1$seq), '[ATCG]')/nchar(a1$seq)) >= 0.95

a1 <- a1[i]

b1 <- tibble(uniqueSample = b$uniqueSample, adriftReadRandomID = b$adriftReadRandomID, vector = b$vectorFastaFile, seq = b$anchorReadSeq)

o <- bind_rows(a1, b1)
rm(a, b, a1, b1)

# Code from mapSiteLeaderSequences.R
#------------------------------------------------------------------------------

clusterExport(cluster, c('opt', 'blast2rearangements', 'buildRearrangementModel'))

m <- rbindlist(lapply(split(o, o$vector), function(x){
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), full.names = TRUE)))
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs'), full.names = TRUE)))
  
  updateLog(paste0('Processing maps for vector: ', x$vector[1], '.'))

  # Clean up tmp database directory.
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs'), full.names = TRUE)))
  
  # Make a blast database for the vector sequence.
  system(paste0('makeblastdb -in ', file.path(opt$softwareDir, 'data', 'vectors', x$vector[1]), 
                ' -dbtype nucl -out ', file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
  
  waitForFile(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs', 'd.nin'))
  
  # Create a data frame of repLeader sequences to map.
  s <- data.frame(seq = unique(x$seq))
  s$id <- paste0('s', 1:nrow(s))
  
  # Bin the repLeader sequences by length using bins incremented by 10 nt.
  s$cut <- cut(nchar(s$seq), c(-Inf, seq(0, max(nchar(s$seq)), by = 10), Inf), labels = FALSE)
  
  # Create a splitting vector that breaks across bins. 
  s <- group_by(s, cut) %>% 
       mutate(n = ntile(1:n(), opt$anchorReadRearrangements_CPUs)) %>% 
       ungroup()
  
  reads <- DNAStringSet(s$seq)
  names(reads) <- s$id
  
  # Chunk the repLeaderSeqs and process each chunk.
  # Here qName stands for a particular repLeaderSeq while in other parts 
  # of the software it typically denotes a single read.

  b <- bind_rows(parLapply(cluster, split(reads, s$n), function(a){
         library(dplyr)
         library(data.table)
         library(Biostrings)
    
         f <- tmpFile()
    
         # Write the chunk out to a tmp file.
         writeXStringSet(a,  file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.fasta')))
    
         # Align the chunk to the vector sequence.
         system(paste0('blastn -dust no -soft_masking false -word_size 5 -evalue 50 -outfmt 6 -query ',
                       file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.fasta')), ' -db ',
                       file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs', 'd'),
                      ' -num_threads 1 -out ', file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.blast'))),
                      ignore.stdout = TRUE, ignore.stderr = TRUE)
    
         waitForFile(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.blast')))
    
         # Catch instances where no alignments were returned.
         if(file.info(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.blast')))$size == 0) return(tibble())
    
         # Parse blastn result.
         b <- read.table(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.blast')), sep = '\t', header = FALSE)
         names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
         b$alignmentLength <- b$qend - b$qstart + 1
         b$strand <- ifelse(b$sstart > b$send, '-', '+')
    
         # Limit repLeaderSeqs in this chunk to those with one or more alignments.
         o <- a[names(a) %in% b$qname]
    
         # Add query lengths to the blastn alignments.
         d <- tibble(qname = names(o), qlen = width(o))
    
         b <- left_join(b, d, by = 'qname')
    
         # Limit alignments based on parameters in the configuration file.
         dplyr::filter(b, pident >= opt$anchorReadRearrangements_seqsMinPercentID, alignmentLength >= opt$anchorReadRearrangements_minAlignmentLength)
       }))
  
  # Create a grouping index for read IDs.
  b$i <- group_by(b, qname) %>% group_indices()
  
  # Create a CPU splitting vector that would not separate the read ID indices. 
  o <- tibble(i2 = 1:n_distinct(b$i))
  o$n <- ntile(1:nrow(o), opt$prepReads_CPUs)
  
  # Bind the CPU splitting vector.
  b <- left_join(b, o, by = c('i' = 'i2'))
  
  # Pass alignment chunks to rearrangement function.
  r <- rbindlist(parallel::parLapply(cluster, split(b, b$n), blast2rearangements, opt$anchorReadRearrangements_maxMissingTailNTs))
  
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), full.names = TRUE)))
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs'), full.names = TRUE)))
  
  z <- left_join(tibble(qname = names(reads), leaderSeq = as.character(reads)), r, by = 'qname')
  x <- left_join(x, select(z, leaderSeq, rearrangement), by = c('seq' = 'leaderSeq'))
  dplyr::select(x, -seq)
}))

m$missingTailAssignment <- grepl('\\[x\\]$', m$rearrangement)
m$numRearrangements <- stringr::str_count(m$rearrangement, ';')

if(opt$anchorReadRearrangements_excludeMissingTailsFromCounts) m[m$missingTailAssignment == TRUE,]$numRearrangements <- m[m$missingTailAssignment == TRUE,]$numRearrangements - 1

m$missingTailAssignment <- NULL

o$sample <- sub('~\\d+$', '', o$uniqueSample)

sampleTotalNTs <- bind_rows(lapply(split(o, o$sample), function(x){
                    tibble(sample = x$sample[1], totalNTs = sum(nchar(x$seq)))
                  }))

# Remove UMIs for which short hands could not be created.
# This should be a small number if the min. alignment length is less than
# the start sequence filters.

m <- m[! is.na(m$rearrangement),]

# Create an alternative rearrangement shorthand that standardized the ends 
# so that shorthands from different length reads can be compared.
rearrangements2 <- bind_rows(lapply(unique(m$rearrangement), function(x){
                      p <- unlist(strsplit(x, ';'))
                      p[length(p)] <- sub('\\.\\.\\d+', '..end', p[length(p)])
                      p[length(p)] <- sub('\\d+\\]$', 'end]', p[length(p)])
                      tibble(rearrangement = x, rearrangement2 = paste0(p, collapse = ';'))
                     }))


m$sample <- sub('~\\d+$', '', m$uniqueSample)

m <- left_join(m, rearrangements2, by = 'rearrangement')

r <- group_by(m, sample) %>%
     summarise(nUMIs = n_distinct(adriftReadRandomID), 
               recombinationEvents = sum(numRearrangements), 
               percentRecombinedUMIs = sprintf("%.2f%%", (sum(numRearrangements > 0)/nUMIs)*100)) %>%
     ungroup() 

r <- left_join(r, sampleTotalNTs, by = 'sample')
r$recombinationEventsPer1K  <- r$recombinationEvents / (r$totalNTs/1000)
r$recombinationEventsPer5K  <- r$recombinationEvents / (r$totalNTs/5000)

readr::write_tsv(r, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'summary.tsv'))
openxlsx::write.xlsx(r, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'summary.xlsx'))

saveRDS(m, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'output.rds'))
