#!/usr/bin/Rscript

# John K. Everett, PhD
# AAVengeR/prepReads.R
#
# This module prepares reads for alignment to a reference genome by removing duplicate
# read pairs, trimming adapter sequences, removing reads with high sequence homology to 
# the vector and selecting anchor reads with the expected leader sequences. 

suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(RMariaDB))
options(stringsAsFactors = FALSE)

set.seed(1)

source(file.path(yaml::read_yaml(commandArgs(trailingOnly=TRUE)[1])$softwareDir, 'lib.R'))
opt <- loadConfig()
optionsSanityCheck()

createOuputDir()
if(! dir.exists(file.path(opt$outputDir, opt$prepReads_outputDir))) dir.create(file.path(opt$outputDir, opt$prepReads_outputDir))
if(! dir.exists(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs'))) dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs'))
if(! dir.exists(file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'))) dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'))
invisible(file.remove(list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'), full.names = TRUE)))

# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, opt$prepReads_outputDir, 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)
write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)

quitOnErorr <- function(msg){
  if(opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()
  updateLog(msg)
  message(msg)
  message(paste0('See log for more details: ', opt$defaultLogFile))
  q(save = 'no', status = 1, runLast = FALSE) 
}

updateLog('Reading sample data.')
samples <- loadSamples()

if(tolower(opt$prepReads_additionalAnchorReadOverReadingSeqs) != 'none'){
  opt$prepReads_additionalAnchorReadOverReadingSeqs <- paste0(' -a ', paste0(unlist(strsplit(opt$prepReads_additionalAnchorReadOverReadingSeqs, ',')), collapse = ' -a '))
} else{
  opt$prepReads_additionalAnchorReadOverReadingSeqs = ''
}

cluster <- makeCluster(opt$prepReads_CPUs)
clusterSetRNGStream(cluster, 1)
clusterExport(cluster, 'opt')

if(! file.exists(file.path(opt$outputDir, opt$prepReads_readsTable))) quitOnErorr('Error - the input data file does not exist.')

# Read in the reads table.
updateLog('Reading in demultiplexed reads.')
reads <- readRDS(file.path(opt$outputDir, opt$prepReads_readsTable))

if(nrow(reads) == 0) quitOnErorr('Error - the input table contained no rows.')

incomingSamples <- unique(reads$uniqueSample)

if(previousSampleDatabaseCheck(tibble(uniqueSample = reads$uniqueSample, refGenome = reads$refGenome) %>% tidyr::separate(uniqueSample, c('trial', 'subject', 'sample', 'replicate'), sep = '~') %>% distinct())) q(save = 'no', status = 1, runLast = FALSE) 

# Trim anchor read over-reading with cutadapt using the RC of the common linker in adrift reads.
# The trim sequence is defined in the reads table created by demultiplex.R.

updateLog('Trimming anchor read over-reading.')

reads <- data.table::rbindlist(parLapply(cluster, split(reads, dplyr::ntile(1:nrow(reads), opt$prepReads_CPUs)), function(x){
           source(file.path(opt$softwareDir, 'lib.R'))
           suppressPackageStartupMessages(library(dplyr))
           suppressPackageStartupMessages(library(data.table))
           suppressPackageStartupMessages(library(Biostrings))
  
           # Split reads in this CPU chunk by their adapter sequences.
           data.table::rbindlist(lapply(split(x, x$adriftReadTrimSeq), function(y){
             
             # Write out reads in this chunk to a tmp file.
             f <- file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp',  tmpFile())
             o <- DNAStringSet(y$anchorReadSeq)
             names(o) <- y$readID
             Biostrings::writeXStringSet(o, f)
         
             # Use cut adapt to trim anchor read over-reading.
             system(paste0('cutadapt -e ', opt$prepReads_cutAdaptErrorRate, ' -a ', y$adriftReadTrimSeq[1], ' ', opt$prepReads_additionalAnchorReadOverReadingSeqs, ' ', f, ' > ', paste0(f, '.cutAdapt')), ignore.stderr = TRUE)
         
             # Read in the cutadapt trimmed sequences.
             t <- readDNAStringSet(paste0(f, '.cutAdapt'))
             invisible(file.remove(f, paste0(f, '.cutAdapt')))
            
             # Remove reads that were trimmed too short.
             t <- t[width(t) >= opt$prepReads_minAnchorReadLength]
             if(length(t) == 0) return(data.table())
         
             # Add back trimmed anchor read sequences.
             o <- subset(y, readID %in% names(t))
             trimmed <- data.table(readID = names(t), anchorReadSeq2 = as.character(t))
             o <- left_join(o, trimmed, by = 'readID')
         
             # Trim off adrift read adapters.
             o$adriftReadSeq2 <- substr(o$adriftReadSeq, o$adriftLinkerSeqEnd+1, nchar(o$adriftReadSeq))
             o <- o[nchar(o$adriftReadSeq2) >= opt$prepReads_minAdriftReadLength,]

             data.table(dplyr::select(o, -adriftReadTrimSeq, -adriftLinkerSeqEnd))
          }))
        }))

stopCluster(cluster)

if(nrow(reads) == 0) quitOnErorr('Error - no reads remaining after trimming.')

# Rename the trimmed reads back to their proper column names.
reads$anchorReadSeq <- NULL
reads$adriftReadSeq <- NULL
reads <- dplyr::rename(reads, anchorReadSeq = anchorReadSeq2, adriftReadSeq = adriftReadSeq2)


# Clean up.
invisible(file.remove(list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'), full.names = TRUE)))


mappings <- tibble()
vectorHits <- tibble()

# Align the ends of anchor reads to the vector to identify vector reads which should be removed.
reads$vectorFastaFile <- file.path(opt$softwareDir, 'data', 'vectors', reads$vectorFastaFile)


# Create a CPU cluster.
cluster <- makeCluster(opt$prepReads_CPUs)
clusterSetRNGStream(cluster, 1)
clusterExport(cluster, c('tmpFile', 'waitForFile', 'opt', 'lpe', 'blastReads', 'alignReadEndsToVector'))


# If requested, align the ends of anchor reads to the vector.
# Reads with high homology to the vector are likely internal vector reads or reads from vector 
# concatemers that will not be able to be mapped to a reference genome. This alignment uses blast 
# for sensitivity and because the queries are short and the database (vector) is small.

if(opt$prepReads_excludeAnchorReadVectorHits | opt$prepReads_excludeAdriftReadVectorHits){
  
  updateLog('Aligning ends of anchor reads to vector sequences.')
  
  vectorHits <- bind_rows(lapply(split(reads, reads$vectorFastaFile), function(x){
            invisible(file.remove(list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs'), full.names = TRUE)))
  
            # Mask lowercase letters in vector file because this is the means by which users signal repeat sequences.
            seq <- readLines(x$vectorFastaFile[1])
            seq <- gsub('[atcg]', 'N', seq)
            f <- tmpFile()
            
            write(seq, file = file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp', f))
            
            system(paste0('makeblastdb -in ', file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp', f), 
                          ' -dbtype nucl -out ', file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd')), ignore.stderr = TRUE, ignore.stdout = TRUE)
            
            waitForFile(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd.nin'))
            
            invisible(file.remove(file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp', f)))
            
            testSeqs <- tibble()

            if(opt$prepReads_excludeAnchorReadVectorHits){
              testSeqs <- bind_rows(testSeqs, tibble(readID = paste0(x$readID, '.1'), 
                                                   testSeq = substr(x$anchorReadSeq, (nchar(x$anchorReadSeq) - opt$prepReads_vectorAlignmentTestLength) + 1, nchar(x$anchorReadSeq)) ))
            }
            
            if(opt$prepReads_excludeAdriftReadVectorHits){
              testSeqs <- bind_rows(testSeqs, tibble(readID = paste0(x$readID, '.2'), 
                                                   testSeq = substr(x$adriftReadSeq, (nchar(x$adriftReadSeq) - opt$prepReads_vectorAlignmentTestLength) + 1, nchar(x$adriftReadSeq)) ))
            }
            
            # remove duplicated anchor reads.
            testSeqs$i <- group_by(testSeqs, testSeq) %>% group_indices() 
            uniqueTestSeqs <- testSeqs[! duplicated(testSeqs$i),] 
            
            r <- bind_rows(parLapply(cluster, split(uniqueTestSeqs, ntile(1:nrow(uniqueTestSeqs), opt$prepReads_CPUs)), alignReadEndsToVector))
            
            # Add the anchorRead sequence index to the blastn result and then use it to expand the result out to include all read ids.
            if(nrow(r) > 0){
              r$qname <- sub('\\.\\d$', '', r$qname)
              testSeqs$readID <- sub('\\.\\d$', '', testSeqs$readID)
              
              return(distinct(left_join(r, select(testSeqs, readID, i), by = c('qname' = 'readID'), relationship = 'many-to-many') %>%
                              left_join(select(testSeqs, readID, i), by = 'i', relationship = 'many-to-many') %>% 
                              dplyr::select(-qname, -i) %>% 
                              dplyr::rename(qname = readID) %>%
                              dplyr::relocate(qname, .before = sseqid)))
            } else {
              return(tibble())
            }
          }))
  
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'), full.names = TRUE)))
}

parallel::stopCluster(cluster)
saveRDS(vectorHits, file.path(opt$outputDir, opt$prepReads_outputDir, 'vectorHits.rds'), compress = opt$compressDataFiles)


# If alignments to the vector are found, create reports and remove them from reads table.

nReadsPreFilter <- n_distinct(reads$readID)

if(nrow(vectorHits) > 0){
  updateLog(paste0(ppNum(n_distinct(vectorHits)), ' reads aligned to the vector and will be exlcuded.'))
  reads <- subset(reads, ! readID %in% vectorHits$qname)
} else {
  updateLog('No reads aligned to the vector.')
}  

if(nrow(reads) == 0) quitOnErorr('Error - no reads remain after filtering for reads aligning to the vector.')
  
nReadsPostFilter <- n_distinct(reads$readID)

updateLog(paste0(sprintf("%.2f%%", 100 - (nReadsPostFilter/nReadsPreFilter)*100), ' unique reads removed because they aligned to the vector.'))


if(! 'leaderSeqHMM' %in% names(reads)){
    # Now align the full anchor reads to the vector excluding those in vectorHits.
    updateLog('Aligning full anchor reads to vector sequences.')
  
    cluster <- makeCluster(opt$prepReads_CPUs)
    clusterSetRNGStream(cluster, 1)
    clusterExport(cluster, c('tmpFile', 'waitForFile', 'opt', 'lpe', 'blastReads', 'blast2rearangements_worker'))
  
    vectorHits2 <- rbindlist(lapply(split(reads, reads$vectorFastaFile), function(x){
      invisible(file.remove(list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs'), full.names = TRUE)))
    
      # No masking here here
      
      system(paste0('makeblastdb -in ', x$vectorFastaFile[1], ' -dbtype nucl -out ', file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd')), ignore.stderr = TRUE, ignore.stdout = TRUE)
      
      waitForFile(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd.nin'))
    
      x$i <- group_by(x, anchorReadSeq) %>% group_indices() 
      x2 <- x[! duplicated(x$i),] 
      
      r <- bind_rows(parLapply(cluster,split(x2, ntile(1:nrow(x2), opt$prepReads_CPUs)), function(y){
             suppressPackageStartupMessages(library(Biostrings))
             suppressPackageStartupMessages(library(data.table))
             suppressPackageStartupMessages(library(dplyr))
        
             s <- DNAStringSet(y$anchorReadSeq)
             names(s) <- y$readID
        
             b <- blastReads(s, tmpDirPath = file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'), dbPath = file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd'))
        
             if(nrow(b) > 0){
               b$alignmentLength <- b$qend - b$qstart + 1
               b <- subset(b, pident >= opt$prepReads_mapLeaderSeqsMinPercentID & 
                           alignmentLength >= opt$prepReads_mapLeaderSeqsMinAlignmentLength &
                           gapopen <= 1)
             }
         
             b
           }))
      
      if(nrow(r) > 0){
        return(distinct(left_join(r, select(x, readID, i), by = c('qname' = 'readID')) %>%
                        left_join(select(x, readID, i), by = 'i') %>% 
                        dplyr::select(-qname, -i) %>% 
                         dplyr::rename(qname = readID) %>%
                        dplyr::relocate(qname, .before = sseqid)))
      } else {
        return(tibble())
      }
    }))
    
    if(nrow(vectorHits2) == 0) quitOnErorr('Error - no reads aligned to the vector for leader seq generation.')
    
    updateLog('Creating read maps from local alignments to the vector file.')
    
    # Split reads into CPU chunks without dividing blast results for a read between chunks.
    # Splitting vector is n.
    vectorHits2$i <- group_by(vectorHits2, qname) %>% group_indices() 
    o <- tibble(i2 = 1:n_distinct(vectorHits2$i))
    o$n <- ntile(1:nrow(o), opt$prepReads_CPUs)

    vectorHits2 <- left_join(vectorHits2, o, by = c('i' = 'i2'))
    
    if(! opt$prepReads_buildReadMaps_blastReconstruction){
      
      b2 <- dplyr::filter(vectorHits2, qstart <= opt$prepReads_buildReadMaps_minMapStartPostion) %>%
            dplyr::group_by(qname) %>% 
            dplyr::slice_min(evalue, n = 1) %>% 
            dplyr::slice(1) %>% 
            dplyr::ungroup()
      
      m <- data.table(id = b2$qname, leaderMapping.qStart = 1, leaderMapping.qEnd = b2$qend, leaderSeqMap = NA)
      
    } else {
      o <- rbindlist(parLapply(cluster, split(vectorHits2, vectorHits2$n), blast2rearangements_worker))
      m <- data.table(id = o$qname, leaderMapping.qStart = 1, leaderMapping.qEnd = o$end, leaderSeqMap = NA)
    }
    
    if(opt$prepReads_limitLeaderSeqsWithQuickAlignFilter){
      m <- left_join(m, dplyr::select(reads, readID, quickFilterStartPos) %>% dplyr::filter(readID %in% m$id), by = c('id' = 'readID'))
      m$leaderMapping.qEnd <- ifelse(m$leaderMapping.qEnd > (m$quickFilterStartPos - 1), (m$quickFilterStartPos - 1), m$leaderMapping.qEnd)
    }
    
    parallel::stopCluster(cluster)
    
} else {
  updateLog('Using leader sequence HMM to define mappings.')
  
  cluster <- parallel::makeCluster(opt$prepReads_CPUs)
  clusterSetRNGStream(cluster, 1)
  clusterExport(cluster, c('opt'))
  
  hmmResults <- rbindlist(parLapply(cluster, split(reads, reads$uniqueSample), function(x){
    suppressPackageStartupMessages(library(Biostrings))
    suppressPackageStartupMessages(library(dplyr))
    source(file.path(opt$softwareDir, 'lib.R'))
    
    seqs <- DNAStringSet(x$anchorReadSeq)
    names(seqs) <- x$readID
    captureHMMleaderSeq(seqs, x$leaderSeqHMM[1], tmpDirPath = file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'))
  }))
  
  parallel::stopCluster(cluster)
  
  if(nrow(hmmResults) > 0){
    m <- tibble(id = hmmResults$id, leaderMapping.qStart = 1, leaderMapping.qEnd = nchar(hmmResults$LTRseq), leaderSeqMap = NA)
  } else {
    quitOnErorr('Error - no reads matched the HMMs.')
  }
}

# If anchorReadStartSeq was included in the sample data file then reads must of had a decent match to this sequence.
# If reads are missing from the mapping object, use this sequence for the mapping.

updateLog(paste0('Leader sequences determined for ', sprintf("%.2f%%", (n_distinct(m$id)/n_distinct(reads$readID))*100), ' of unqiue read pairs.'))

if(opt$prepReads_forceAnchorReadStartSeq & 'anchorReadStartSeq' %in% names(reads)){
  a <- subset(reads, readID %in% m$id)
  b <- subset(reads, ! readID %in% m$id)
  
  if(nrow(b) > 0){
    updateLog('Some reads failed to create a leader sequence map; using anchorReadStartSeq sequences from sample data for missing entries.')
  
    d <- distinct(select(reads, uniqueSample, anchorReadStartSeq))
    
    m2 <- tibble(id = b$readID, leaderMapping.qStart = 1, uniqueSample = b$uniqueSample) %>% left_join(d, by = 'uniqueSample')
    m2$leaderMapping.qEnd <- nchar(m2$anchorReadStartSeq)
    m2$leaderSeqMap <- paste0('1..', m2$leaderMapping.qEnd, '[00+00]')
    m2 <- dplyr::select(m2, -uniqueSample, -anchorReadStartSeq)
    m <- bind_rows(m, m2)
  }
}

reads <- subset(reads, readID %in% m$id)

if(nrow(reads) == 0) quitOnErorr('Error - no reads remaining after filtering for leader sequence hits.')

reads <- left_join(reads, dplyr::select(m, id, leaderMapping.qStart, leaderMapping.qEnd), by = c('readID' = 'id'))
reads$leaderSeq = substr(reads$anchorReadSeq, 1, reads$leaderMapping.qEnd)

saveRDS(m, file.path(opt$outputDir, opt$prepReads_outputDir, 'leaderSeqMaps.rds'), compress = opt$compressDataFiles)

updateLog('Removing indentified leader sequences from anchor reads.')

reads$anchorReadSeq2 <- substr(reads$anchorReadSeq, reads$leaderMapping.qEnd+1, nchar(reads$anchorReadSeq))
reads <- dplyr::select(reads, -leaderMapping.qStart, -leaderMapping.qEnd, -anchorReadSeq)
reads <- dplyr::rename(reads, anchorReadSeq = anchorReadSeq2)

nReadsPreFilter <- n_distinct(reads$readID)
reads <- dplyr::filter(reads, nchar(anchorReadSeq) >= opt$prepReads_minAnchorReadLength)

if(nrow(reads) == 0) quitOnErorr('Error - no reads remaining after trimming for min anchor read length.')

nReadsPostFilter <- n_distinct(reads$readID)

msg <- paste0(sprintf("%.2f%%", (1 - n_distinct(reads$readID) / nReadsPreFilter)*100), 
               ' of reads removed because they were less than ', opt$prepReads_minAnchorReadLength, ' NTs after trimming.')

updateLog(msg)

# Create 15 NT adapter sequences by taking the reverse complement of identified leader sequences.
updateLog('Creating adrift read over-reading trim sequences.')

reads$adriftReadTrimSeq <- as.character(reverseComplement(DNAStringSet(substr(reads$leaderSeq, nchar(reads$leaderSeq)-14, nchar(reads$leaderSeq)))))

closeAllConnections()
cluster <- makeCluster(opt$prepReads_CPUs)
clusterSetRNGStream(cluster, 1)
clusterExport(cluster, c('opt'))

updateLog('Triming adrift read over-reading.')

nReadsPreFilter <- n_distinct(reads$readID)

# Trim adrift reads with adapter sequences determined by the RC of their anchor read leader sequences.
reads <- data.table::rbindlist(parLapply(cluster, split(reads, dplyr::ntile(1:nrow(reads), opt$prepReads_CPUs)), function(x){
  source(file.path(opt$softwareDir, 'lib.R'))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(Biostrings))
  
  data.table::rbindlist(lapply(split(x, x$adriftReadTrimSeq), function(y){
    f <- file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp',  tmpFile())
    o <- DNAStringSet(y$adriftReadSeq)
    names(o) <- y$readID
    Biostrings::writeXStringSet(o, f)
    
    # Use cut adapt to trim anchor read over-reading.
    system(paste0('cutadapt -e ', opt$prepReads_cutAdaptErrorRate, ' -a ', y$adriftReadTrimSeq[1], ' ',
                  f, ' > ', paste0(f, '.cutAdapt')), ignore.stderr = TRUE)
    
    waitForFile(paste0(f, '.cutAdapt'))
    
    t <- readDNAStringSet(paste0(f, '.cutAdapt'))
    invisible(file.remove(f, paste0(f, '.cutAdapt')))
    t <- t[width(t) >= opt$prepReads_minAdriftReadLength]
    if(length(t) == 0) return(data.table())
    
    # Add back trimmed adrift read sequences.
    o <- subset(y, readID %in% names(t))
    trimmed <- data.table(readID = names(t), adriftReadSeq2 = as.character(t))
    o <- left_join(o, trimmed, by = 'readID')
    
    data.table(dplyr::select(o, -adriftReadTrimSeq))
  }))
}))

if(nrow(reads) == 0) quitOnErorr('Error - no reads remaining after adrift read over-reading trimming.')

updateLog(paste0(sprintf("%.2f%%", (1-n_distinct(reads$readID)/nReadsPreFilter)*100), 
                 ' of reads removed because their trimmed lengths were less than ', opt$prepReads_minAdriftReadLength, ' NTs.'))

reads <- dplyr::select(reads, -adriftReadSeq) %>% 
         dplyr::rename(adriftReadSeq = adriftReadSeq2)

invisible(file.remove(list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'), full.names = TRUE)))

reads <- bind_rows(lapply(split(reads, reads$vectorFastaFile), function(x){
           x$vectorFastaFile <- lpe(x$vectorFastaFile[1])
           x
          }))

if('leaderSeqHMM' %in% names(reads)){
  reads <- bind_rows(lapply(split(reads, reads$leaderSeqHMM), function(x){
             x$leaderSeqHMM <- lpe(x$leaderSeqHMM[1])
             x
           }))
}

unlink(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs'), recursive = TRUE) 

if('quickFilterStartPos' %in% names(reads)) reads$quickFilterStartPos <- NULL

# If UMI processing is disabled, set all random UMI sequences to poly-A and updated
# adrift reads by replacing the random UMI sequences with poly-A.

if(! opt$processAdriftReadLinkerUMIs){
  reads$adriftReadRandomID <- 'AAAAAAAAAAAA'
} else {
  if('n' %in% names(reads)) reads$n <- NULL
  
  if(n_distinct(reads$adriftReadRandomID) == 1){
    tab <- data.frame('Var1' = unique(reads$adriftReadRandomID), 'Freq' = nrow(reads)) 
  } else{
    tab <- data.frame(sort(table(reads$adriftReadRandomID), decreasing = TRUE))
  }
  
  tab$Var1 <- as.character(tab$Var1)
  names(tab) <- c('adriftReadRandomID', 'n')
  tab <- tab[! duplicated(tab$adriftReadRandomID),]
  reads <- left_join(reads, tab, by = 'adriftReadRandomID')
  reads <- arrange(reads, desc(n))
  
  # UMIs read three or more times are considered abundant and are used as truth for correcting less read UMIs.
  
  a <- reads[reads$n >= 3,]
  u <- unique(a$adriftReadRandomID)
  b <- reads[reads$n < 3,]
  
  clusterExport(cluster, c('u'))
  b <- data.table(b)
  
  b <- rbindlist(parLapply(cluster, split(b, ntile(1:nrow(b), opt$demultiplex_CPUs)), function(k){
         suppressPackageStartupMessages(library(dplyr))
         suppressPackageStartupMessages(library(data.table))
         suppressPackageStartupMessages(library(stringdist))
    
         t <- n_distinct(k$adriftReadRandomID)
         n <- 1
         f <- paste0('UMI_correction_log_', paste0(paste0(stringi::stri_rand_strings(8, 1, '[A-Za-z0-9_\\-\\!\\%]'), collapse = '')))
    
         logFile <- file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp', f)
    
         rbindlist(lapply(split(k, k$adriftReadRandomID), function(x){
           if(n %% 100 == 0) write(sprintf("%.2f%%", (n/t)*100), file = logFile)
           n <<- n + 1
      
           d <- stringdist(x$adriftReadRandomID[1], u, nthread = 1)
      
           if(any(which(d == 1))){
             o <- u[which(d == 1)]
        
             if(length(o) == 1){
               x$adriftReadRandomID <- o
             } else {
               # If an UMI is 1 away from two or more unique UMIs then we let it go.
              x$adriftReadRandomID <- 'x'
             }
          }
          x
         }))
  }))
  
  b <- b[b$adriftReadRandomID != 'x']
  reads <- rbindlist(list(a, b))
  
  tab <- lazy_dt(reads) %>% 
         group_by(adriftReadRandomID) %>%
         summarise(nSamples = n_distinct(uniqueSample)) %>%
         ungroup() %>%
         as.data.table()
  
  a <- reads[reads$adriftReadRandomID %in% tab[tab$nSamples == 1]$adriftReadRandomID]
  b <- reads[reads$adriftReadRandomID %in% tab[tab$nSamples > 1]$adriftReadRandomID]
  
  if(nrow(b) > 0){
    b <- rbindlist(lapply(split(b, b$adriftReadRandomID), function(x){
           tab <- sort(table(x$uniqueSample), decreasing = TRUE)
    
           if(tab[1] > tab[2]){
             x$uniqueSample <- names(tab)[1]
           } else {
             x$uniqueSample <- 'x'
          }
          x
          }))
  
     b <- b[b$uniqueSample != 'x']
  }
  
  reads <- rbindlist(list(a, b))
}

stopCluster(cluster)

saveRDS(reads, file.path(opt$outputDir, opt$prepReads_outputDir, 'reads.rds'), compress = opt$compressDataFiles)

updateLog('prepReads completed.')

if(any(! incomingSamples %in% reads$uniqueSample) & opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()

q(save = 'no', status = 0, runLast = FALSE) 
