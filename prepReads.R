# John K. Everett, PhD
# AAVengeR/prepReads.R
#
# This module prepares reads for alignment to a reference genome by removing duplicate
# read pairs, trimming adapter sequences, removing reads with high sequence homology to 
# the vector and selecting anchor reads with the expected leader sequences. 

library(ShortRead)
library(dplyr)
library(parallel)
library(lubridate)
library(Biostrings)
library(data.table)
options(stringsAsFactors = FALSE)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)

invisible(file.remove(list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'), full.names = TRUE)))

source(file.path(opt$softwareDir, 'lib.R'))
setMissingOptions()
setOptimalParameters()

dir.create(file.path(opt$outputDir, opt$prepReads_outputDir))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs'))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'))

write(c(paste(now(), '   Reading in demultiplexed reads.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = FALSE)
reads <- readRDS(file.path(opt$outputDir, opt$prepReads_readsTable))

incomingSamples <- unique(reads$uniqueSample)

cluster <- makeCluster(opt$prepReads_CPUs)
clusterExport(cluster, 'opt')

# Trim anchor read over-reading with cutadapt using the RC of the common linker in adrift reads.
# The trim sequence is created by demultiplex.R.

write(c(paste(now(), '   Trimming anchor read over-reading.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)

reads <- data.table::rbindlist(parLapply(cluster, split(reads, dplyr::ntile(1:nrow(reads), opt$prepReads_CPUs)), function(x){
           source(file.path(opt$softwareDir, 'lib.R'))
           library(dplyr)
           library(data.table)
           library(Biostrings)
  
           data.table::rbindlist(lapply(split(x, x$adriftReadTrimSeq), function(y){
             f <- file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp',  tmpFile())
             o <- DNAStringSet(y$anchorReadSeq)
             names(o) <- y$readID
             Biostrings::writeXStringSet(o, f)
         
            # Use cut adapt to trim anchor read over-reading.
            system(paste0('cutadapt -e 0.15 -a ', y$adriftReadTrimSeq[1], ' ',
                       f, ' > ', paste0(f, '.cutAdapt')), ignore.stderr = TRUE)
         
            t <- readDNAStringSet(paste0(f, '.cutAdapt'))
            invisible(file.remove(f, paste0(f, '.cutAdapt')))
            
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

if(nrow(reads) == 0){
  write(c(paste(now(), "Error - no reads remaining after trimming.")), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
  if(opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()
  q(save = 'no', status = 1, runLast = FALSE) 
}

trimReport <- mutate(reads, sample = sub('~\\d+$', '', uniqueSample)) %>% 
              group_by(sample) %>% 
              summarise(reads = n_distinct(readID), percentTrimmed = sprintf("%.2f%%", (sum(anchorReadSeq != anchorReadSeq2)/n())*100)) %>% 
              ungroup()
                                                    
trimReport$sample <- paste0('                       ', trimReport$sample)
write(paste0(now(),         '    Anchor reads trimmed with over-read adapter (sample, reads, percent trimmed):'), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
readr::write_tsv(trimReport, file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE, col_names = FALSE)

reads$anchorReadSeq <- NULL
reads$adriftReadSeq <- NULL
reads <- dplyr::rename(reads, anchorReadSeq = anchorReadSeq2, adriftReadSeq = adriftReadSeq2)


# Identify and remove identical read pairs.
write(c(paste(now(), '   Identifying duplicate read pairs.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
reads$i <- group_by(reads, uniqueSample, adriftReadRandomID, anchorReadSeq, adriftReadSeq) %>% group_indices() 
reads <- group_by(reads, i) %>% mutate(n = n()) %>% ungroup()

a <- subset(reads, n == 1) # Non-duplicated read pairs.
b <- subset(reads, n > 1)  # Duplicated read pairs.

write(c(paste(now(), '   Removing duplicate read pairs.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)

# Create a table of duplicate read pairs where one is chosen (id) to move forward and the others are logged (id2).
c <- group_by(b, i) %>%
      summarise(id = readID[1], n = n() - 1, id2 = list(readID[2:n()])) %>%
      ungroup() %>% select(-i) %>% tidyr::unnest(id2)

saveRDS(c, file.path(opt$outputDir, opt$prepReads_outputDir, 'duplicateReads.rds'))


# Exclude duplicate read pairs.
reads <- data.table(bind_rows(a, subset(b, ! readID %in% c$id2)) %>% dplyr::select(-i, -n))
saveRDS(reads, file.path(opt$outputDir, opt$prepReads_outputDir, 'uniqueReadPairs.rds'))


# Clean up.
invisible(file.remove(list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'), full.names = TRUE)))
rm(a, b, c)
gc()

mappings <- tibble()
vectorHits <- tibble()

# Align the ends of anchor reads to the vector to identify vector reads which should be removed.

reads$vectorFastaFile <- file.path(opt$softwareDir, 'data', 'vectors', reads$vectorFastaFile)

# save.image('~/CART_12418_prepReadsDev.RData')



alignReadEndsToVector <- function(y){
  library(Biostrings)
  library(data.table)
  library(dplyr)
  library(lubridate)
  
  s <- DNAStringSet(y$testSeq)
  names(s) <- y$readID
  
  b <- blastReads(s, tmpDirPath = file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'), dbPath = file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd'))
  
  if(nrow(b) > 0){
    b$alignmentLength <- b$qend - b$qstart + 1
    b <- dplyr::filter(b, pident >= opt$prepReads_vectorAlignmentTest_minPercentID, alignmentLength >= floor(opt$prepReads_vectorAlignmentTestLength * (opt$prepReads_vectorAlignmentTest_minPercentCoverage/100)), gapopen <= 1)
  }
  
  if(nrow(b) > 0){
    b <- left_join(b, data.table(qname = names(s), testSeq = as.character(s)), by = 'qname')
    b$start  <- ifelse(b$sstart > b$send, b$send, b$sstart)
    b$end    <- ifelse(b$sstart > b$send, b$sstart, b$send)
    b$strand <- ifelse(b$sstart > b$send, '-', '+')
    b <- group_by(b, qname) %>% dplyr::top_n(1, wt = bitscore) %>% dplyr::arrange(desc(strand)) %>% dplyr::slice(1) %>% ungroup()
  }
  
  b
}


cluster <- makeCluster(opt$prepReads_CPUs)
clusterExport(cluster, c('tmpFile', 'waitForFile', 'opt', 'lpe', 'blastReads', 'alignReadEndsToVector'))


if(opt$prepReads_excludeAnchorReadVectorHits | opt$prepReads_excludeAdriftReadVectorHits){
  
  write(c(paste(now(), '   Aligning ends of anchor reads to vector sequences.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
  
  vectorHits <- bind_rows(lapply(split(reads, reads$vectorFastaFile), function(x){
            invisible(file.remove(list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs'), full.names = TRUE)))
  
            # Mask lowercase letters in vector file because this is the means by which users signal repeat sequences.
            seq <- readLines(x$vectorFastaFile[1])
            seq <- gsub('[atcg]', 'N', seq)
            f <- tmpFile()
            
            write(seq, file = file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp', f))
            
            system(paste0('makeblastdb -in ', file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp', f), 
                          ' -dbtype nucl -out ', file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
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
              
              return(distinct(left_join(r, select(testSeqs, readID, i), by = c('qname' = 'readID')) %>%
                     left_join(select(testSeqs, readID, i), by = 'i') %>% 
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
saveRDS(vectorHits, file.path(opt$outputDir, opt$prepReads_outputDir, 'vectorHits.rds'))

if(nrow(vectorHits) > 0){
  readsPerSample <- mutate(reads, sample = sub('~\\d+$', '', uniqueSample)) %>%
                    group_by(sample) %>%
                    summarise(nReads = n_distinct(readID)) %>%
                    ungroup()

  vectorHitsReport <- left_join(vectorHits, select(reads, uniqueSample, readID), by = c('qname' = 'readID')) %>%
                      mutate(sample = sub('~\\d+$', '', uniqueSample)) %>%
                      left_join(readsPerSample, by = 'sample') %>%
                      group_by(sample) %>%  
                      summarise(reads = n_distinct(qname), 
                                percentSampleReads = sprintf("%.2f%%", (n_distinct(qname)/nReads[1])*100)) %>% 
                      ungroup()
  
  vectorHitsReport$sample <- paste0('                       ', vectorHitsReport$sample)
  write(paste0(now(),         '    Anchor reads aligning to the vector (sample, reads, percent sample reads):'), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
  readr::write_tsv(vectorHitsReport, file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE, col_names = FALSE)
} else {
  write(c(paste(now(), '   No anchor reads aligned to the vector.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
}  


nReadsPreFilter <- n_distinct(reads$readID)
reads <- subset(reads, ! readID %in% vectorHits$qname)
nReadsPostFilter <- n_distinct(reads$readID)
write(c(paste0(now(), '    ', sprintf("%.2f%%", 100 - (nReadsPostFilter/nReadsPreFilter)*100), ' unique reads removed because they aligned to the vector.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)


if(! 'leaderSeqHMM' %in% names(reads)){
    # Now align the full anchor reads to the vector excluding those in vectorHits.
    write(c(paste(now(), '   Aligning full anchor reads to vector sequences.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
  
    cluster <- makeCluster(opt$prepReads_CPUs)
    clusterExport(cluster, c('tmpFile', 'waitForFile', 'opt', 'lpe', 'blastReads', 'blast2rearangements_worker'))
  
    vectorHits2 <- rbindlist(lapply(split(reads, reads$vectorFastaFile), function(x){
      invisible(file.remove(list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs'), full.names = TRUE)))
    
      # No masking here here
      
      system(paste0('makeblastdb -in ', x$vectorFastaFile[1], ' -dbtype nucl -out ', file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
      waitForFile(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd.nin'))
    
      x$i <- group_by(x, anchorReadSeq) %>% group_indices() 
      x2 <- x[! duplicated(x$i),] 
      
      r <- bind_rows(parLapply(cluster,split(x2, ntile(1:nrow(x2), opt$prepReads_CPUs)), function(y){
      #r <- bind_rows(lapply(split(x2, ntile(1:nrow(x2), opt$prepReads_CPUs)), function(y){
             library(Biostrings)
             library(data.table)
             library(dplyr)
        
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
    
    if(nrow(vectorHits2) == 0){
      write(c(paste(now(), "Error - reads aligned to the vector for leader seq generation.")), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
      if(opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()
      q(save = 'no', status = 1, runLast = FALSE) 
    }
    
    write(c(paste(now(), '   Creating read maps from local alignments to the vector file.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
    
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
    
    parallel::stopCluster(cluster)
    
} else {
  write(c(paste(now(), '   Using leader sequence HMM to define mappings.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
  
  cluster <- parallel::makeCluster(opt$prepReads_CPUs)
  clusterExport(cluster, c('opt'))
  
  hmmResults <- rbindlist(parLapply(cluster, split(reads, reads$uniqueSample), function(x){
  #hmmResults <- rbindlist(lapply(split(reads, reads$uniqueSample), function(x){
    library(Biostrings)
    library(dplyr)
    source(file.path(opt$softwareDir, 'lib.R'))
    
    seqs <- DNAStringSet(x$anchorReadSeq)
    names(seqs) <- x$readID
    captureHMMleaderSeq(seqs, x$leaderSeqHMM[1], tmpDirPath = file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'))
  }))
  
  parallel::stopCluster(cluster)
 
  if(nrow(hmmResults) > 0){
    m <- tibble(id = hmmResults$id, leaderMapping.qStart = 1, leaderMapping.qEnd = nchar(hmmResults$LTRseq), leaderSeqMap = NA)
  } else {
    write(c(paste(now(), "   Error - no reads matched the HMMs.")), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
    if(opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()
    q(save = 'no', status = 1, runLast = FALSE) 
  }
}

# If anchorReadStartSeq was included in the sample data file then reads must of had a decent match to this sequence.
# If reads are missing from the mapping object, use this sequence for the mapping.

write(c(paste0(now(), '    Leader sequences determined for ', sprintf("%.2f%%", (n_distinct(m$id)/n_distinct(reads$readID))*100), 
              ' of unqiue read pairs.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)

if('anchorReadStartSeq' %in% names(reads)){
  a <- subset(reads, readID %in% m$id)
  b <- subset(reads, ! readID %in% m$id)
  
  if(nrow(b) > 0){
    write(c(paste(now(), '   Some reads failed to create a leader sequence map; using anchorReadStartSeq sequences from sample data for missing entries.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
    
    d <- distinct(select(reads, uniqueSample, anchorReadStartSeq))
    
    m2 <- tibble(id = b$readID, leaderMapping.qStart = 1, uniqueSample = b$uniqueSample) %>% left_join(d, by = 'uniqueSample')
    m2$leaderMapping.qEnd <- nchar(m2$anchorReadStartSeq)
    m2$leaderSeqMap <- paste0('1..', m2$leaderMapping.qEnd, '[00+00]')
    m2 <- dplyr::select(m2, -uniqueSample, -anchorReadStartSeq)
    m <- bind_rows(m, m2)
  }
}

reads <- subset(reads, readID %in% m$id)

if(nrow(reads) == 0){
  write(c(paste(now(), "Error - no reads remaining after trimming.")), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
  if(opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()
  q(save = 'no', status = 1, runLast = FALSE) 
}

reads <- left_join(reads, dplyr::select(m, id, leaderMapping.qStart, leaderMapping.qEnd), by = c('readID' = 'id'))
reads$leaderSeq = substr(reads$anchorReadSeq, 1, reads$leaderMapping.qEnd)

saveRDS(m, file.path(opt$outputDir, opt$prepReads_outputDir, 'leaderSeqMaps.rds'), compress = FALSE)

write(c(paste(now(), '   Removing indentified leader sequences from anchor reads.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)

reads$anchorReadSeq2 <- substr(reads$anchorReadSeq, reads$leaderMapping.qEnd+1, nchar(reads$anchorReadSeq))
reads <- dplyr::select(reads, -leaderMapping.qStart, -leaderMapping.qEnd, -anchorReadSeq)
reads <- dplyr::rename(reads, anchorReadSeq = anchorReadSeq2)

nReadsPreFilter <- n_distinct(reads$readID)
reads <- dplyr::filter(reads, nchar(anchorReadSeq) >= opt$prepReads_minAnchorReadLength)

if(nrow(reads) == 0){
  write(c(paste(now(), "Error - no reads remaining after trimming.")), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
  if(opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()
  q(save = 'no', status = 1, runLast = FALSE) 
}

nReadsPostFilter <- n_distinct(reads$readID)

write(c(paste0(now(), '    ', sprintf("%.2f%%", (1 - n_distinct(reads$readID) / nReadsPreFilter)*100), 
               ' of reads removed because they were less than ', opt$prepReads_minAnchorReadLength, 
               ' NTs after trimming.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)


# Create 15 NT adapter sequences by taking the reverse complement of identified leader sequences.
write(c(paste(now(), '   Creating adrift read over-reading trim sequences.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
reads$adriftReadTrimSeq <- as.character(reverseComplement(DNAStringSet(substr(reads$leaderSeq, nchar(reads$leaderSeq)-14, nchar(reads$leaderSeq)))))

closeAllConnections()
cluster <- makeCluster(opt$prepReads_CPUs)
clusterExport(cluster, c('opt'))

write(c(paste(now(), '   Triming adrift read over-reading.')), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)

nReadsPreFilter <- n_distinct(reads$readID)


cluster <- parallel::makeCluster(opt$prepReads_CPUs)
clusterExport(cluster, c('opt'))

# Trim adrift reads with adapter sequences determined by the RC of their anchor read leader sequences.
reads <- data.table::rbindlist(parLapply(cluster, split(reads, dplyr::ntile(1:nrow(reads), opt$prepReads_CPUs)), function(x){
  source(file.path(opt$softwareDir, 'lib.R'))
  library(dplyr)
  library(data.table)
  library(Biostrings)
  
  data.table::rbindlist(lapply(split(x, x$adriftReadTrimSeq), function(y){
    f <- file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp',  tmpFile())
    o <- DNAStringSet(y$adriftReadSeq)
    names(o) <- y$readID
    Biostrings::writeXStringSet(o, f)
    
    # Use cut adapt to trim anchor read over-reading.
    system(paste0('cutadapt -e 0.15 -a ', y$adriftReadTrimSeq[1], ' ',
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

stopCluster(cluster)

if(nrow(reads) == 0){
  write(c(paste(now(), "Error - no reads remaining after trimming.")), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
  if(opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()
  q(save = 'no', status = 1, runLast = FALSE) 
}

write(c(paste0(now(), paste0('    ', sprintf("%.2f%%", (1-n_distinct(reads$readID)/nReadsPreFilter)*100), 
                             ' of reads removed because their trimmed lengths were less than ', opt$prepReads_minAdriftReadLength, ' NTs.'))), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
  
trimReport <- mutate(reads, sample = sub('~\\d+$', '', uniqueSample)) %>% 
              group_by(sample) %>% 
              summarise(reads = n_distinct(readID), percentTrimmed = sprintf("%.2f%%", (sum(adriftReadSeq != adriftReadSeq2)/n())*100)) %>% 
              ungroup()

trimReport$sample <- paste0('                       ', trimReport$sample)
write(paste0(now(),         '    Adrift reads trimmed with over-read adapters from leader sequences, reads below min. length threshold post trim removed, (sample, reads, percent trimmed):'), file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE)
readr::write_tsv(trimReport, file = file.path(opt$outputDir, opt$prepReads_outputDir, 'log'), append = TRUE, col_names = FALSE)

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

saveRDS(reads, file.path(opt$outputDir, opt$prepReads_outputDir, 'reads.rds'), compress = TRUE)

if(any(! incomingSamples %in% reads$uniqueSample) & opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()

q(save = 'no', status = 0, runLast = FALSE) 
