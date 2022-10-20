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

invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$prepReads_outputDir))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs'))

reads <-  readRDS(file.path(opt$outputDir, opt$prepReads_readsTable))

samples <- loadSamples()

if('vectorFastaFile' %in% names(samples)){
  samples$vectorFastaFile <- file.path(opt$softwareDir, 'data', 'vectors', samples$vectorFastaFile)

  if(! all(sapply(unique(samples$vectorFastaFile), file.exists))){
    write(c(paste(now(), "Error - one or more vector FASTA files could not be found in AAVengeR's data/vectors directory")), file = file.path(opt$outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE) 
  }
} else {
  samples$vectorFastaFile <- NA
}


if('leaderSeqHMM' %in% names(samples)){
  samples$leaderSeqHMM <- file.path(opt$softwareDir, 'data', 'hmms', samples$leaderSeqHMM)
  
  if(! all(sapply(unique(samples$leaderSeqHMM), file.exists))){
    write(c(paste(now(), "Error - one or more leader sequence HMM files could not be found in AAVengeR's data/hmms directory")), file = file.path(opt$outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE) 
  }
}


cluster <- makeCluster(opt$prepReads_CPUs)
clusterExport(cluster, c('samples', 'tmpFile', 'waitForFile', 'opt', 'lpe', 'blastReads'))

# Trim anchor read over-reading and trim off adrift read linkers.
reads <- data.table::rbindlist(parLapply(cluster, split(reads, dplyr::ntile(1:nrow(reads), opt$prepReads_CPUs)), function(x){
           source(file.path(opt$softwareDir, 'lib.R'))
           library(dplyr)
           library(data.table)
           library(Biostrings)
  
           data.table::rbindlist(lapply(split(x, x$adriftReadTrimSeq), function(y){
             f <- file.path(opt$outputDir, 'tmp',  tmpFile())
             o <- DNAStringSet(y$anchorReadSeq)
             names(o) <- y$readID
             Biostrings::writeXStringSet(o, f)
         
            # Use cut adapt to trim anchor read over-reading.
            system(paste0(opt$command_cutadapt, ' -e 0.15 -a ', y$adriftReadTrimSeq[1], ' --overlap 2 ',
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

reads$anchorReadSeq <- NULL
reads$adriftReadSeq <- NULL
reads <- dplyr::rename(reads, anchorReadSeq = anchorReadSeq2, adriftReadSeq = adriftReadSeq2)


reads$i <- group_by(reads, uniqueSample, anchorReadSeq, adriftReadSeq) %>% group_indices() 
reads <- group_by(reads, i) %>% mutate(n = n()) %>% ungroup()

a <- subset(reads, n == 1)
b <- subset(reads, n > 1)

c <- rbindlist(lapply(split(b, b$i), function(x){
       data.table(id = x[1,]$readID, n = nrow(x) - 1, id2 = x[2:nrow(x),]$readID)
     }))

saveRDS(c, file.path(opt$outputDir, opt$prepReads_outputDir, 'duplicateReads.rds'))

reads <- data.table(bind_rows(a, subset(b, ! readID %in% c$id2)) %>% dplyr::select(-i, -n))

invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))
rm(a, b, c)
gc()


mappings <- tibble()
vectorHits <- tibble()

cluster <- makeCluster(opt$prepReads_CPUs)
clusterExport(cluster, c('samples', 'tmpFile', 'waitForFile', 'opt', 'lpe', 'blastReads'))

  
# Align the ends of anchor reads to the vector to identify vector reads which should be removed.

if(opt$prepReads_excludeAnchorReadVectorHits | opt$prepReads_excludeAdriftReadVectorHits){
  
  vectorHits <- rbindlist(lapply(split(reads, reads$vectorFastaFile), function(x){
            invisible(file.remove(list.files(file.path(opt$outputDir, opt$vectorFilter_outputDir, 'dbs'), full.names = TRUE)))
  
            system(paste0(opt$command_makeblastdb, ' -in ', file.path(opt$softwareDir, 'data', 'vectors', x$vectorFastaFile[1]), 
                          ' -dbtype nucl -out ', file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
            waitForFile(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd.nin'))
  
            rbindlist(parLapply(cluster, split(x, ntile(1:nrow(x), opt$prepReads_CPUs)), function(y){
            #rbindlist(lapply(split(x, ntile(1:nrow(x), opt$prepReads_CPUs)), function(y){   
              library(Biostrings)
              library(data.table)
              library(dplyr)
              library(lubridate)
              
              b1 <- data.table()
              b2 <- data.table()
              
              if(opt$prepReads_excludeAnchorReadVectorHits){
                s <- DNAStringSet(y$anchorReadSeq)
                names(s) <- y$readID
                s <- subseq(s, (width(s) - opt$prepReads_vectorAlignmentTestLength) + 1 , width(s))
              
                b1 <- data.table(blastReads(s))
                if(nrow(b1) > 0){
                    b1$source <- 'anchor'
                    b1$alignmentLength <- b1$qend - b1$qstart + 1
                    b1 <- dplyr::filter(b1, pident >= opt$prepReads_vectorAlignmentTest_minPercentID, alignmentLength >= floor(opt$prepReads_vectorAlignmentTestLength * (opt$prepReads_vectorAlignmentTest_minPercentCoverage/100)), gapopen <= 1)
                }
              }
              
              if(opt$prepReads_excludeAdriftReadVectorHits){
                s <- DNAStringSet(y$adriftReadSeq)
                names(s) <- y$readID
                s <- subseq(s, (width(s) - opt$prepReads_vectorAlignmentTestLength) + 1 , width(s))
                
                b2 <- data.table(blastReads(s))
                if(nrow(b2) > 0){
                  b2$source <- 'adrift'
                  b2$alignmentLength <- b2$qend - b2$qstart + 1
                  b2 <- dplyr::filter(b2, pident >= opt$prepReads_vectorAlignmentTest_minPercentID, alignmentLength >= floor(opt$prepReads_vectorAlignmentTestLength * (opt$prepReads_vectorAlignmentTest_minPercentCoverage/100)), gapopen <= 1)
                }
              }
              
              b <- rbindlist(list(b1, b2))
              
              if(nrow(b) > 0){
                b$start  <- ifelse(b$sstart > b$send, b$send, b$sstart)
                b$end    <- ifelse(b$sstart > b$send, b$sstart, b$send)
                b$strand <- ifelse(b$sstart > b$send, '-', '+')
                b <- group_by(b, qname) %>% dplyr::top_n(1, wt = bitscore) %>% dplyr::arrange(desc(strand)) %>% dplyr::slice(1) %>% ungroup() %>% data.table()
              }
              
              b
            }))
          }))
  
  invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))
}

saveRDS(vectorHits, file.path(opt$outputDir, opt$prepReads_outputDir, 'vectorHits.rds'))

reads <- subset(reads, ! readID %in% vectorHits$qname)

if(! 'leaderSeqHMM' %in% names(samples)){
    # Now align the full anchor reads to the vector excluding those in vectorHits.
    write(c(paste(now(), '   Aligning full length anchor reads to the vector sequence.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
    vectorHits2 <- rbindlist(lapply(split(reads, reads$vectorFastaFile), function(x){
      invisible(file.remove(list.files(file.path(opt$outputDir, opt$vectorFilter_outputDir, 'dbs'), full.names = TRUE)))
    
      system(paste0(opt$command_makeblastdb, ' -in ', file.path(opt$softwareDir, 'data', 'vectors', x$vectorFastaFile[1]), ' -dbtype nucl -out ', file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
      waitForFile(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd.nin'))
    
      rbindlist(parLapply(cluster,split(x, ntile(1:nrow(x), opt$prepReads_CPUs)), function(y){
        library(Biostrings)
        library(data.table)
        library(dplyr)
        
        s <- DNAStringSet(y$anchorReadSeq)
        names(s) <- y$readID
        
        b <- data.table(blastReads(s))
        if(nrow(b) > 0){
          b$alignmentLength <- b$qend - b$qstart + 1
          b <- subset(b, pident >= opt$prepReads_minAlignmentPercentID & 
                         alignmentLength >= opt$prepReads_minAlignmentLength &
                         gapopen <= 1)
        }
        
        b
      }))
    }))
    
    write(c(paste(now(), '   Creating read maps from local alignments to the vector file.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
    
    # Limit reads to those with hits to the vector since we expect the start of anchor reads to begin with vector sequences.
    reads <- subset(reads, readID %in% vectorHits2$qname)
    
    # Split reads into CPU chunks without dividing blast results for a read between chunks.
    vectorHits2$i <- group_by(vectorHits2, qname) %>% group_indices() 
    o <- tibble(n = 1:n_distinct(vectorHits2$i))
    o$i2 <- ntile(1:nrow(o), opt$prepReads_CPUs)
    vectorHits2 <- left_join(vectorHits2, o, by = c('i' = 'n'))
    
    m <- rbindlist(parLapply(cluster, split(vectorHits2, vectorHits2$i2), function(b){
           library(data.table)
           library(dplyr)
      
           if(! opt$prepReads_buildReadMaps_blastReconstruction){
             b2 <- filter(b, qstart <= opt$prepReads_buildReadMaps_minMapStartPostion) %>%
                          group_by(qname) %>% 
                          slice_min(evalue, n = 1) %>% 
                          dplyr::slice(1) %>% 
                          ungroup()
        
              mappings <- data.table(id = b2$qname, leaderMapping.qStart = 1, leaderMapping.qEnd = b2$qend,
                                 leaderMapping.sStart = NA, leaderMapping.sEnd = NA)
           } else {
             library(GenomicRanges)
          
             mappings <- rbindlist(lapply(split(b, b$qname), function(a){
               g <- makeGRangesFromDataFrame(a, ignore.strand = TRUE, seqnames.field = 'sseqid',  start.field = 'qstart', end.field = 'qend')
               if(length(g) == 0) return(tibble())
               g <- GenomicRanges::reduce(g, min.gapwidth = 4, ignore.strand = TRUE) # Allow merging if ranges separated by <= 3 NTs.
               g <- g[start(g) <= opt$prepReads_buildReadMaps_minMapStartPostion]
               if(length(g) == 0) return(tibble())
               g <- g[width(g) == max(width(g))][1]
               return(data.table(id = a$qname[1], leaderMapping.qStart = 1, leaderMapping.qEnd = end(g),
                            leaderMapping.sStart = NA, leaderMapping.sEnd = NA))
            }))
           }
      
       mappings
    }))
    
} else {
  write(c(paste(now(), '   Using leader sequence HMM to define mappings.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
  # (!) UPDATE ME (!)
  
  
  hmmResults <- bind_rows(lapply(split(d, 1:nrow(d)), function(x){ 
                 source(file.path(opt$softwareDir, 'lib.R'))
                 library(dplyr)
                 library(Biostrings)
                 message(x$file)
                 captureLTRseqsLentiHMM(readDNAStringSet(x$file), subset(samples, uniqueSample == x$uniqueSample)$leaderSeqHMM)
              }))
 
 if(nrow(hmmResults) > 0){
   m <- tibble(id = hmmResults$id, leaderMapping.qStart = 1, leaderMapping.qEnd = nchar(hmmResults$LTRseq), 
                      leaderMapping.sStart = NA, leaderMapping.sEnd = NA)
 } else {
   write(c(paste(now(), "   Error - no reads matched the HMMs.")), file = file.path(opt$outputDir, 'log'), append = TRUE)
   q(save = 'no', status = 1, runLast = FALSE) 
 }
}

stopCluster(cluster)


reads <- subset(reads, readID %in% m$id)
reads <- left_join(reads, dplyr::select(m, id, leaderMapping.qStart, leaderMapping.qEnd), by = c('readID' = 'id'))
reads$leaderSeq = substr(reads$anchorReadSeq, 1, reads$leaderMapping.qEnd)


reads$anchorReadSeq2 <- substr(reads$anchorReadSeq, reads$leaderMapping.qEnd+1, nchar(reads$anchorReadSeq))
reads <- dplyr::select(reads, -leaderMapping.qStart, -leaderMapping.qEnd, -anchorReadSeq)
reads <- dplyr::rename(reads, anchorReadSeq = anchorReadSeq2)

reads <- dplyr::filter(reads, nchar(anchorReadSeq) >= opt$prepReads_minAnchorReadLength)

# Create 15 NT adapter sequences. 
reads$adriftReadTrimSeq <- as.character(reverseComplement(DNAStringSet(substr(reads$leaderSeq, nchar(reads$leaderSeq)-14, nchar(reads$leaderSeq)))))


cluster <- makeCluster(opt$prepReads_CPUs)
clusterExport(cluster, c('opt'))

reads <- data.table::rbindlist(parLapply(cluster, split(reads, dplyr::ntile(1:nrow(reads), opt$prepReads_CPUs)), function(x){
  source(file.path(opt$softwareDir, 'lib.R'))
  library(dplyr)
  library(data.table)
  library(Biostrings)
  
  data.table::rbindlist(lapply(split(x, x$adriftReadTrimSeq), function(y){
    f <- file.path(opt$outputDir, 'tmp',  tmpFile())
    o <- DNAStringSet(y$adriftReadSeq)
    names(o) <- y$readID
    Biostrings::writeXStringSet(o, f)
    
    # Use cut adapt to trim anchor read over-reading.
    system(paste0(opt$command_cutadapt, ' -e 0.15 -a ', y$adriftReadTrimSeq[1], ' --overlap 2 ',
                  f, ' > ', paste0(f, '.cutAdapt')), ignore.stderr = TRUE)
    
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

reads <- dplyr::select(reads, -adriftReadSeq) %>% 
         dplyr::rename(adriftReadSeq = adriftReadSeq2) %>%
         dplyr::filter(nchar(adriftReadSeq) >= opt$prepReads_minAdriftReadLength)

invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

saveRDS(reads, file.path(opt$outputDir, opt$prepReads_outputDir, 'reads.rds'), compress = FALSE)

q(save = 'no', status = 0, runLast = FALSE) 
