#!/usr/bin/Rscript

# AAVengeR/anchorReadRearrangements.R
# John K. Everett, Ph.D.

for (p in c('ShortRead', 'parallel', 'dplyr', 'data.table', 'Biostrings')) suppressPackageStartupMessages(library(p, character.only = TRUE))

args <- commandArgs(trailingOnly=TRUE)

if(length(args) == 0) stop('Expected at least one command line argument')
source(file.path(yaml::read_yaml(args[1])$softwareDir, 'lib', 'main.R'))
opt <- startModule(args)

createOuputDir()
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir), showWarnings = FALSE)
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), showWarnings = FALSE)
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences'))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences'), showWarnings = FALSE)
invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), full.names = TRUE)))
invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences'), full.names = TRUE)))

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


# Split apart list parameters.
seeds <- as.integer(unlist(strsplit(opt$anchorReadRearrangements_seed, '\\s*,\\s*')))
rarifactionLevels <- as.integer(unlist(strsplit(opt$anchorReadRearrangements_rarifactionLevels, '\\s*,\\s*')))
windows <- as.integer(unlist(strsplit(opt$anchorReadRearrangements_windows, '\\s*,\\s*')))
windowsTypes <- unlist(strsplit(opt$anchorReadRearrangements_windowsTypes, '\\s*,\\s*'))

# Read in the reads table.
updateLog('Reading in demultiplexed reads.')
if(! file.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_inputFile))) quitOnErorr('Error - the input data file does not exist.')
reads <- readRDS(file.path(opt$outputDir, opt$anchorReadRearrangements_inputFile))
if(nrow(reads) == 0) quitOnErorr('Error - the input table contained no rows.')


# Test for raw reads, vector files, and 500mers in config file.
if(! all(unique(reads$vectorFastaFile) %in% names(opt$anchorReadRearrangements_expectedSeqs))) quitOnErorr('Error: one or more vector expected sequences are missing from the configuration file.')
if(! all(sapply(unique(reads$vectorFastaFile), function(x) file.exists(file.path(opt$softwareDir, 'data', 'vectors', x))))) quitOnErorr('Error: one or more vector sequences are missing from the AAVengeR data/vectors folder.')


updateLog('Reading raw sequencing data.')
updateMasterLog()


R1 <- readFastq(opt$demultiplex_adriftReadsFile)
R1@id <- BStringSet(sub('\\s+.+$', '', as.character(R1@id)))

R2 <- readFastq(opt$demultiplex_anchorReadsFile)
R2@id <- BStringSet(sub('\\s+.+$', '', as.character(R2@id)))

updateLog('Writing out sequencing data.')
updateMasterLog()

writeFastq(R1, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'adriftReadSeqs.fastq'), compress = FALSE, width = NA)
writeFastq(R2, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'anchorReadSeqs.fastq'), compress = FALSE, width = NA)

updateLog('Running PEAR.')
updateMasterLog()

system(paste0('pear',
              ' -q ', opt$anchorReadRearrangements_overlapQualTrimLevel,
              ' -t ', opt$anchorReadRearrangements_overlapTrimmedReadMinLen, 
              ' -f ', file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'adriftReadSeqs.fastq'), 
              ' -r ', file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'anchorReadSeqs.fastq'), 
              ' -o ', file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'merged'), 
              ' -p ', opt$anchorReadRearrangements_overlapMaxPval, 
              ' -v ', opt$anchorReadRearrangements_overlapMinReadOverlap, 
              ' -n ', opt$anchorReadRearrangements_minAnchorReadLength, 
              ' -j ', opt$anchorReadRearrangements_CPUs))

m <- readFastq(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'merged.assembled.fastq'))
m <- m[! grepl('NNN', as.character(m@sread))]
m <- m[as.character(m@id) %in% reads$readID]

o <- left_join(tibble(readID = as.character(m@id)), select(reads, readID, adriftLinkerSeqEnd), by = 'readID')

m <- narrow(m, (o$adriftLinkerSeqEnd + 1), width(m))
m <- reverseComplement(m)

# File clean up.
invisible(file.remove(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'merged.assembled.fastq')))
invisible(file.remove(c(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'anchorReadSeqs.fastq'),
                        file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'adriftReadSeqs.fastq'))))

updateLog('Processing sample replicates.')
updateMasterLog()

# Short circuit
#### m <- ShortRead::ShortReadQ(sread = DNAStringSet('ATG'),  id = BStringSet('read'), quality = BStringSet('???'))

reads <- bind_rows(lapply(split(reads, reads$uniqueSample), function(x){
  updateLog(paste0('Processing ', x$uniqueSample[1], '.'))
  updateLog(paste0('  ', ppNum(n_distinct(x$readID)), ' demultiplexed reads.'))
  
  paths <- opt[['anchorReadRearrangements_expectedSeqs']][[x$vectorFastaFile[1]]]
  testSeqs <- unique(substr(paths, 1, opt$anchorReadRearrangements_vectorLeaderSeqFilterLength))

  a <- subset(x, readID %in% as.character(m@id))
  b <- subset(x, ! readID %in% as.character(m@id))
  
  updateLog(paste0('  ', round((n_distinct(a$readID)/n_distinct(x$readID))*100, 2), '% of reads found in overlapped collection.'))
  
  if(nrow(a) > 0){
    mm <- m[as.character(m@id) %in% a$readID]
    
    tt <- sapply(testSeqs, function(xx) vcountPattern(xx, subseq(mm@sread, 1, opt$anchorReadRearrangements_vectorLeaderSeqFilterLength), max.mismatch = 1), simplify = FALSE)
    zz <- apply(bind_cols(tt), 1, function(row) any(row == 1))
    
    m <- m[zz]
    
    if(length(m) > 0){
      o <- tibble(readID = as.character(mm@id), anchorReadSeq = as.character(mm@sread))
      o <- o[match(a$readID, o$readID),]
      a$anchorReadSeq <- NULL
      a <- left_join(a, o, by = 'readID')
    }
  }
  
  x <- bind_rows(a, b)
  
  x <- x[nchar(x$anchorReadSeq) >= opt$anchorReadRearrangements_minAnchorReadLength,]
  
  updateLog(paste0('  ', ppNum(n_distinct(x$readID)), ' reads available after ', opt$anchorReadRearrangements_minAnchorReadLength, ' nt filter.'))
  
  if(nrow(x) == 0) return(tibble())

  tests <- as.data.frame(lapply(testSeqs, function(t){
    vcountPattern(t, DNAStringSet(substr(x$anchorReadSeq, 1, nchar(t))),
                  max.mismatch = opt$anchorReadRearrangements_vectorLeaderSeqFilterMaxMismatch) == 1
  }), fix.empty.names = FALSE)
  
  i <- apply(tests, 1, any)
  updateLog(paste0('  ', round((sum(i == TRUE) / length(i))*100, 2), '% of reads start with expected sequence(s).'))
  x <- x[i,]
  
  if(nrow(x) == 0) return(tibble())
  
  x$lastNTs <- substr(x$anchorReadSeq, nchar(x$anchorReadSeq) - (opt$anchorReadRearrangements_readEndAlignmentTestLength - 1), nchar(x$anchorReadSeq))
  x$i <- paste0('seq', 1:nrow(x))
  
  system(paste0('faToTwoBit ',
                file.path(opt$softwareDir, 'data', 'vectors', x$vectorFastaFile[1]), ' ',
                file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'TMP_DB.2bit')))
  
  d <- DNAStringSet(x$lastNTs)
  names(d) <- x$i
  
  writeXStringSet(d,  file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'TMP_QUERY.fasta'))
  
  system(paste0('blat ', file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'TMP_DB.2bit'), ' ',
                file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'TMP_QUERY.fasta'), ' ',
                file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'TMP_QUERY.psl'),
                ' -tileSize=6 -stepSize=3 -repMatch=3000 -out=psl -t=dna -q=dna -minScore=0 -minIdentity=0 -noHead -noTrimA'))
  
  b <- parseBLAToutput(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'TMP_QUERY.psl'))
  b <- subset(b, b$matches >= opt$anchorReadRearrangements_readEndAlignmentTestMinMatch)
  
  invisible(file.remove(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'TMP_DB.2bit'),
                        file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'TMP_QUERY.fasta'),
                        file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', 'TMP_QUERY.psl')))
  
  updateLog(paste0('  ', round((sum(x$i %in% b$qName) / nrow(x))*100, 2), '% reads have ends that align to the vector plasmid.'))
  updateMasterLog()
  
  if(nrow(b) == 0) return(tibble())
  
  x <- subset(x, i %in% b$qName)
  
  updateLog(paste0('  ', ppNum(n_distinct(x$readID)), ' reads available.'))
  
  tibble(readID = x$readID,
         uniqueSample = x$uniqueSample[1],
         vectorFastaFile = x$vectorFastaFile[1],
         seq = x$anchorReadSeq)
}))

rm(m); gc()


buildSeqModel <- function(seq, db, id, tmpDir){
  write(c('>seq', seq), file.path(tmpDir, paste0(id, '.fasta2')))
  
  system(paste0('blastn -num_threads 1 -penalty -4 -max_target_seqs 10000 -gapopen 10 -gapextend 5 -dust no -soft_masking false -word_size 5 -evalue 100 -outfmt 6 ',
                '-query ', file.path(tmpDir, paste0(id, '.fasta2')), ' ', 
                '-db ', db, ' ',
                '-out ', file.path(tmpDir, paste0(id, '.blast'))), ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  altStructSchema <- NA
  
  if(file.exists(file.path(tmpDir, paste0(id, '.blast')))){
    if(file.info(file.path(tmpDir, paste0(id, '.blast')))$size != 0){
      b <- read.table(file.path(tmpDir, paste0(id, '.blast')), sep = '\t', header = FALSE)
      names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
      b$alignmentLength <- b$qend - b$qstart + 1
      b$strand <- ifelse(b$sstart > b$send, '-', '+')
      b <- dplyr::filter(b, pident >= opt$anchorReadRearrangements_minAlignmentPercentID, alignmentLength >= opt$anchorReadRearrangements_minAlignmentLength)
      
      if(nrow(b) > 0){
        b$qlen <- nchar(seq)
        altStructSchema <- blast2rearangements(b, opt$anchorReadRearrangements_maxMissingTailNTs, opt$anchorReadRearrangements_minAlignmentLength)$rearrangement
      }
    } 
  }
  
  altStructSchema
}


worker <- function(x){
  library(dplyr)
  library(readr)
  library(ShortRead)
  library(data.table)
  
  db <- tmpFile()
  
  file.copy(file.path(opt$softwareDir, 'data', 'vectors', x$vectorFastaFile[1]),  
            file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(db, '.db.fasta')))
  system(paste0('makeblastdb -in ',  file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(db, '.db.fasta')),
                ' -dbtype nucl -out ', file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', db)), ignore.stderr = FALSE)
  
  o <- rbindlist(lapply(windowsTypes, function(windowType){
    rbindlist(lapply(rarifactionLevels, function(rarifactionLevel){
      rbindlist(lapply(windows, function(w){
        rbindlist(lapply(seeds, function(seed){
          
          # No need to sample additional seeds if rarifactionLevel == 0.    
          if(rarifactionLevel == 0 & which(seeds == seed) > 1) return(tibble())
          
          f <- tmpFile()
          
          # Read in anchor reads for this sample.
          d <- DNAStringSet(x$seq)
          names(d) <- x$readID
          
          d1 <- d[width(d) >= w]
          if(length(d1) > 0) d1 <- subseq(d1, 1, w)
          
          d2 <- d[width(d) < w]
          
          # Truncate sequence to this window width.
          if(windowType == 'soft'){
            dw <- c(d1, d2)
          } else {
            dw <- d1
          }
          
          if(length(dw) == 0){
            message(x$uniqueSample[1], ' - no reads meet or exceed the window (seed: ', seed, ', rarifactionLevel: ', rarifactionLevel, ', window: ', w)
            invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), pattern = f, full.names = TRUE)))
            return(data.table())
          }
          
          if(rarifactionLevel != 0){
            if(length(dw) <= rarifactionLevel){
              message(x$uniqueSample[1], ' - not enough reads in window (seed: ', seed, ', rarifactionLevel: ', rarifactionLevel, ', window: ', w)
              invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), pattern = f, full.names = TRUE)))
              return(data.table())
            }
            
            dw <- {set.seed(seed); sample(dw, rarifactionLevel)}
          }
          
          # Read in the expected vector sequence for this sample.
          # The vector may have more than one expected sequence when reading inward.
          d2 <- DNAStringSet(opt[['anchorReadRearrangements_expectedSeqs']][[unique(x$vectorFastaFile)]])
          
          # Force expected sequences to be cd-hit rep sequences by making them a little longer.
          d2 <- subseq(d2, 1, w + 5)
          names(d2) <- paste0('expectedSeq', 1:length(d2))
          
          # Write out inward read sequences and expected vector sequences to the same file.
          writeXStringSet(c(dw, d2), file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.fasta')))
          
          
          system(paste0("cd-hit-est -i ", file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.fasta')), 
                        " -o ", file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', f), 
                        " -T 2 ", opt$anchorReadRearrangements_clusterParams), ignore.stdout = TRUE, ignore.stderr = TRUE)
          
          if(! file.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.clstr')))){
            message(x$uniqueSample[1], ' - cd-hit-est did not return a result. (seed: ', seed, ', rarifactionLevel: ', rarifactionLevel, ', window: ', w)
            invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), pattern = f, full.names = TRUE)))
            return(data.table())
          }
          
          if(file.info(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.clstr')))$size == 0){
            message(x$uniqueSample[1], ' - cd-hit-est returned an empty file. (seed: ', seed, ', rarifactionLevel: ', rarifactionLevel, ', window: ', w)
            invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), pattern = f, full.names = TRUE)))
            return(data.table())
          }
          
          # Parse cluster output.
          r <- paste0(readLines(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.clstr'))), collapse = '')
          
          # Create a table where each read is a row and has a representative read id (rep)
          # For sequences that match the expected sequence, this will be 'expectedSeq' (one expected seq for each path).
          
          o <- bind_rows(lapply(unlist(strsplit(r, '>Cluster')), function(x){
            e <- sub('>', '', unlist(stringr::str_extract_all(x, '>[^\\.]+')))
            
            if(length(e) > 0){
              repSeq <- sub('^>', '', stringr::str_extract(stringr::str_extract(x, '>[^\\.]+\\.+\\s+\\*'), '>[^\\.]+'))
              return(tibble(readID = e, rep = repSeq))
            } else {
              return(tibble())
            }
          }))
          
          # Identify reads to remove.
          readsToRemove <- dplyr::group_by(o, rep) %>% 
            dplyr::mutate(nReads = n_distinct(readID)) %>% 
            dplyr::ungroup() %>% 
            dplyr::filter(nReads < opt$anchorReadRearrangements_minReadsPerAltStruct) %>%
            dplyr::pull(readID)
          
          if(length(readsToRemove) > 0 ) o <- o[! o$readID %in% readsToRemove,]
          
          if(nrow(o) == 0){
            message(x$uniqueSample[1], ' - no reads remain after removing reads with too few reads to support alt structures (seed: ', seed, ', rarifactionLevel: ', rarifactionLevel, ', window: ', w)
            invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), pattern = f, full.names = TRUE)))
            return(data.table())
          }
          
          invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), pattern = f, full.names = TRUE)))
          
          # Remove dummy expected sequence(s).
          o <- o[! grepl('expectedSeq', o$readID),]
          
          if(nrow(o) > 0){
            expectedSeqReadsIDs   <- unique(o[grepl('expectedSeq', o$rep),]$readID)
            unexpectedSeqReadsIDs <- unique(o[! grepl('expectedSeq', o$rep),]$readID)
            
            expectedStructReadIDs   <- unique(o[grepl('expectedSeq', o$rep),]$readID) 
            unExpectedStructReadIDs <- unique(o[! grepl('expectedSeq', o$rep),]$readID)
            
            altStructNum <- n_distinct(o[! grepl('expectedSeq', o$rep),]$rep)
            
            expectedStructReadNum   <- ifelse(length(expectedStructReadIDs) == 0, 0, n_distinct(expectedStructReadIDs))
            unExpectedStructReadNum <- ifelse(length(unExpectedStructReadIDs) == 0, 0, n_distinct(unExpectedStructReadIDs))
          } else {
            expectedSeqReadsIDs   <- NA
            unexpectedSeqReadsIDs <- NA
            expectedStructReadIDs <- NA
            unExpectedStructReadIDs <- NA
            altStructNum <- NA
            
            expectedStructReadNum <- NA
            unExpectedStructReadNum <- NA
          }
          
          # Write our representative structures.
          rtab <- data.frame(table(o$rep))
          
          rtab <- rtab[! grepl('expected', rtab$Var1),]
          
          if(nrow(rtab) > 0){
            rtab$uniqueSample = x$uniqueSample[1]
            rtab$sample = x$sample[1]
            rtab$windowType = windowType
            rtab$rarifactionLevel = rarifactionLevel
            rtab$window = w
            rtab$totalReads = expectedStructReadNum + unExpectedStructReadNum
            rtab$percent = round((rtab$Freq / rtab$totalReads)*100, 2)
            rtab$seed = seed
            rtab$vector = x$vectorFastaFile[1]
            rtab <- left_join(rtab, select(x, readID, seq), by = c('Var1' = 'readID'))
            rtab$seq <- substr(rtab$seq, 1, w)
            rtab$seqModel <- sapply(rtab$seq, buildSeqModel, 
                                    db = file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', db),
                                    id = f,
                                    tmpDir = file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'))
            
            write_tsv(rtab, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences', paste0(x$uniqueSample[1], '_', windowType, '_', rarifactionLevel, '_',  w, '_', seed)))
          }
          
          
          
          data.table(sampleReplicate          = x$uniqueSample[1], 
                     sample                   = x$sample[1],
                     seed                     = seed,
                     rarifactionLevel         = rarifactionLevel,
                     vector                   = x$vectorFastaFile[1],
                     vectorPaths              = length(opt$anchorReadRearrangements_expectedSeqs[[x$vectorFastaFile[1]]]),
                     window                   = w,
                     windowType               = windowType,
                     totalReads               = expectedStructReadNum + unExpectedStructReadNum, 
                     expectedStructReads      = expectedStructReadNum,
                     unExpectedStructReads    = unExpectedStructReadNum,
                     altStructs               = altStructNum,
                     percentAltStructReads    = (unExpectedStructReadNum / (expectedStructReadNum + unExpectedStructReadNum))*100)
        }))
      }))
    }))
  }))
  
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), pattern = db, full.names = TRUE)))
  o
}

updateLog('Requesting CPU cluster and exporting objects.')
cluster <- makeCluster(opt$anchorReadRearrangements_CPUs)
clusterExport(cluster, c('opt', 'waitForFile', 'tmpFile', 'windowsTypes', 'rarifactionLevels', 'windows', 'seeds', 'buildSeqModel', 'buildRearrangementModel', 'blast2rearangements'))

reads$sample <- sub('~\\d+$', '', reads$uniqueSample)

if(opt$anchorReadRearrangements_calcReplicateLevelStats){
  updateLog('Calculating rearrangements for sample replicates.')
  updateMasterLog()
  
  r.replicates <- rbindlist(parLapply(cluster, split(reads, paste(reads$vectorFastaFile, reads$uniqueSample)), worker))
  
  r.replicates$altStructsPerKB <- (r.replicates$altStructs / (r.replicates$window * r.replicates$vectorPaths))*1000
  r.replicates$sample <- NULL
  
  f <- list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences'), full.names = TRUE)
  r.replicates.altStructs <- bind_rows(lapply(f, readr::read_tsv, show_col_types = FALSE))
  r.replicates.altStructs$sample <- NULL
  
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), full.names = TRUE)))
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences'), full.names = TRUE)))
  
  updateLog('Writing outputs.')
  saveRDS(r.replicates, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'result_replicates.rds'))
  r.replicates.altStructs <- dplyr::rename(r.replicates.altStructs, repReadID = Var1, count = Freq)
  saveRDS(r.replicates.altStructs, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'result_replicates_altStructs.rds'))
}

if(opt$anchorReadRearrangements_calcSampleLevelStats){
  updateLog('Calculating rearrangements for samples.')
  updateMasterLog()
  
  r.samples <- rbindlist(parLapply(cluster, split(reads, paste(reads$vectorFastaFile, reads$sample)), worker))
  r.samples$altStructsPerKB <- (r.samples$altStructs / (r.samples$window * r.samples$vectorPaths))*1000
  r.samples$sampleReplicate <- NULL
  
  f <- list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences'), full.names = TRUE)
  r.samples.altStructs <- bind_rows(lapply(f, readr::read_tsv, show_col_types = FALSE))
  r.samples.altStructs$uniqueSample <- NULL
  
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), full.names = TRUE)))
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences'), full.names = TRUE)))
  
  updateLog('Writing outputs.')
  saveRDS(r.samples, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'result_samples.rds'))
  r.samples.altStructs <- dplyr::rename(r.samples.altStructs, repReadID = Var1, count = Freq)
  saveRDS(r.samples.altStructs, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'result_samples_altStructs.rds'))
}

updateLog('anchorReadRearrangements completed.')
updateMasterLog()

q(save = 'no', status = 0, runLast = FALSE) 