#!/usr/bin/Rscript

# AAVengeR/anchorReadRearrangements.R
# John K. Everett, Ph.D.
# 
# This script accepts the a demultiplexed read table from the demultiplex module
# and identifies instances of vector rearrangements. This analysis focus on anchor
# reads and it is often applied to sequencing experiments where anchor reads lengths
# are prioritized over adrift read lengths.

suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) stop('Expected at least one command line argument')
source(file.path(yaml::read_yaml(args[1])$softwareDir, 'lib.R'))
opt <- startModule(args)

createOuputDir()
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir))
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs'))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs'))
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'))
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences'))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences'))
invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), full.names = TRUE)))

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
clusterExport(cluster, c('opt', 'quitOnErorr', 'buildRearrangementModel'))


# Read in the reads table.
updateLog('Reading in demultiplexed reads.')
if(! file.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_inputFile))) quitOnErorr('Error - the input data file does not exist.')
reads <- readRDS(file.path(opt$outputDir, opt$anchorReadRearrangements_inputFile))
if(nrow(reads) == 0) quitOnErorr('Error - the input table contained no rows.')


# (!) Hot fix.
# (!) -------------
reads[grepl('pTBG_HC_VCN20', reads$uniqueSample)]$vectorFastaFile <- 'Sabatino_PMC7855056_heavyChain_plasmid.fasta'
# (!) -------------


# Test for raw reads, vector files, and 500mers in config file.
if(! all(unique(reads$vectorFastaFile) %in% names(opt$anchorReadRearrangements_expectedSeqs))) quitOnErorr('Error: one or more vector expected sequences are missing from the configuration file.')
if(! all(sapply(unique(reads$vectorFastaFile), function(x) file.exists(file.path(opt$softwareDir, 'data', 'vectors', x))))) quitOnErorr('Error: one or more vector sequences are missing from the AAVengeR data/vectors folder.')

# Read in the raw read data since the demultiplex module has already quality trimmed data. 
updateLog('Reading in raw reads.')
rawAnchorReads    <- readFastq(file.path(opt$demultiplex_anchorReadsFile))
rawAdriftReads    <- readFastq(file.path(opt$demultiplex_adriftReadsFile))
rawAnchorReads@id <- BStringSet(sub('\\s+.+$', '', as.character(rawAnchorReads@id)))
rawAdriftReads@id <- BStringSet(sub('\\s+.+$', '', as.character(rawAdriftReads@id)))


updateLog('Building overlapped reads.')
reads <- bind_rows(lapply(split(reads, reads$uniqueSample), function(x){
  f <- paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')
  
  updateLog(paste0('Processing ', x$uniqueSample[1], '.'))
  
  # Retrieve raw read sequences.
  R1 <- rawAdriftReads[rawAdriftReads@id %in% x$readID]
  R2 <- rawAnchorReads[rawAnchorReads@id %in% x$readID]
  
  # Quality trim read tails.
  updateLog('Quality trimming tails.')
  R1 <- trimTailw(R1, opt$anchorReadRearrangements_qualtrim_events, opt$anchorReadRearrangements_qualtrim_code, opt$anchorReadRearrangements_qualtrim_halfWidth)
  R2 <- trimTailw(R2, opt$anchorReadRearrangements_qualtrim_events, opt$anchorReadRearrangements_qualtrim_code, opt$anchorReadRearrangements_qualtrim_halfWidth)
  
  # Limit R1 reads to those long enough to trim linker sequences.
  R1 <- R1[width(R1) > max(x$adriftLinkerSeqEnd)+1]
  if(length(R1) == 0) return(tibble())
  
  # Trim R1 linker sequences using position information provided by demultiplex module.
  R1 <- narrow(R1, (x[match(R1@id, x$readID),]$adriftLinkerSeqEnd + 1), width(R1))
  
  # Trim R2 over-reading into the common linker sequence on R1.
  updateLog('Using cutAdapt to trim anchor read over-reading.')
  writeFastq(R2, paste0(f, '.R2.fastq'), compress = FALSE)
  system(paste0('cutadapt -e ', opt$anchorReadRearrangements_cutAdaptErrorRate, ' -a ', x$adriftReadTrimSeq[1], ' ', paste0(f, '.R2.fastq'), ' > ', paste0(f, '.R2.cutAdapt')), ignore.stderr = TRUE)
  R2 <- readFastq(paste0(f, '.R2.cutAdapt'))
  invisible(file.remove(list.files(pattern = f)))
  
  # Select reads that are at least 50 nt to move forward.
  R1 <- R1[width(R1) >= 50]
  R2 <- R2[width(R2) >= 50]
  
  # Limit reads to those with read IDs in both R1 and R2.
  i <- base::intersect(as.character(R1@id), as.character(R2@id))
  updateLog('Limiting reads to those with both R1 and R2 reads.')
  if(length(i) == 0) return(tibble())
  
  R1 <- R1[R1@id %in% i]
  R2 <- R2[R2@id %in% i]
  if(! all(R1@id == R2@id)) stop('Read id order error.')
  
  updateLog(paste0(length(i), ' read pairs progressing to Pear.'))
  
  # Merge R1 and R2 where possible.
  writeFastq(R1, paste0(f, '.R1.fastq'), compress = FALSE)
  writeFastq(R2, paste0(f, '.R2.fastq'), compress = FALSE)
  
  updateLog('Starting Pear.')
  system(paste0('pear ',
                ' -j ', opt$anchorReadRearrangements_CPUs, 
                ' -p ', opt$anchorReadRearrangements_maxReadOverlapPval,
                ' -v ', opt$anchorReadRearrangements_minReadOverlapNTs,
                ' -f ', paste0(f, '.R1.fastq'), 
                ' -r ', paste0(f, '.R2.fastq'), 
                ' -o ', f))
  
  updateLog('Pear completed.')
  a <- readFastq(paste0(f, '.assembled.fastq'))
  b <- readFastq(paste0(f, '.unassembled.reverse.fastq'))
  
  # Include only reads that overlap if anchorReadRearrangements_requireReadPairOverlap = TRUE.
  if(! opt$anchorReadRearrangements_requireReadPairOverlap){
    o <- Reduce('append', list(a, b))
  } else {
    o <- a
  }
  
  updateLog(paste0(length(o), ' reads acquired.'))
  
  # Cleanup tmp files.
  invisible(file.remove(list.files(pattern = f)))
  
  r <- tibble()
  
  if(length(o) > 0){
    r <- tibble(readID = as.character(o@id),
                uniqueSample = x$uniqueSample[1],
                vectorFastaFile = x$vectorFastaFile[1],
                seq = as.character(reverseComplement(o@sread)))
    
    r$mergedRead <- r$readID %in% a@id
  }
  
  r
}))


reads$w <- nchar(reads$seq)

# Limit reads to those that start with the expected sequences.
updateLog('Limiting reads to select anchor read start sequences.')

reads <- rbindlist(lapply(split(reads, reads$vectorFastaFile), function(x){
  v <- unique(x$vectorFastaFile)
  updateLog(paste0('Limiting ', ppNum(nrow(x)), ' reads for vector:  ', v))
  
  testSeqs <- unique(substr(opt[['anchorReadRearrangements_expectedSeqs']][[v]], 1, opt$anchorReadRearrangements_vectorLeaderSeqFilterLength))
  
  rbindlist(lapply(testSeqs, function(k){
    updateLog(paste0('  Limiting reads for this vector to those starting with ', k))
    r <- x[vcountPattern(k,
                         DNAStringSet(substr(x$seq, 1, nchar(k))),
                         max.mismatch = opt$anchorReadRearrangements_vectorLeaderSeqFilterMaxMismatch) == 1,]
    updateLog(paste0('    ', ppNum(nrow(r)), ' reads matched the leader sequence.'))
    r
  }))
}))

if(nrow(reads) == 0) quitOnErorr('Error - not reads matched any of the provided leader sequences.')

reads$sample <- sub('~\\d+$', '', reads$uniqueSample)

nReadsPreFilter <- nrow(reads)

# Remove reads whoes ends do not align with their vector plasmid.
reads <- bind_rows(lapply(split(reads, reads$vectorFastaFile), function(x){
  f <- tmpFile()
  x$lastNTs <- substr(x$seq, nchar(x$seq) - (opt$anchorReadRearrangements_readEndAlignmentTestLength - 1), nchar(x$seq))
  
  system(paste0('faToTwoBit ',
                file.path(opt$softwareDir, 'data', 'vectors', x$vectorFastaFile[1]), ' ',
                file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', f)))
  
  d <- DNAStringSet(unique(x$lastNTs))
  names(d) <- as.character(d)
  writeXStringSet(d,  file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.fasta')))
  
  system(paste0('blat ', file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', f), ' ',
                file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.fasta')), ' ',
                file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.psl')),
                ' -tileSize=6 -stepSize=3 -repMatch=3000 -out=psl -t=dna -q=dna -minScore=0 -minIdentity=0 -noHead -noTrimA'))
  
  o <- parseBLAToutput(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp', paste0(f, '.psl')))
  o <- subset(o, matches >= opt$anchorReadRearrangements_readEndAlignmentTestMinMatch)
  
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'), pattern = f, full.names = TRUE)))
  
  # add alignment failure test...
  subset(x, lastNTs %in% o$qName)
}))

updateLog(paste0(sprintf("%.2f%%", (n_distinct(reads$readID) / nReadsPreFilter)*100), ' reads remain after removing reads whoes tails do not align to the vector.'))

samplesProcessed <- 1

r <- bind_rows(lapply(split(reads, paste(reads$vectorFastaFile, reads$uniqueSample)), function(x){
  updateLog(paste0('Processing sample ', samplesProcessed, ' / ', n_distinct(reads$uniqueSample), ' - replicate name: ', x$uniqueSample[1]))
  samplesProcessed <<- samplesProcessed  + 1
  
  # Read in anchor reads for this sample.
  d <- DNAStringSet(x$seq)
  names(d) <- x$readID
  
  o <- bind_rows(lapply(c(100, 200, 300, 400), function(w){
    
    # Truncate sequence to this window width.
    dw1 <- d[width(d) >= w]
    dw1 <- subseq(dw1, 1, w)
    
    dw2 <- d[width(d) <= w]
    
    if(opt$anchorReadRearrangements_useHardWindows){
      dw <- dw1
    } else {
      dw <- c(dw1, dw2)
    }
    
    # Read in the expected vector sequence for this sample.
    # The vector may have more than one expected sequence when reading inward.
    d2 <- DNAStringSet(opt[['anchorReadRearrangements_expectedSeqs']][[unique(x$vectorFastaFile)]])
    
    # Force expected sequences to be cd-hit rep sequences by making them a little longer.
    d2 <- subseq(d2, 1, w + 5)
    names(d2) <- paste0('expectedSeq', 1:length(d2))
    
    # Write out inward read sequences and expected vector sequences to the same file.
    writeXStringSet(c(dw, d2), 'tmpDNAseq.fasta')
    
    # Cluster sequences with cd-hit-est.
    system(paste0("cd-hit-est -i tmpDNAseq.fasta -o tmpDNAseq -T ", 
                  opt$anchorReadRearrangements_CPUs, ' ',
                  opt$anchorReadRearrangements_clusterParams))
    
    # Parse cluster output.
    r <- paste0(readLines('tmpDNAseq.clstr'), collapse = '')
    invisible(file.remove(list.files(pattern = 'tmpDNAseq', full.names = TRUE)))
    
    o <- bind_rows(lapply(unlist(strsplit(r, '>Cluster')), function(x){
      e <- sub('>', '', unlist(stringr::str_extract_all(x, '>[^\\.]+')))
      
      if(length(e) > 0){
        repSeq <- sub('^>', '', stringr::str_extract(stringr::str_extract(x, '>[^\\.]+\\.+\\s+\\*'), '>[^\\.]+'))
        return(tibble(readID = e, rep = repSeq))
      } else {
        return(tibble())
      }
    }))
    
    # Identify reads.
    readsToRemove <- dplyr::group_by(o, rep) %>% 
      dplyr::mutate(nReads = n_distinct(readID)) %>% 
      dplyr::ungroup() %>% 
      dplyr::filter(nReads < opt$anchorReadRearrangements_minReadsPerAltStruct) %>%
      dplyr::pull(readID)
    
    if(length(readsToRemove) > 0 ){
      updateLog(paste0('Removing ', length(readsToRemove), ' reads since they support alt structs with less than ', opt$anchorReadRearrangements_minReadsPerAltStruct, ' reads.'))
      o <- o[! o$readID %in% readsToRemove,]
    }
    
    # Write out clustered sequences named by the cluster rep sequence read.
    invisible(lapply(split(o, o$rep), function(clust){
      k <- subset(x, readID %in% clust$readID)
      d <- DNAStringSet(k$seq)
      names(d) <- k$readID
      if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences', x$uniqueSample[1]))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences', x$uniqueSample[1]))
      writeXStringSet(d, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences', x$uniqueSample[1], paste0('Clust_', clust$rep[1], '.fasta')))
    }))
    
    # Remove dummy expected sequence(s).
    o <- o[! grepl('expectedSeq', o$readID),]
    
    expectedSeqReadsIDs <- unique(o[grepl('expectedSeq', o$rep),]$readID)
    expectedSeqReadsNTs <- sum(nchar(subset(reads, readID %in% expectedSeqReadsIDs)$seq))
    
    unexpectedSeqReadsIDs <- unique(o[! grepl('expectedSeq', o$rep),]$readID)
    unexpectedSeqReadsNTs <- sum(nchar(subset(reads, readID %in% unexpectedSeqReadsIDs)$seq))
    
    allReadNTs <- sum(nchar(x$seq))
    allReadIDs <- unique(o$readID)
    
    altStructBreaks <- NA 
    
    # Remove all reads associated with expected sequences.
    o <- o[! grepl('expectedSeq', o$rep),]
    
    if(nrow(o) > 0){
      # Make a blast database for the vector sequence.
      files <- list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs'), full.names = TRUE)
      if(length(files) > 0) invisible(file.remove(files))
      
      system(paste0('/home/ubuntu/software/ncbi-blast-2.12.0+/bin/makeblastdb -in ', file.path(opt$softwareDir, 'data', 'vectors', x$vector[1]), 
                    ' -dbtype nucl -out ', file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs', 'd')), ignore.stderr = FALSE)
      
      # Here we cluster representatives to control for wiggle between recombined forms.
      repReads <- DNAStringSet(x[x$readID %in% o$rep,]$seq)
      names(repReads) <- x[x$readID %in% o$rep,]$readID
      
      writeXStringSet(repReads, 'repReads.fasta')
      
      # Align the chunk to the vector sequence.  repReads.fasta
      #
      # (!) Here we increase the penalty from the default -3 to -4 to prevent runs of mismatches to come through if followed by a string of matches.
      #
      system(paste0('/home/ubuntu/software/ncbi-blast-2.12.0+/bin/blastn -penalty -4 -max_target_seqs 10000 -gapopen 10 -gapextend 5 -dust no -soft_masking false -word_size 5 -evalue 100 -outfmt 6 -query repReads.fasta -db ',
                    file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs', 'd'),
                    ' -num_threads ', opt$anchorReadRearrangements_CPUs  ,' -out repReads.blast'), ignore.stdout = TRUE, ignore.stderr = TRUE)
      
      invisible(file.remove('repReads.fasta'))
      
      if(file.info('repReads.blast')$size != 0){
        # Parse blastn result.
        b <- read.table('repReads.blast', sep = '\t', header = FALSE)
        invisible(file.remove('repReads.blast'))
        names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
        b$alignmentLength <- b$qend - b$qstart + 1
        b$strand <- ifelse(b$sstart > b$send, '-', '+')
        
        b <- dplyr::filter(b, pident >= opt$anchorReadRearrangements_minAlignmentPercentID, alignmentLength >= opt$anchorReadRearrangements_minAlignmentLength)
        
        # Create a grouping index for read IDs.
        b$i <- group_by(b, qname) %>% group_indices()
        
        # Create a CPU splitting vector that would not separate the read ID indices. 
        k <- tibble(i2 = 1:n_distinct(b$i))
        k$n <- ntile(1:nrow(k), opt$anchorReadRearrangements_CPUs)
        
        # Bind the CPU splitting vector.
        b <- left_join(b, k, by = c('i' = 'i2'))
        
        # qlen needed by blast2rearangements.
        b <- left_join(b, tibble(qname = names(dw), qlen = width(dw)), by = 'qname')
        
        r <- rbindlist(parallel::parLapply(cluster, split(b, b$n), blast2rearangements, opt$anchorReadRearrangements_maxMissingTailNTs, opt$anchorReadRearrangements_minAlignmentLength))
        
        r <- left_join(r, data.frame(table(o$rep)), by = c('qname' = 'Var1'))
        names(r) <- c('clusterRepReadID', 'rearrangement', 'readsInCluster')
        r <- arrange(r, desc(readsInCluster))
        
        r <- mutate(r, sample = x$uniqueSample[1], .before = 'clusterRepReadID')
        readr::write_tsv(r, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'repLeaderSeqStructs.tsv'), append = TRUE)
        
        altStructBreaks <- sum(stringr::str_count(r$rearrangement, ';'))
      }
    }   
    
    tibble(sample                   = x$uniqueSample[1], 
           window                   = w,
           totalReads               = n_distinct(allReadIDs), 
           totalNTs                 = allReadNTs,
           expectedStructReads      = n_distinct(expectedSeqReadsIDs),
           expectedStructNTs        = expectedSeqReadsNTs,
           altStructs               = n_distinct(o$rep), 
           altStructReads           = n_distinct(unexpectedSeqReadsIDs),
           altStructNTs             = unexpectedSeqReadsNTs,
           altStructBreaks          = altStructBreaks,
           percentAltStructReads    = (n_distinct(unexpectedSeqReadsIDs)/n_distinct(allReadIDs))*100)
  }))
}))

stopCluster(cluster)

saveRDS(r, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'result.rds'))
openxlsx::write.xlsx(r, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'result.xlsx'))

o <- mutate(r, sample2 = sub('~\\d+$', '', sample)) %>%
     group_by(sample2, window) %>%
     summarise(avgRepReads = mean(totalReads),
              avgAltStructs = mean(altStructs),
              avgPercentAltStructReads = mean(percentAltStructReads),
              sd = sd(percentAltStructReads)) %>%
    ungroup() %>%
    mutate(altStructBin = cut(avgAltStructs,
                              breaks = c(0, 10, 50, 100, 250, 500, Inf), labels = FALSE, right = FALSE))

ggplot(o, aes(window, avgPercentAltStructReads, group = sample2, color= sample2, fill = sample2)) +
theme_bw() +
geom_line(position = position_dodge(width = 30), guide = 'none') +
geom_jitter(aes(size = altStructBin),
            shape = 21,
            color = 'black',
            position = position_dodge(width = 30)) +
scale_fill_manual(name = 'Samples',
                  values = grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_distinct(o$sample2))) +
scale_color_manual(values = grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_distinct(o$sample2))) +
scale_size(name = 'Alternative structures',
           range = c(2, 6),
           breaks = c(1:6),
           labels = c('0 - 10', '10 - 50', '50 - 100', '100 - 250', '250 - 500', '> 500'),
           guide = "legend") +
scale_y_continuous(limits = c(0, 30)) +
geom_errorbar(aes(ymin = avgPercentAltStructReads - sd,
                  ymax = avgPercentAltStructReads + sd),
                  position = position_dodge(width = 30), width = 30) +
labs(x = 'Read width', y = 'Percent rearranged') +
theme(text = element_text(size = 14)) +
guides(color = "none", fill = guide_legend(override.aes = list(size = 6)))
