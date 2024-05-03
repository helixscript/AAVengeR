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

# Read the configuration file and set additional parameters.
configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')

opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

createOuputDir()
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir))
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs'))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs'))
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'tmp'))
if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences'))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences'))
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
clusterExport(cluster, c('opt', 'quitOnErorr', 'buildRearrangementModel'))

if(! file.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_inputFile))) quitOnErorr('Error - the input data file does not exist.')

# Read in the reads table.
updateLog('Reading in demultiplexed reads.')

reads <- readRDS(file.path(opt$outputDir, opt$anchorReadRearrangements_inputFile))
if(nrow(reads) == 0) quitOnErorr('Error - the input table contained no rows.')

# Update
reads[grepl('_LC~', reads$uniqueSample),]$vectorFastaFile <- 'Sabatino_PMC7855056_lightChain_plasmid.fasta'
reads[grepl('_HC~', reads$uniqueSample),]$vectorFastaFile <- 'Sabatino_PMC7855056_heavyChain_plasmid.fasta'

rawAnchorReads <- readFastq(file.path(opt$demultiplex_anchorReadsFile))
rawAdriftReads <- readFastq(file.path(opt$demultiplex_adriftReadsFile))

rawAnchorReads@id <- BStringSet(sub('\\s+.+$', '', as.character(rawAnchorReads@id)))
rawAdriftReads@id <- BStringSet(sub('\\s+.+$', '', as.character(rawAdriftReads@id)))


z <- rawAnchorReads[rawAnchorReads@id %in% subset(reads, vectorFastaFile == 'pAAV-GFP-plasmid.fasta')$readID]
i1 <- vcountPattern('TTCCTTGTAGTTAATGATTAACCCGCCATGCTACTTATCTACGTAGCCATGCTCTAGCGGCCTCGGCC', z@sread, max.mismatch = 1) == 1
i2 <- vcountPattern('TTCCTTGTAGTTAATGATTAACCCGCCATGCTACTTATCTACGTAGCCATGCTCTAGATGTCGAGGGA', z@sread, max.mismatch = 1) == 1



# r <- readFasta('readsToRemove')
# reads <- reads[! reads$readID %in% as.character(r@id),]
# rawAnchorReads <- rawAnchorReads[! rawAnchorReads@id %in% as.character(r@id)]
# rawAdriftReads <- rawAdriftReads[! rawAdriftReads@id %in% as.character(r@id)]


reads <- bind_rows(lapply(split(reads, reads$uniqueSample), function(x){
           f <- paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')
  
           # Retrieve raw read sequences.
           R1 <- rawAdriftReads[rawAdriftReads@id %in% x$readID]
           R2 <- rawAnchorReads[rawAnchorReads@id %in% x$readID]
           
           # Quality trim read tails.
           R1 <- trimTailw(R1, opt$anchorReadRearrangements_qualtrim_events, opt$anchorReadRearrangements_qualtrim_code, opt$anchorReadRearrangements_qualtrim_halfWidth)
           R2 <- trimTailw(R2, opt$anchorReadRearrangements_qualtrim_events, opt$anchorReadRearrangements_qualtrim_code, opt$anchorReadRearrangements_qualtrim_halfWidth)
           
           # Limit R1 reads to those long enough to trim linker sequences.
           R1 <- R1[width(R1) > max(x$adriftLinkerSeqEnd)+1]
           if(length(R1) == 0) return(tibble())
 
           # Trim R1 linker sequences using position information provided by demultiplex module.
           R1 <- narrow(R1, (x[match(R1@id, x$readID),]$adriftLinkerSeqEnd + 1), width(R1))

           # Trim R2 over-reading into the common linker sequence on R1.
           writeFastq(R2, paste0(f, '.R2.fastq'), compress = FALSE)
           system(paste0('cutadapt -e ', opt$anchorReadRearrangements_cutAdaptErrorRate, ' -a ', x$adriftReadTrimSeq[1], ' ', paste0(f, '.R2.fastq'), ' > ', paste0(f, '.R2.cutAdapt')), ignore.stderr = TRUE)
           R2 <- readFastq(paste0(f, '.R2.cutAdapt'))
           invisible(file.remove(list.files(pattern = f)))
           
           # Select reads that are at least 50 nt to move forward.
           R1 <- R1[width(R1) >= 50]
           R2 <- R2[width(R2) >= 50]
           
           # Limit reads to those with read IDs in both R1 and R2.
           i <- base::intersect(as.character(R1@id), as.character(R2@id))
           if(length(i) == 0) return(tibble())
           
           R1 <- R1[R1@id %in% i]
           R2 <- R2[R2@id %in% i]
           if(! all(R1@id == R2@id)) stop('Read id order error.')
           
           # Merge R1 and R2 where possible.
           writeFastq(R1, paste0(f, '.R1.fastq'), compress = FALSE)
           writeFastq(R2, paste0(f, '.R2.fastq'), compress = FALSE)
           
           message(x$uniqueSample[1])
           system(paste0('pear ',
                         ' -j ', opt$anchorReadRearrangements_CPUs, 
                         ' -p ', opt$anchorReadRearrangements_maxReadOverlapPval,
                         ' -v ', opt$anchorReadRearrangements_minReadOverlapNTs,
                         ' -f ', paste0(f, '.R1.fastq'), 
                         ' -r ', paste0(f, '.R2.fastq'), 
                         ' -o ', f))
  
           a <- readFastq(paste0(f, '.assembled.fastq'))
           b <- readFastq(paste0(f, '.unassembled.reverse.fastq'))
  
           # Include only reads that overlap if anchorReadRearrangements_requireReadPairOverlap = TRUE.
           if(! opt$anchorReadRearrangements_requireReadPairOverlap){
             o <- Reduce('append', list(a, b))
           } else {
             o <- a
           }
           
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
reads$sample <- sub('~\\d+$', '', reads$uniqueSample)
p <- ggplot(reads, aes(w)) + 
     theme_bw() +
     geom_histogram(bins = 80) + 
     geom_vline(xintercept = 200, color = 'red') +
      facet_wrap(sample~., scales = 'free_y', ncol = 3)
ggsave(p, units = 'in', height = 20, width = 20, file = 'p.pdf')
saveRDS(reads, 'myReads.rds')

reads <- reads[reads$w >= opt$anchorReadRearrangements_anchorReadCutPos,]
reads$seq <- substr(reads$seq, 1, opt$anchorReadRearrangements_anchorReadCutPos)

nReadsPreFilter <- nrow(reads)

if(! all(unique(reads$vectorFastaFile) %in%  names(opt$anchorReadRearrangements_expectedSeqs))){
  updateLog('Error -- all the vector names found in the read data are not defined in the anchorReadRearrangements_expectedSeqs section of the config file.')
  updateLog(paste0('These vectors are missing: ', paste0(unique(reads$vectorFastaFile)[! unique(reads$vectorFastaFile) %in%  names(opt$anchorReadRearrangements_expectedSeqs)], collapse = ',')))
  q()
}

# Limit reads to those that start wit the expected sequences.
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

updateLog(paste0(sprintf("%.2f%%", (nrow(reads) / nReadsPreFilter)*100), ' reads retained.'))
if(nrow(reads) == 0) quitOnErorr('Error - not reads matched any of the provided leader sequences.')


reads$sample <- sub('~\\d+$', '', reads$uniqueSample)

nReadsPreFilter <- nrow(reads)

# Remove reads whoes ends do not align with their vector plasmid.
reads <- bind_rows(lapply(split(reads, reads$vectorFastaFile), function(x){
           f <- tmpFile()
           x$lastNTs <- substr(x$seq, nchar(x$seq) - (opt$anchorReadRearrangements_readEndAlignmentTestLength - 1), nchar(x$seq))

           system(paste0('/home/ubuntu/software/faToTwoBit ',
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

r <- bind_rows(lapply(split(reads, paste(reads$vectorFastaFile, reads$sample)), function(x){
       updateLog(paste0('Processing sample ', samplesProcessed, ' / ', n_distinct(reads$sample), ' - sample name: ', x$sample[1]))
       samplesProcessed <<- samplesProcessed  + 1
       
       # Read in anchor reads for this sample.
       d <- DNAStringSet(x$seq)
       names(d) <- x$readID
       
       # Read in the expected vector sequence for this sample.
       # The vector may have more than one expected sequence when reading inward.
       d2 <- DNAStringSet(opt[['anchorReadRearrangements_expectedSeqs']][[unique(x$vectorFastaFile)]])
       d2 <- unique(subseq(d2, 1, max(nchar(x$seq)) + 5))
       names(d2) <- paste0('expectedSeq', 1:length(d2))
       message('vector: ', unique(x$vectorFastaFile), ' nSeqs: ', length(d2))
       
       # Write out inward read sequences and expected vector sequences to the same file.
       writeXStringSet(c(d, d2), 'tmpDNAseq.fasta')
       
       system(paste0("cd-hit-est -i tmpDNAseq.fasta -o tmpDNAseq -T ", 
                     opt$anchorReadRearrangements_CPUs, ' ',
                     opt$anchorReadRearrangements_clusterParams))
       
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
       
       readsToRemove <- dplyr::group_by(o, rep) %>% 
                        dplyr::summarise(nReads = n_distinct(readID)) %>% 
                        dplyr::ungroup() %>% 
                        dplyr::filter(nReads < opt$anchorReadRearrangements_minReadsPerAltStruct) %>%
                        dplyr::pull(rep)
       
       if(length(readsToRemove) > 0 ){
         message('Removing ', length(readsToRemove), ' reads.')
         message(nrow(o))
         o <- o[! o$readID %in% readsToRemove,]
         message(nrow(o))
       }

       invisible(lapply(split(o, o$rep), function(clust){
         k <- subset(x, readID %in% clust$readID)
         d <- DNAStringSet(k$seq)
         names(d) <- k$readID
         if(! dir.exists(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences', x$sample[1]))) dir.create(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences', x$sample[1]))
         writeXStringSet(d, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'clusterSequences', x$sample[1], paste0('Clust_', clust$rep[1], '.fasta')))
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
       
       o <- o[! grepl('expectedSeq', o$rep),]
       
       if(nrow(o) > 0){
         # Make a blast database for the vector sequence.
         files <- list.files(file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs'), full.names = TRUE)
         if(length(files) > 0) invisible(file.remove(files))
         system(paste0('makeblastdb -in ', file.path(opt$softwareDir, 'data', 'vectors', x$vector[1]), 
                       ' -dbtype nucl -out ', file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
         
         # Here we cluster representatives to control for wiggle between recombined forms.
         repReads <- DNAStringSet(x[x$readID %in% o$rep,]$seq)
         names(repReads) <- x[x$readID %in% o$rep,]$readID
         
         writeXStringSet(repReads, 'repReads.fasta')
         
         # Align the chunk to the vector sequence.  repReads.fasta
         #
         # (!) Here we increase the penalty from the default -3 to -4 to prevent runs of mismatches to come through if followed by a string of matches.
         #
         system(paste0('blastn -penalty -4 -max_target_seqs 10000 -gapopen 10 -gapextend 5 -dust no -soft_masking false -word_size 5 -evalue 100 -outfmt 6 -query repReads.fasta -db ',
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
           
           b$qlen <- opt$anchorReadRearrangements_anchorReadCutPos
           
           b <- dplyr::filter(b, pident >= opt$anchorReadRearrangements_minAlignmentPercentID, alignmentLength >= opt$anchorReadRearrangements_minAlignmentLength)
           
           # Create a grouping index for read IDs.
           b$i <- group_by(b, qname) %>% group_indices()
           
           # Create a CPU splitting vector that would not separate the read ID indices. 
           k <- tibble(i2 = 1:n_distinct(b$i))
           k$n <- ntile(1:nrow(k), opt$anchorReadRearrangements_CPUs)
           
           # Bind the CPU splitting vector.
           b <- left_join(b, k, by = c('i' = 'i2'))
           
           r <- rbindlist(parallel::parLapply(cluster, split(b, b$n), blast2rearangements, opt$anchorReadRearrangements_maxMissingTailNTs, opt$anchorReadRearrangements_minAlignmentLength))
           
           r <- left_join(r, data.frame(table(o$rep)), by = c('qname' = 'Var1'))
           names(r) <- c('clusterRepReadID', 'rearrangement', 'readsInCluster')
           r <- arrange(r, desc(readsInCluster))
           
           r <- mutate(r, sample = x$sample[1], .before = 'clusterRepReadID')
           readr::write_tsv(r, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'repLeaderSeqStructs.tsv'), append = TRUE)
           
           altStructBreaks <- sum(stringr::str_count(r$rearrangement, ';'))
         }
       }   
       
       tibble(sample                   = x$sample[1], 
              totalReads               = n_distinct(allReadIDs), 
              totalNTs                 = allReadNTs,
              expectedStructReads      = n_distinct(expectedSeqReadsIDs),
              expectedStructNTs        = expectedSeqReadsNTs,
              altStructs               = n_distinct(o$rep), 
              altStructReads           = n_distinct(unexpectedSeqReadsIDs),
              altStructNTs             = unexpectedSeqReadsNTs,
              altStructBreaks          = altStructBreaks,
              percentAltStructReads    = (n_distinct(unexpectedSeqReadsIDs)/n_distinct(allReadIDs))*100,
              altStructsPer1KB         = altStructs / (totalNTs/1000),
              altStructsPer10KB        = altStructs / (totalNTs/10000),
              altStructBreaksPer1KB    = altStructBreaks / (totalNTs/5000),    
              altStructBreaksPer10KB   = altStructBreaks / (totalNTs/10000))
}))

saveRDS(r, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'result.rds'))
openxlsx::write.xlsx(r, file.path(opt$outputDir, opt$anchorReadRearrangements_outputDir, 'result.xlsx'))

