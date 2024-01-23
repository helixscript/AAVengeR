library(lubridate)
library(dplyr)
library(parallel)
library(data.table)
library(Biostrings)
options(stringsAsFactors = FALSE)

# Read in configuration file and set parameters.
configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))
setMissingOptions()
setOptimalParameters()
set.seed(1)


# Create the required directory structure.
if(! dir.exists(file.path(opt$outputDir, opt$alignReads_outputDir))) dir.create(file.path(opt$outputDir, opt$alignReads_outputDir))
if(! dir.exists(file.path(opt$outputDir, opt$alignReads_outputDir, 'tmp'))) dir.create(file.path(opt$outputDir, opt$alignReads_outputDir, 'tmp'))
if(! dir.exists(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1'))) dir.create(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1'))
if(! dir.exists(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2'))) dir.create(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2'))


# Read in sequencing data.
reads <- readRDS(file.path(opt$outputDir, opt$alignReads_inputFile))
incomingSamples <- unique(reads$uniqueSample)


# Create a CPU cluster.
cluster <- makeCluster(opt$alignReads_CPUs)
clusterSetRNGStream(cluster, 1)
clusterExport(cluster, c('opt', 'tmpFile'))

alignReads <- function(r, refGenome, minPercentSeqID, maxQstart, dir){

  write(c(paste(lubridate::now(), 'Started running aligner')), file = file.path(opt$outputDir, opt$alignReads_outputDir, 'log'), append = TRUE)
  
  if(opt$alignReads_aligner == 'bwa2'){
        refGenomePath <-  file.path(opt$softwareDir, 'data', 'referenceGenomes', 'bwa2', refGenome)
        f <- file.path(dir, tmpFile())
        write(paste0('>', r$id, '\n', r$seq), f)
        system(paste0('bwa-mem2 mem -t ', opt$alignReads_CPUs  ,' -c 100000 -a ', refGenomePath, ' ', f, ' > ', f, '.sam'))
        system(paste0(file.path(opt$softwareDir, 'bin', 'sam2psl.py'), ' -i ', f, '.sam -o ', f, '.psl'))
        invisible(file.remove(paste0(f, '.sam')))
  }

  # If blat is requested as the aligner, run blat on sequence chunks with parLapply().
  if(opt$alignReads_aligner == 'blat'){

    refGenomePath <-  file.path(opt$softwareDir, 'data', 'referenceGenomes', 'blat', paste0(refGenome, '.2bit'))
      
    # Create an OOC file if requested
    if(opt$alignReads_blatUseOocFile){
      #system(paste0('blat ', refGenomePath, ' /dev/null /dev/null -repMatch=',
      system(paste0('/home/ubuntu/software/blat_37x1/blat ', refGenomePath, ' /dev/null /dev/null -repMatch=',
                    opt$alignReads_genomeAlignment_blatRepMatch, ' -makeOoc=',
                    file.path(opt$outputDir, opt$alignReads_outputDir, paste0(opt$alignReads_genomeAlignment_blatTileSize, '.ooc'))))
    }

    # Create a split vector that attempts to equally distributes different length sequences between CPUs.
    r$cut <- cut(nchar(r$seq), c(-Inf, seq(0, max(nchar(r$seq)), by = 10), Inf), labels = FALSE)
    r <- group_by(r, cut) %>% mutate(n = ntile(1:n(), opt$alignReads_CPUs)) %>% ungroup()
    invisible(parLapply(cluster, split(r, r$n), blat, refGenomePath, dir))
  }

  write(c(paste(lubridate::now(), 'Finished running aligner')), file = file.path(opt$outputDir, opt$alignReads_outputDir, 'log'), append = TRUE)
  
  # Read in and parse the psl files created by either bwa2 or blat.
  write(c(paste(lubridate::now(), 'Parsing BLAT PSL')), file = file.path(opt$outputDir, opt$alignReads_outputDir, 'log'), append = TRUE)

  b <- tibble()
  f <- list.files(dir, pattern = '*.psl', full.names = TRUE)
  
  if(length(f) > 0){
    b <- rbindlist(lapply(f, function(x){
    
           if(opt$alignReads_aligner == 'bwa2'){
             # bwa2 will yield a single psl file (converted from sam) to parse.
             b <- data.table(parseBLAToutput(x, convertToBlatPSL = TRUE))
           } else {
             b <- data.table(parseBLAToutput(x))
           }

           if(nrow(b) == 0) return(data.table())
         
          dplyr::filter(b, alignmentPercentID >= minPercentSeqID, tNumInsert <= 1, qNumInsert <= 1, 
                           tBaseInsert <= 2, qBaseInsert <= 2, qStart <= maxQstart) %>%
           dplyr::select(qName, strand, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, queryPercentID, tAlignmentWidth, queryWidth, alignmentPercentID, percentQueryCoverage)
        }))
  }

  write(c(paste(lubridate::now(), 'Finished running parsing BLAT PSL')), file = file.path(opt$outputDir, opt$alignReads_outputDir, 'log'), append = TRUE)

  # Files need to be removed otherwise the join below will keep joining files from previous genomes.
  invisible(file.remove(list.files(dir, full.names = TRUE)))

  b
}


# Align anchor reads.
anchorReadAlignments <- rbindlist(lapply(split(reads, reads$refGenome), function(x){
  s <- data.table(seq = unique(x$anchorReadSeq))
  s$id <- paste0('s', 1:nrow(s))
  
  b <- alignReads(s, x$refGenome[1], opt$alignReads_genomeAlignment_minPercentID, opt$alignReads_genomeAlignment_anchorRead_maxStartPos,  file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1'))
  
  if(nrow(b) > 0){
    # Use the unique sequences to map alignments back to sequences
    b <- left_join(b, distinct(dplyr::select(s, id, seq)), by = c('qName' = 'id'))
    
    # Select reads that have an alignment.
    x2 <- dplyr::select(x, uniqueSample, readID, refGenome, anchorReadSeq) %>% dplyr::filter(anchorReadSeq %in% b$seq)
    
    # Join alignments to all reads using aligned sequences as keys.
    b <- left_join(x2, b, by = c('anchorReadSeq' = 'seq'), relationship = 'many-to-many') %>% dplyr::select(-anchorReadSeq, -qName)
  }
  
  b
}))


write(c(paste0(lubridate::now(), '    ', sprintf("%.2f%%", (n_distinct(anchorReadAlignments$readID)/n_distinct(reads$readID))*100), 
               ' of prepped anchor reads aligned to the reference genome.')), file = file.path(opt$outputDir, opt$alignReads_outputDir, 'log'), append = TRUE)

# Select anchor reads where the ends align to the genome.
i <- (anchorReadAlignments$qSize - anchorReadAlignments$qEnd) <= opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned

if(sum(i) == 0){
  write(c(paste(lubridate::now(), 'Error - no anchor read alignments remain after alignReads_genomeAlignment_anchorReadEnd_maxUnaligned filter.')), file = file.path(opt$outputDir, opt$alignReads_outputDir, 'log'), append = TRUE)
  if(opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()
  q(save = 'no', status = 1, runLast = FALSE)
}

alignedReadIDsBeforeFilter <- n_distinct(anchorReadAlignments$readID)
anchorReadAlignments <- anchorReadAlignments[i,]

write(c(paste0(now(), '    ', sprintf("%.2f%%", (1 - n_distinct(anchorReadAlignments$readID) / alignedReadIDsBeforeFilter)*100), 
               ' of anchor reads removed because their alignments ended more than ',
               opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned,
               ' NTs from the end of reads.')), file = file.path(opt$outputDir, opt$alignReads_outputDir, 'log'), append = TRUE)


# Subset reads to those with good anchor read alignments.
reads <- subset(reads, readID %in% anchorReadAlignments$readID)


# Align adrift reads.
adriftReadAlignments <- rbindlist(lapply(split(reads, reads$refGenome), function(x){
  s <- data.table(seq = unique(x$adriftReadSeq))
  s$id <- paste0('s', 1:nrow(s))
  
  b <- alignReads(s, x$refGenome[1], opt$alignReads_genomeAlignment_minPercentID, opt$alignReads_genomeAlignment_adriftRead_maxStartPos,  file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2'))
  
  if(nrow(b) > 0){
    # Use the unique sequences to map alignments back to sequences
    b <- left_join(b, distinct(dplyr::select(s, id, seq)), by = c('qName' = 'id'))
    
    # Select reads that have an alignment.
    x2 <- dplyr::select(x, uniqueSample, readID, refGenome, adriftReadSeq) %>% dplyr::filter(adriftReadSeq %in% b$seq)
    
    # Join alignments to all reads using aligned sequences as keys.
    b <- left_join(x2, b, by = c('adriftReadSeq' = 'seq'), relationship = 'many-to-many') %>% dplyr::select(-adriftReadSeq, -qName)
  }
  
  b
}))


# Select adrift reads where the ends align to the genome.
i <- (adriftReadAlignments$qSize - adriftReadAlignments$qEnd) <= opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned

if(sum(i) == 0){
  write(c(paste(lubridate::now(), 'Error - no adrift read alignments remain after alignReads_genomeAlignment_adriftReadEnd_maxUnaligned filter.')), file = file.path(opt$outputDir, opt$alignReads_outputDir, 'log'), append = TRUE)
  if(opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()
  q(save = 'no', status = 1, runLast = FALSE)
}


alignedReadIDsBeforeFilter <- n_distinct(adriftReadAlignments$readID)
adriftReadAlignments <- adriftReadAlignments[i,]

write(c(paste0(now(), '    ', sprintf("%.2f%%", (1 - n_distinct(adriftReadAlignments$readID) / alignedReadIDsBeforeFilter)*100), 
               ' of anchor reads removed because their alignments ended more than ',
               opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned,
               ' NTs from the end of reads.')), file = file.path(opt$outputDir, opt$alignReads_outputDir, 'log'), append = TRUE)


# Select adrift reads with alignments that start at the beginning of reads.
i <- adriftReadAlignments$qStart <= opt$alignReads_genomeAlignment_adriftRead_maxStartPos

if(sum(i) == 0){
  write(c(paste(lubridate::now(), 'Error - no adrift read alignments remain after alignReads_genomeAlignment_adriftRead_maxStartPos filter.')), file = file.path(opt$outputDir, opt$alignReads_outputDir, 'log'), append = TRUE)
  if(opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()
  q(save = 'no', status = 1, runLast = FALSE)
}

alignedReadIDsBeforeFilter <- n_distinct(adriftReadAlignments$readID)
adriftReadAlignments <- adriftReadAlignments[i,]

write(c(paste0(now(), '    ', sprintf("%.2f%%", (1 - n_distinct(adriftReadAlignments$readID) / alignedReadIDsBeforeFilter)*100), 
               ' of anchor reads removed because their alignments ended more than ',
               opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned,
               ' NTs from the end of reads.')), file = file.path(opt$outputDir, opt$alignReads_outputDir, 'log'), append = TRUE)


# Find the intersection between retained anchor and adrift reads and limit read pairs to pairs
# where both mates aligned within parameters.

write(c(paste(now(), '   Finding read pairs where both anchor and adrift read pair mates aligned well.')), file = file.path(opt$outputDir, opt$alignReads_outputDir, 'log'), append = TRUE)


# Select read pairs that have both a good anchor and adrift alignment.
i <- base::intersect(anchorReadAlignments$readID, adriftReadAlignments$readID)
anchorReadAlignments <- anchorReadAlignments[anchorReadAlignments$readID %in% i,]
adriftReadAlignments <- adriftReadAlignments[adriftReadAlignments$readID %in% i,]


# Add refGenome, vector, and flags.
anchorReadAlignments <- left_join(anchorReadAlignments, distinct(select(reads, readID, vectorFastaFile, seqRunID, flags)), by = 'readID')


# Expand predicted leaderSeq sequences by extending with delayed alignment sequences. 
anchorReadAlignments <- left_join(anchorReadAlignments, select(reads, readID, anchorReadSeq, leaderSeq), by = 'readID')
anchorReadAlignments$leaderSeq <- paste0(anchorReadAlignments$leaderSeq, substr(anchorReadAlignments$anchorReadSeq, 1, anchorReadAlignments$qStart))
anchorReadAlignments <- dplyr::select(anchorReadAlignments, -anchorReadSeq)


# Save anchor and adrift read alignments.
saveRDS(dplyr::distinct(anchorReadAlignments), file.path(opt$outputDir, opt$alignReads_outputDir, 'anchorReadAlignments.rds'), compress = opt$compressDataFiles)  
saveRDS(dplyr::distinct(adriftReadAlignments), file.path(opt$outputDir, opt$alignReads_outputDir, 'adriftReadAlignments.rds'), compress = opt$compressDataFiles) 

unlink(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1'), recursive = TRUE)
unlink(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2'), recursive = TRUE)

if(any(! incomingSamples %in% anchorReadAlignments$uniqueSample) & opt$core_createFauxFragDoneFiles) core_createFauxFragDoneFiles()

q(save = 'no', status = 0, runLast = FALSE) 
