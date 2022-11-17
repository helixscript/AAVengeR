library(lubridate)
library(ShortRead)
library(dplyr)
library(parallel)
library(readr)
library(data.table)
options(stringsAsFactors = FALSE)

# Read in configuration file.
configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))


# Create the required directory structure.
write(c(paste(lubridate::now(), '   Creating required directories.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
dir.create(file.path(opt$outputDir, opt$alignReads_outputDir))
dir.create(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1'))
dir.create(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2'))


# Load sample data and check to see if all requested blat DBs are available.
write(c(paste(now(), '   Reading sample data.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
samples <- loadSamples()

if(! all(file.exists(file.path(opt$softwareDir, 'data', 'blatDBs', paste0(unique(samples$refGenome.id), '.2bit'))))){
  write(c(paste(now(), 'Errror -- all reference genomes are not available.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}

# Create a CPU cluster.
cluster <- makeCluster(opt$alignReads_CPUs)
clusterExport(cluster, c('opt', 'samples', 'tmpFile'))

write(c(paste(now(), '   Reading in prepped reads.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

reads <- readRDS(file.path(opt$outputDir, opt$alignReads_inputFile))


blat <- function(y, ref, dir){
  f <- file.path(dir, tmpFile())
  write(paste0('>', y$id, '\n', y$seq), f)
  
  db <- file.path(opt$softwareDir, 'data', 'blatDBs', paste0(ref, '.2bit'))
  
  system(paste0(opt$command_blat, ' ', db, ' ', f, ' ', paste0(f, '.psl'), 
                ' -tileSize=', opt$alignReads_genomeAlignment_blatTileSize, 
                ' -stepSize=', opt$alignReads_genomeAlignment_blatStepSize, 
                ' -repMatch=', opt$alignReads_genomeAlignment_repMatch,
                ' -minIdentity=', (opt$alignReads_genomeAlignment_minPercentID - 1), 
                ' -out=psl -noHead'))
  
  write('done', paste0(f, '.done'))
}

anchorReadAlignments <- rbindlist(lapply(split(reads, reads$refGenome), function(x){
  s <- data.table(seq = unique(x$anchorReadSeq))
  s$id <- paste0('s', 1:nrow(s))
  
  s$cut <- cut(nchar(s$seq), c(-Inf, seq(0, max(nchar(s$seq)), by = 10), Inf), labels = FALSE)
  s <- group_by(s, cut) %>% mutate(n = ntile(1:n(), opt$alignReads_CPUs)) %>% ungroup()
  
  write(c(paste0(now(), '    Aligning ', nrow(s), ' anchor reads against ', x$refGenome[1], '.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
  dir <- file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1')
  invisible(parLapply(cluster, split(s, s$n), blat, x$refGenome[1], dir))
  
  b <- rbindlist(lapply(list.files(dir, pattern = '*.psl', full.names = TRUE), function(x){
         b <- data.table(parseBLAToutput(x))
         if(nrow(b) == 0) return(data.table())
    
         dplyr::filter(b, alignmentPercentID >= opt$alignReads_genomeAlignment_minPercentID, tNumInsert <= 1, 
                       qNumInsert <= 1, tBaseInsert <= 2, qBaseInsert <= 2, qStart <= opt$alignReads_genomeAlignment_anchorRead_maxStartPos) %>%
         dplyr::select(qName, strand, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, queryPercentID, tAlignmentWidth, queryWidth, alignmentPercentID, percentQueryCoverage)
  }))
  
  # Files need to be removed otherwise the join below will keep joining files from previous genomes.
  invisible(file.remove(list.files(dir, full.names = TRUE)))
  
  if(nrow(b) > 0){
    # Use the unique sequences to map alignments back to sequences 
    b <- left_join(b, distinct(dplyr::select(s, id, seq)), by = c('qName' = 'id'))
    x2 <- dplyr::select(x, uniqueSample, readID, refGenome, anchorReadSeq) %>% dplyr::filter(anchorReadSeq %in% b$seq)
    b <- left_join(x2, b, by = c('anchorReadSeq' = 'seq')) %>% dplyr::select(-anchorReadSeq, -qName)
  } 
  
  b
}))

# save.image('~/alignReadsDev.RData')


# Select anchor reads where the ends align to the genome.
i <- (anchorReadAlignments$qSize - anchorReadAlignments$qEnd) <= opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned

if(sum(i) == 0){
  write(c(paste(lubridate::now(), 'Error - no anchor read alignments remain after alignReads_genomeAlignment_anchorReadEnd_maxUnaligned filter.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE)
}

anchorReadAlignments <- anchorReadAlignments[i,]

reads <- subset(reads, readID %in% anchorReadAlignments$readID)

adriftReadAlignments <- rbindlist(lapply(split(reads, reads$refGenome), function(x){
  s <- data.table(seq = unique(x$adriftReadSeq))
  s$id <- paste0('s', 1:nrow(s))
  
  s$cut <- cut(nchar(s$seq), c(-Inf, seq(0, max(nchar(s$seq)), by = 10), Inf), labels = FALSE)
  s <- group_by(s, cut) %>% mutate(n = ntile(1:n(), opt$alignReads_CPUs)) %>% ungroup() %>% data.table()
  
  write(c(paste0(now(), '    Aligning ', nrow(s), ' adrift reads against ', x$refGenome[1], '.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
  dir <- file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2')
  invisible(parLapply(cluster, split(s, s$n), blat, x$refGenome[1], dir))
  
  b <- rbindlist(lapply(list.files(dir, pattern = '*.psl', full.names = TRUE), function(x){
    b <- data.table(parseBLAToutput(x))
    if(nrow(b) == 0) return(data.table())
    
    dplyr::filter(b, alignmentPercentID >= opt$alignReads_genomeAlignment_minPercentID, tNumInsert <= 1, 
                  qNumInsert <= 1, tBaseInsert <= 2, qBaseInsert <= 2, qStart <= opt$alignReads_genomeAlignment_adriftReadMaxStart) %>%
      dplyr::select(qName, strand, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, queryPercentID, tAlignmentWidth, queryWidth, alignmentPercentID, percentQueryCoverage)
  }))
  
  invisible(file.remove(list.files(dir, full.names = TRUE)))
  
  # Use the unique sequences to map alignments back to sequences 
  if(nrow(b) > 0){
    b <- left_join(b, distinct(dplyr::select(s, id, seq)), by = c('qName' = 'id'))
    x2 <- dplyr::select(x, uniqueSample, readID, refGenome, adriftReadSeq) %>% dplyr::filter(adriftReadSeq %in% b$seq)
    b <- left_join(x2, b, by = c('adriftReadSeq' = 'seq')) %>% dplyr::select(-adriftReadSeq, -qName)
  } 
  
  b
}))

# Select adrift reads where the ends align to the genome.
i <- (adriftReadAlignments$qSize - adriftReadAlignments$qEnd) <= opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned

if(sum(i) == 0){
  write(c(paste(lubridate::now(), 'Error - no adrift read alignments remain after alignReads_genomeAlignment_adriftReadEnd_maxUnaligned filter.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE)
}

adriftReadAlignments <- adriftReadAlignments[i,]


# Select adrift reads with alignments that start at the beginning of reads.
i <- adriftReadAlignments$qStart <= opt$alignReads_genomeAlignment_adriftReadMaxStart

if(sum(i) == 0){
  write(c(paste(lubridate::now(), 'Error - no adrift read alignments remain after alignReads_genomeAlignment_adriftReadMaxStart filter.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE)
}

adriftReadAlignments <- adriftReadAlignments[i,]


# Find the intersection between retained anchor and adrift reads and limit read pairs to pairs
# where both mates aligned within parameters.

write(c(paste(now(), '   Finding read pairs where both anchor and adrift read pair mates aligned well.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

i <- base::intersect(anchorReadAlignments$readID, adriftReadAlignments$readID)
anchorReadAlignments <- anchorReadAlignments[anchorReadAlignments$readID %in% i,]
adriftReadAlignments <- adriftReadAlignments[adriftReadAlignments$readID %in% i,]

anchorReadAlignments <- left_join(anchorReadAlignments, select(reads, readID, anchorReadSeq, leaderSeq), by = 'readID')
anchorReadAlignments$leaderSeq <- paste0(anchorReadAlignments$leaderSeq, substr(anchorReadAlignments$anchorReadSeq, 1, anchorReadAlignments$qStart))
anchorReadAlignments <- dplyr::select(anchorReadAlignments, -anchorReadSeq)

# Save anchor and adrift read alignments.
saveRDS(dplyr::distinct(anchorReadAlignments), file.path(opt$outputDir, opt$alignReads_outputDir, opt$alignReads_anchorReadAlignmentsOutputFile))  
saveRDS(dplyr::distinct(adriftReadAlignments), file.path(opt$outputDir, opt$alignReads_outputDir, opt$alignReads_adriftReadAlignmentsOutputFile)) 

unlink(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1'), recursive = TRUE)
unlink(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2'), recursive = TRUE)

q(save = 'no', status = 0, runLast = FALSE) 
