library(lubridate)
library(ShortRead)
library(dplyr)
library(parallel)
library(readr)
library(data.table)
options(stringsAsFactors = FALSE)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

write(c(paste(lubridate::now(), 'Creating required directories.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
dir.create(file.path(opt$outputDir, opt$alignReads_outputDir))
dir.create(file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReads'))
dir.create(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1'))
dir.create(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2'))
opt$inputFastaDir <- file.path(opt$outputDir, opt$alignReads_inputDir)

samples <- loadSamples()

if(! all(file.exists(file.path(opt$softwareDir, 'data', 'blatDBs', paste0(unique(samples$refGenome.id), '.2bit'))))){
  write(c(paste(now(), 'Errror -- all reference genomes are not available.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}


# Create a mapping of read ids to sample ids.
readID_to_sampleID <- bind_rows(lapply(list.files(opt$inputFastaDir, full.names = FALSE, pattern = 'anchor'), function(f){
                        r <- Biostrings::readDNAStringSet(file.path(opt$inputFastaDir, f))  
                        tibble(sample = unlist(strsplit(f, '\\.'))[1], readID = names(r))
                      }))

# Create a CPU cluster.
cluster <- makeCluster(opt$alignReads_CPUs)
clusterExport(cluster, c('opt', 'samples'))


# Read in anchor reads.
write(c(paste(lubridate::now(), 'Reading anchor reads.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
anchorReads <- Reduce('append', parLapply(cluster, list.files(opt$inputFastaDir, full.names = FALSE, pattern = 'anchor'), function(r){
  Biostrings::readDNAStringSet(file.path(opt$inputFastaDir, r))  
}))


# Read in adrift reads.
adriftReads <- shortRead2DNAstringSet(Reduce('append', lapply(list.files(opt$inputFastaDir, full.names = FALSE, pattern = 'adrift'), function(r){
  readFasta(file.path(opt$inputFastaDir, r))
})))


# Remove recognizable leader sequences from anchorReads so that they do not 
# interfere with alignments and save the original sequences so that they can 
# be reconstructed.


write(c(paste(lubridate::now(), 'Removing leader sequences.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

m <- readRDS(file.path(opt$outputDir, opt$alignReads_mapLeaderSequencesOutputFile))
  
# Select anchor reads with ids found in the leader sequence table.
a <- anchorReads[names(anchorReads) %in% m$id]

# Present...

# Arrange leader sequence table to match anchor read ids.
m <- m[match(names(a), m$id),]
if(! all(names(a) == m$id)) stop('Error - could not align anchor reads to leader sequence mapping data.')

# Select anchor reads that are long enough to trim and still have enough NTs to align.
a <- a[width(a) >= (m$leaderMapping.qEnd + opt$alignReads_minAnchorReadLengthPostTrim)]
m <- m[m$id %in% names(a),]
if(! all(names(a) == m$id)) stop('Error - could not align anchor reads to leader sequence mapping data.')

# Trim anchor reads for alignment.
anchorReads_leaderSeqsTrimmed <- subseq(a, m$leaderMapping.qEnd+1, width(a))
anchorReads_leaderSeqs <- subseq(a, 1, m$leaderMapping.qEnd)
if(length(anchorReads_leaderSeqsTrimmed) == 0) stop('Error - no anchor reads remain after trimming leader sequences.')


# Limit reads to those with anchor mates that can be aligned.
anchorReads <- anchorReads[names(anchorReads) %in% names(anchorReads_leaderSeqsTrimmed)]
adriftReads <- adriftReads[names(adriftReads) %in% names(anchorReads_leaderSeqsTrimmed)]

if(! all(names(anchorReads) == names(adriftReads))) stop('Read organization Error 1')
if(! all(names(anchorReads) == m$id)) stop('Read organization Error 2')

rm(a); gc()

# Trim adrift read for over-reading.
# This trimming may not be complete since we can only predict LTR/ITR rearrangements to a degree.
# We do not use parLapply here because the adriftReads is large and it is not optimal to export it to each node.
o <- tibble(id = names(anchorReads), seq = as.character(reverseComplement(subseq(anchorReads, m$leaderMapping.qStart, m$leaderMapping.qEnd))))
o$seqLen <- nchar(o$seq)
o$seq2 <- substr(o$seq, 1, ifelse(o$seqLen > 15, 15, o$seqLen))

adriftReads_overReadingTrimmed <- Reduce('append', lapply(split(o, o$seq2), function(x){
                 f <- tmpFile()
                 adapter <- paste0(' -a ', x$seq2[1])
  
                 ShortRead::writeFasta(adriftReads[names(adriftReads) %in% x$id], file = file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReads', f))
  
                 system(paste0(opt$command_cutadapt, ' -e 0.15 ', adapter, ' --overlap 2 ',
                        file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReads', f), ' > ', 
                        file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReads', paste0(f, '.trimmed'))), ignore.stderr = TRUE)
  
                 r <- Biostrings::readDNAStringSet(file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReads', paste0(f, '.trimmed')))
                 r[width(r) >= opt$alignReads_minAdriftReadLengthPostTrim]
               }))

# Sync the trimmed anchor reads to the trimmed adrift reads.
anchorReads_leaderSeqsTrimmed  <- anchorReads_leaderSeqsTrimmed[names(anchorReads_leaderSeqsTrimmed) %in% names(adriftReads_overReadingTrimmed)]
adriftReads_overReadingTrimmed <- adriftReads_overReadingTrimmed[match(names(anchorReads_leaderSeqsTrimmed), names(adriftReads_overReadingTrimmed))]
anchorReads_leaderSeqs <- anchorReads_leaderSeqs[names(anchorReads_leaderSeqs) %in% names(anchorReads_leaderSeqsTrimmed)]
unlink(file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReads'), recursive = TRUE)

rm(adriftReads); gc()

# Create a mapping of read ids to samples and reference genome.
# The readSampleMap object contains sample names, read ids, read sequences, 
# and reference genomes to align against.

write(c(paste(lubridate::now(), 'Creating anchor read map.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

readSampleMap <- bind_rows(lapply(split(readID_to_sampleID, readID_to_sampleID$sample), function(r){
                   reads <- anchorReads_leaderSeqsTrimmed[names(anchorReads_leaderSeqsTrimmed) %in% r$readID]
                   if(length(reads) == 0) return(tibble())
                   tibble(sample = r$sample[1],
                          id = names(reads), 
                          seq = as.character(reads)) 
                 })) %>% left_join(select(samples, uniqueSample, refGenome.id), by = c('sample' = 'uniqueSample'))



# Create anchor read BLAT chunks such that each chunk contains equal number of reads and 
# equal numbers of short and long sequences. 

write(c(paste(lubridate::now(), 'Creating mixed anchor read chunks.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

o <- unlist(lapply(split(readSampleMap, readSampleMap$refGenome.id), function(x){
       x <- x[! duplicated(x$seq),]
       x <- x[order(nchar(x$seq)),]
       a <- round(nrow(x) / 2)
       b <- suppressWarnings(unique(c(rbind(c(1:a),  rev(c((a+1):nrow(x)))))))
       x <- x[b,]
       
       chunkSize <- floor(nrow(x) / opt$alignReads_CPUs)
       if(opt$alignReads_alignmentChunkSize > 0) chunkSize <- opt$alignReads_alignmentChunkSize
       
       split(x, (seq(nrow(x))-1) %/% chunkSize)
    }), recursive = FALSE)


invisible(parLapply(cluster, o, function(x){
       f <-  file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1', 
                       paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')))

       write(paste0('>', x$id, '\n', x$seq), file = f)
       db <- file.path(opt$softwareDir, 'data', 'blatDBs', paste0(x[1,]$refGenome.id, '.2bit'))

       rm(x)
       gc()
       
       system(paste0(opt$command_blat, ' ', db, ' ', f, ' ', paste0(f, '.psl'), 
                      ' -tileSize=', opt$alignReads_genomeAlignment_blatTileSize, 
                      ' -stepSize=', opt$alignReads_genomeAlignment_blatStepSize, 
                      ' -repMatch=', opt$alignReads_genomeAlignment_repMatch,
                      ' -minIdentity=', (opt$alignReads_genomeAlignment_minPercentID - 1), 
                      ' -out=psl -noHead'))
}))


# #---
# b <- bind_rows(lapply(list.files(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1'), pattern = '*.psl', full.names = TRUE), function(x){
#        parseBLAToutput(x)
# }))
# 
# #---

# Parse and colate BLAT results.
# Here we filter on alignmentPercentID rather than % query alignment because we may not of removed all the not-genomic NTs from the anchor read.
b <- bind_rows(lapply(list.files(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1'), pattern = '*.psl', full.names = TRUE), function(x){
        b <- parseBLAToutput(x)
        if(nrow(b) == 0) return(tibble())

        dplyr::filter(b, alignmentPercentID >= opt$alignReads_genomeAlignment_minPercentID, tNumInsert <= 1, 
                      qNumInsert <= 1, tBaseInsert <= 2, qBaseInsert <= 2, qStart <= opt$alignReads_genomeAlignment_anchorRead_maxStartPos) %>%
        dplyr::select(qName, strand, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, queryPercentID, tAlignmentWidth, queryWidth, alignmentPercentID, percentQueryCoverage)
     }))


# Expand the unique sequence blat result to include duplicate reads from readSampleMap.
# This step-wise expansion is slow -- designed not to overwhelm other approaches such as joins.
b2 <- left_join(b, select(readSampleMap, id, seq, refGenome.id), by = c('qName' = 'id'))

anchorReadAlignments <- dplyr::distinct(bind_rows(lapply(split(b2, b2$refGenome.id), function(x){
  o <- dplyr::filter(readSampleMap, refGenome.id == x$refGenome.id[1])
  z  <- full_join(x, o, by = 'seq')
  z$qName <- z$id
  z$refGenome.id <- z$refGenome.id.x
  tidyr::drop_na(dplyr::select(z, -sample, -id, -refGenome.id.y, -refGenome.id.x))
})))
     


# Select anchor reads where the ends align to the genome.
i <- (anchorReadAlignments$qSize - anchorReadAlignments$qEnd) <= opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned

if(sum(i) == 0){
  write(c(paste(lubridate::now(), 'Error - no anchor read alignments remain after alignReads_genomeAlignment_anchorReadEnd_maxUnaligned filter.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE)
}

anchorReadAlignments <- anchorReadAlignments[i,]


anchorReads_leaderSeqsTrimmed  <- anchorReads_leaderSeqsTrimmed[names(anchorReads_leaderSeqsTrimmed) %in% anchorReadAlignments$qName]
adriftReads_overReadingTrimmed <- adriftReads_overReadingTrimmed[names(adriftReads_overReadingTrimmed) %in% anchorReadAlignments$qName]


write(c(paste(lubridate::now(), 'Creating adrift read map.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

readSampleMap <- bind_rows(lapply(split(readID_to_sampleID, readID_to_sampleID$sample), function(r){
  reads <- adriftReads_overReadingTrimmed[names(adriftReads_overReadingTrimmed) %in% r$readID]
  if(length(reads) == 0) return(tibble())
  tibble(sample = r$sample[1],
         id = names(reads), 
         seq = as.character(reads)) 
})) %>% left_join(select(samples, uniqueSample, refGenome.id), by = c('sample' = 'uniqueSample'))


o <- unlist(lapply(split(readSampleMap, readSampleMap$refGenome.id), function(x){
  x <- x[! duplicated(x$seq),]
  x <- x[order(nchar(x$seq)),]
  a <- round(nrow(x) / 2)
  b <- suppressWarnings(unique(c(rbind(c(1:a),  rev(c((a+1):nrow(x)))))))
  x <- x[b,]
  
  chunkSize <- floor(nrow(x) / opt$alignReads_CPUs)
  if(opt$alignReads_alignmentChunkSize > 0) chunkSize <- opt$alignReads_alignmentChunkSize
  
  split(x, (seq(nrow(x))-1) %/% chunkSize)
}), recursive = FALSE)


invisible(parLapply(cluster, o, function(x){
  f <-  file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2', 
                  paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')))
  
  write(paste0('>', x$id, '\n', x$seq), file = f)
  db <- file.path(opt$softwareDir, 'data', 'blatDBs', paste0(x[1,]$refGenome.id, '.2bit'))
  
  rm(x)
  gc()
  
  system(paste0(opt$command_blat, ' ', db, ' ', f, ' ', paste0(f, '.psl'), 
                ' -tileSize=', opt$alignReads_genomeAlignment_blatTileSize, 
                ' -stepSize=', opt$alignReads_genomeAlignment_blatStepSize, 
                ' -repMatch=', opt$alignReads_genomeAlignment_repMatch,
                ' -minIdentity=', (opt$alignReads_genomeAlignment_minPercentID - 1), 
                ' -out=psl -noHead'))
}))


b <- bind_rows(lapply(list.files(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2'), pattern = '*.psl', full.names = TRUE), function(x){
  b <- parseBLAToutput(x)
  if(nrow(b) == 0) return(tibble())
  
  dplyr::filter(b, alignmentPercentID >= opt$alignReads_genomeAlignment_minPercentID, tNumInsert <= 1, 
                qNumInsert <= 1, tBaseInsert <= 2, qBaseInsert <= 2, qStart <= opt$alignReads_genomeAlignment_anchorRead_maxStartPos) %>%
    dplyr::select(qName, strand, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, queryPercentID, tAlignmentWidth, queryWidth, alignmentPercentID, percentQueryCoverage)
}))


# Expand the unique sequence blat result to include duplicate reads from readSampleMap.
# This step-wise expansion is slow -- designed not to overwhelm other approaches such as joins.
b2 <- left_join(b, select(readSampleMap, id, seq, refGenome.id), by = c('qName' = 'id'))

adriftReadAlignments <- dplyr::distinct(bind_rows(lapply(split(b2, b2$refGenome.id), function(x){
  o <- dplyr::filter(readSampleMap, refGenome.id == x$refGenome.id[1])
  z  <- full_join(x, o, by = 'seq')
  z$qName <- z$id
  z$refGenome.id <- z$refGenome.id.x
  tidyr::drop_na(dplyr::select(z, -sample, -id, -refGenome.id.y, -refGenome.id.x))
})))


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

i <- base::intersect(anchorReadAlignments$qName, adriftReadAlignments$qName)

if(length(i) == 0){
  write(c(paste(lubridate::now(), 'Error - no shared read ids found between anchor and adrift read alignments.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE)
}

anchorReadAlignments <- anchorReadAlignments[anchorReadAlignments$qName %in% i,]
adriftReadAlignments <- adriftReadAlignments[adriftReadAlignments$qName %in% i,]

anchorReadAlignments <- dplyr::select(anchorReadAlignments, -seq)
adriftReadAlignments <- dplyr::select(adriftReadAlignments, -seq)

anchorReadAlignments <- left_join(anchorReadAlignments, readID_to_sampleID, by = c('qName' = 'readID'))
adriftReadAlignments <- left_join(adriftReadAlignments, readID_to_sampleID, by = c('qName' = 'readID'))

anchorReadAlignments <- left_join(anchorReadAlignments, select(m, id, leaderMapping.qEnd), by = c('qName' = 'id'))

s <- anchorReads[names(anchorReads) %in% anchorReadAlignments$qName]
x <- tibble(qName = names(s), seq = as.character(s))
rm(s, anchorReads)
gc()

anchorReadAlignments <- left_join(anchorReadAlignments, x, by = 'qName')
anchorReadAlignments$leaderSeq <- substr(anchorReadAlignments$seq, 1, anchorReadAlignments$leaderMapping.qEnd)

saveRDS(dplyr::distinct(dplyr::select(anchorReadAlignments, -seq)), file.path(opt$outputDir, opt$alignReads_outputDir, opt$alignReads_anchorReadAlignmentsOutputFile))  
saveRDS(dplyr::distinct(adriftReadAlignments), file.path(opt$outputDir, opt$alignReads_outputDir, opt$alignReads_adriftReadAlignmentsOutputFile)) 

if(! opt$alignReads_keepIntermediateFiles){
  unlink(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1'), recursive = TRUE)
  unlink(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2'), recursive = TRUE)
}

q(save = 'no', status = 0, runLast = FALSE) 
