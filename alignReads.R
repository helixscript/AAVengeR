library(lubridate)
library(ShortRead)
library(dplyr)
library(parallel)
library(readr)
options(stringsAsFactors = FALSE)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)

source(file.path(opt$softwareDir, 'lib.R'))

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

f <- list.files(opt$inputFastaDir, full.names = FALSE) 

cluster <- makeCluster(opt$alignReads_CPUs)
clusterExport(cluster, c('opt', 'samples'))

# Read in anchor reads.
anchorReads <- Reduce('append', parLapply(cluster, f[grepl('anchorReads', f)], function(r){
  Biostrings::readDNAStringSet(file.path(opt$inputFastaDir, r))  
}))


# Remove recognizable leader sequences from anchorReads so that they do not interfer with alignments
# and save the original sequences so that they can be reconstructed.

anchorReads_preLeaderTrim <- anchorReads

message('Trimming leader sequences')
m <- readRDS(file.path(opt$outputDir, opt$alignReads_mapLeaderSequencesOutputFile))
  
a <- anchorReads[names(anchorReads) %in% m$id]
b <- anchorReads[! names(anchorReads) %in% m$id]

m <- m[match(names(a), m$id),]
a <- a[width(a) >= (m$leaderMapping.qEnd + 1)]
m <- m[match(names(a), m$id),]
  
if(! all(names(a) == m$id)) stop('Error - could not align anchor reads to leader sequence mapping data.')

if(opt$alignReads_includeAnchorReadsWithoutMappings){
  anchorReads <- Reduce('append', list(subseq(a, (m$leaderMapping.qEnd + 1)), b))
}else{
  write(c(paste(now(), sprintf("%.2f%%", (length(b) / length(anchorReads))*100), 'reads without mappings removed.'), file = file.path(opt$outputDir, 'log'), append = TRUE))
  anchorReads <- subseq(a, (m$leaderMapping.qEnd + 1))
}


# Select anchor reads which are still as long as the required minimum length 
# after removing over read sequences and removing recognizable leader sequences.
anchorReads <- anchorReads[width(anchorReads) >= opt$alignReads_minAnchorReadLengthPostTrim]
  
if(length(anchorReads) == 0) stop('Error - no anchor reads remain after trimming leader sequences.')

  
  
# Create a mapping of read ids to samples and reference genome.
readSampleMap <- bind_rows(lapply(f[grepl('anchorReads', f)], function(r){
                   reads <- readDNAStringSet(file.path(opt$inputFastaDir, r))
                   reads <- anchorReads[names(anchorReads) %in% names(reads)]
                   
                   tibble(sample = unlist(strsplit(r, '\\.'))[1], 
                          id = names(reads), 
                          seq = as.character(reads)) 
                 })) %>% 
                 left_join(select(samples, uniqueSample, refGenome.id), by = c('sample' = 'uniqueSample'))


# Create Anchor read BLAT chunks such that each chunk contains equal number of reads and 
# equal numbers of short and long sequences. 

o <- unlist(lapply(split(readSampleMap, readSampleMap$refGenome.id), function(x){
       x <- x[! duplicated(x$seq),]
       x <- x[order(nchar(x$seq)),]
       a <- round(nrow(x) / 2)
       b <- suppressWarnings(unique(c(rbind(c(1:a),  rev(c((a+1):nrow(x)))))))
       x <- x[b,]
       
       chunkSize <- floor(nrow(x) / opt$alignReads_CPUs)
       if(opt$alignReads_alignmentChunkSize > 0) chunkSize <- opt$alignReads_alignmentChunkSize
       # if('alignReads_alignmentChunkSize' %in% names(opt)) chunkSize <- opt$alignReads_alignmentChunkSize
       
       split(x, (seq(nrow(x))-1) %/% chunkSize)
    }), recursive = FALSE)


# For each BLAT chunk, create FASTA files and run BLAT.
invisible(parLapply(cluster, o, function(x){
       f <-  file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1', 
                       paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')))

       write(paste0('>', x$id, '\n', x$seq), file = f)
       db <- file.path(opt$softwareDir, 'data', 'blatDBs', paste0(x[1,]$refGenome.id, '.2bit'))
       rm(x); gc()
       system(paste0(opt$command_blat, ' ', db, ' ', f, ' ', paste0(f, '.psl'), ' -tileSize=', opt$alignReads_genomeAlignment_blatTileSize, 
                     ' -stepSize=', opt$alignReads_genomeAlignment_blatStepSize, ' -minIdentity=', (opt$alignReads_genomeAlignment_minPercentID - 1), 
                     ' -out=psl -noHead'))
}))


# Parse and colate BLAT results.
# Here we filter on alignmentPercentID rather than % query alignment because we may not of removed all the not-genomic NTs from the anchor read.
b <- bind_rows(lapply(list.files(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat1'), pattern = '*.psl', full.names = TRUE), function(x){
       b <- parseBLAToutput(x)
       if(nrow(b) == 0) return(tibble())

       ### dplyr::filter(b, alignmentPercentID >= opt$alignReads_genomeAlignment_minPercentID, tNumInsert  <= 1, qNumInsert <= 1, tBaseInsert <= 2, qBaseInsert <= 2) %>%
       dplyr::filter(b, alignmentPercentID >= opt$alignReads_genomeAlignment_minPercentID, tNumInsert  <= 1, qNumInsert <= 1, tBaseInsert <= 2, qBaseInsert <= 2, qStart <= opt$alignReads_genomeAlignment_anchorRead_maxStartPos) %>%
       dplyr::select(qName, strand, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, queryPercentID, tAlignmentWidth, queryWidth, alignmentPercentID, percentQueryCoverage)
     }))

b <- left_join(b, select(readSampleMap, id, seq), by = c('qName' = 'id'))

anchorReadAlignments <- left_join(select(readSampleMap, id, seq), b, by = 'seq') 
anchorReadAlignments <- anchorReadAlignments[! is.na(anchorReadAlignments$qName),]
anchorReadAlignments <- select(anchorReadAlignments, -qName, -seq) %>% dplyr::rename(qName = id) %>% distinct()


# Select anchor reads where the ends align to the genome.
anchorReadAlignments$endDiff <- anchorReadAlignments$qSize - anchorReadAlignments$qEnd


i <- (anchorReadAlignments$qSize - anchorReadAlignments$qEnd) <= opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned
anchorReadAlignments <- anchorReadAlignments[i,]




# We need to reconstruct the leader sequences by determining what we removed before aligning and adding on additional
# NTs for those alignments which do not start at qStart = 1. Not All reads were trimmed going in because some mappings may of failed.

anchorReadAlignments <- left_join(anchorReadAlignments, tibble(id = names(anchorReads_preLeaderTrim), readSeq = as.character(anchorReads_preLeaderTrim)), by = c('qName' = 'id'))

anchorReadAlignments <- left_join(anchorReadAlignments, select(m, id, leaderMapping.qEnd), by = c('qName' = 'id')) 

a <- anchorReadAlignments[is.na(anchorReadAlignments$leaderMapping.qEnd),]
b <- anchorReadAlignments[! is.na(anchorReadAlignments$leaderMapping.qEnd),]

a$leaderSeq <- substr(a$readSeq, 1, a$qStart)
b$leaderSeq <- substr(b$readSeq, 1, (b$leaderMapping.qEnd + b$qStart))
  
anchorReadAlignments <- bind_rows(a, b)
anchorReadAlignments$leaderMapping.qEnd <- NULL

anchorReadAlignments$readSeq <-NULL


#~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

# Read in adrift reads and limit reads to found in the adrift reads alignment
adriftReads <- shortRead2DNAstringSet(Reduce('append', lapply(f[grepl('adriftReads', f)], function(r){
                 readFasta(file.path(opt$inputFastaDir, r))
               })))

# Limit adrift reads to those ids in the anchor reads.
adriftReads <- adriftReads[names(adriftReads) %in% anchorReadAlignments$qName]

# Limit anchor read alignments to those with leader sequences >= 10 NT because we need to create adrift over-read adapters.
# (!) This removed the possibility for keeping alignments that start the near the beginning of reads. 
anchorReadAlignments <- anchorReadAlignments[nchar(anchorReadAlignments$leaderSeq) >= opt$alignReads_genomeAlignment_anchorReadMinStart,]


# Retrieve the 10NTs before anchor read alignments begin which will be used for adrift read over-reading trimming. 
# Need to map by both read id and query start because a read can align with multiple start positions
adriftReadsAdapters <- tibble(id = anchorReadAlignments$qName, 
                              adapter = as.character(reverseComplement(DNAStringSet(substr(anchorReadAlignments$leaderSeq, nchar(anchorReadAlignments$leaderSeq)-9, nchar(anchorReadAlignments$leaderSeq))))))


# An anchor read may align to multiples positions in the genome each with a different start position (qStart).
# We need to trim all possibilities which will reduce the length of the adrift genomic letters in some cases. 

dupReadIds <- adriftReadsAdapters$id[duplicated(adriftReadsAdapters$id)]
adriftReadsAdapters.dups <- subset(adriftReadsAdapters, id %in% dupReadIds)

o <- bind_rows(lapply(split(adriftReadsAdapters.dups, adriftReadsAdapters.dups$id), function(x){
       tibble(id = x$id[1], adapter = paste0(sort(unique(x$adapter)), collapse = ','))
     }))

adriftReadsAdapters <- bind_rows(subset(adriftReadsAdapters, ! id %in% dupReadIds), o)

# Trim adrift read over-reading.
clusterExport(cluster, c('opt', 'samples', 'adriftReads', 'parseCutadaptLog', 'tmpFile'))

adriftReads <- Reduce('append', parLapply(cluster, split(adriftReadsAdapters, adriftReadsAdapters$adapter), function(x){
#adriftReads <- Reduce('append', lapply(split(adriftReadsAdapters, adriftReadsAdapters$adapter), function(x){
  f <- tmpFile()
  adapter <- paste0(paste0(' -a ', unlist(strsplit(x$adapter[1], ',')), ' '), collapse = ' ')
  
  ShortRead::writeFasta(adriftReads[names(adriftReads) %in% x$id], file = file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReads', f))
  
  system(paste0(opt$command_cutadapt, ' -f fasta  -e 0.15 ', adapter, ' --overlap 2 ',
               ### '--info-file=', file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReadsLogs', paste0(x$adapter[1], '.log')), ' ',
                file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReads', f), ' > ', file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReads', paste0(f, '.trimmed'))), ignore.stderr = TRUE)

  ### parseCutadaptLog(file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReadsLogs', paste0(x$adapter[1], '.log')))
  Biostrings::readDNAStringSet(file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReads', paste0(f, '.trimmed')))
}))


# Remove post-trimmed adrift reads that fall below the min. length threshold.
adriftReads <- adriftReads[width(adriftReads) >= opt$alignReads_minAdriftReadLengthPostTrim]


# Create a mapping of read ids to samples and reference genome.
readSampleMap <- bind_rows(lapply(f[grepl('adriftReads', f)], function(r){
  reads <- readDNAStringSet(file.path(opt$inputFastaDir, r))
  reads <- adriftReads[names(adriftReads) %in% names(reads)]
  
  tibble(sample = unlist(strsplit(r, '\\.'))[1], 
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
  # if('alignReads_alignmentChunkSize' %in% names(opt)) chunkSize <- opt$alignReads_alignmentChunkSize
  
  split(x, (seq(nrow(x))-1) %/% chunkSize)
}), recursive = FALSE)


# For each BLAT chunk, create FASTA files and run BLAT.
# Orignally developed with -tileSize=11 -stepSize=9
invisible(parLapply(cluster, o, function(x){
  f <-  file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2', 
                  paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')))
  
  write(paste0('>', x$id, '\n', x$seq), file = f)
  db <- file.path(opt$softwareDir, 'data', 'blatDBs', paste0(x[1,]$refGenome.id, '.2bit'))
  rm(x); gc()
  system(paste0(opt$command_blat, ' ', db, ' ', f, ' ', paste0(f, '.psl'), ' -tileSize=', opt$alignReads_genomeAlignment_blatTileSize, 
              ' -stepSize=', opt$alignReads_genomeAlignment_blatStepSize, ' -minIdentity=', (opt$alignReads_genomeAlignment_minPercentID - 1), 
              ' -out=psl -noHead'))
}))


stopCluster(cluster)

# Parse and colate BLAT results.
b <- bind_rows(lapply(list.files(file.path(opt$outputDir, opt$alignReads_outputDir, 'blat2'), pattern = '*.psl', full.names = TRUE), function(x){
  b <- parseBLAToutput(x)
  if(nrow(b) == 0) return(tibble())
  dplyr::filter(b, alignmentPercentID >= opt$alignReads_genomeAlignment_minPercentID, tNumInsert  <= 1, qNumInsert <= 1, tBaseInsert <= 2, qBaseInsert <= 2) %>%
  dplyr::select(qName, strand, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, queryPercentID, tAlignmentWidth, queryWidth, alignmentPercentID, percentQueryCoverage)
}))


b <- left_join(b, select(readSampleMap, id, seq), by = c('qName' = 'id'))
adriftReadAlignments <- left_join(select(readSampleMap, id, seq), b, by = 'seq') 
adriftReadAlignments <- adriftReadAlignments[! is.na(adriftReadAlignments$qName),]
adriftReadAlignments <- select(adriftReadAlignments, -qName, -seq) %>% dplyr::rename(qName = id) %>% distinct()


# Select adrift reads where the ends align to the genome.
i <- (adriftReadAlignments$qSize - adriftReadAlignments$qEnd) <= opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned
adriftReadAlignments <- adriftReadAlignments[i,]


# Select adrift reads with alignments that start at the beginning of reads.
i <- adriftReadAlignments$qStart <= opt$alignReads_genomeAlignment_adriftReadMaxStart
adriftReadAlignments <- adriftReadAlignments[i,]


i <- base::intersect(anchorReadAlignments$qName, adriftReadAlignments$qName)
anchorReadAlignments <- anchorReadAlignments[anchorReadAlignments$qName %in% i,]
adriftReadAlignments <- adriftReadAlignments[adriftReadAlignments$qName %in% i,]


anchorReadAlignments <- left_join(anchorReadAlignments, select(readSampleMap, sample, id), by = c('qName' = 'id'))
adriftReadAlignments <- left_join(adriftReadAlignments, select(readSampleMap, sample, id), by = c('qName' = 'id'))


saveRDS(anchorReadAlignments, file.path(opt$outputDir, opt$alignReads_outputDir, opt$alignReads_anchorReadAlignmentsOutputFile))  
saveRDS(adriftReadAlignments, file.path(opt$outputDir, opt$alignReads_outputDir, opt$alignReads_adriftReadAlignmentsOutputFile)) 

q(save = 'no', status = 0, runLast = FALSE) 
