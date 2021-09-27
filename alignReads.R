library(ShortRead)
library(dplyr)
library(parallel)
library(readr)
options(stringsAsFactors = FALSE)

opt <- yaml::read_yaml('config.yml')
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$alignReads_outputDir))
dir.create(file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReads'))

opt$inputFastaDir <- file.path(opt$outputDir, opt$alignReads_inputDir)

samples <- loadSamples()

if(! all(file.exists(file.path(opt$softwareDir, 'data', 'blatDBs', paste0(unique(samples$refGenome.id), '.2bit'))))) stop('All reference genomes are not available.')

f <- list.files(opt$inputFastaDir, full.names = FALSE) 

cluster <- makeCluster(opt$alignReads_CPUs)
clusterExport(cluster, c('opt', 'samples'))

# Trim the reverse complement of the last 10 NT of the static linker sequence from the end of anchor reads.
anchorReads <- Reduce('append', parLapply(cluster, f[grepl('anchorReads', f)], function(r){
  Biostrings::readDNAStringSet(file.path(opt$inputFastaDir, r))  
}))

stopCluster(cluster)


# Remove reads that aligned to the vector if vector alignments were performed.

if('alignReads_vectorFilterOutputFile' %in% names(opt)){ 
  if(file.exists(file.path(opt$outputDir, opt$alignReads_vectorFilterOutputFile))){
    v <- readRDS(file.path(opt$outputDir, opt$alignReads_vectorFilterOutputFile))
    if(nrow(v) > 0) anchorReads <- anchorReads[! names(anchorReads) %in% v$qname]
  }
}


# Remove recognizable leader sequences from anchorReads so that they do not interfer with alignments.

anchorReads_preLeaderTrim <- anchorReads

if('alignReads_mapLeaderSequencesOutputFile' %in% names(opt) & file.exists(file.path(opt$outputDir, opt$alignReads_mapLeaderSequencesOutputFile))){
  message('Trimming leader sequences')
  m <- readRDS(file.path(opt$outputDir, opt$alignReads_mapLeaderSequencesOutputFile))
  
  a <- anchorReads[names(anchorReads) %in% m$id]
  b <- anchorReads[! names(anchorReads) %in% m$id]

  m <- m[match(names(a), m$id),]
  a <- a[width(a) >= (m$leaderMapping.qEnd + 1)]
  m <- m[match(names(a), m$id),]
  
  if(! all(names(a) == m$id)) stop('Error - could not align anchor reads to leader sequence mapping data.')
  
  anchorReads <- Reduce('append', list(subseq(a, (m$leaderMapping.qEnd + 1)), b))

  # Select anchor reads which are still as long as the required minimum length 
  # after removing over read sequences and removing recognizable leader sequences.
  anchorReads <- anchorReads[width(anchorReads) >= opt$alignReads_minAnchorReadLengthPostTrim]
  
  if(length(anchorReads) == 0) stop('Error - no anchor reads remain after trimming leader sequences.')
}


# Create a mapping of read ids to samples and reference genome.
readSampleMap <- bind_rows(lapply(f[grepl('anchorReads', f)], function(r){
                   tibble(sample = unlist(strsplit(r, '\\.'))[1], 
                          id = as.character(readFasta(opt$inputFastaDir, r)@id)) 
                 })) %>% 
                 dplyr::filter(id %in% names(anchorReads)) %>%
                 left_join(select(samples, uniqueSample, refGenome.id), by = c('sample' = 'uniqueSample'))



# Align the anchor reads to their reference genome(s).
align <- function(r, refGenome, CPUs, chunkSize){
  cluster <- makeCluster(CPUs)
  clusterExport(cluster, c('opt', 'samples', 'refGenome'), envir = environment())
  
  # Create short tmp read ids and then group reads by sequence and assign grouped reads comma delimited ids.
  # Grouped tmp ids will be expanded after alignments are complete and reset to their true values.
  d <- tibble(id = names(r), seq = as.character(r)) %>%
       group_by(seq) %>% 
       summarise(idList = list(id)) %>% 
       ungroup() %>%
       mutate(id2 = paste0('s', 1:n()))
  
  r <- DNAStringSet(d$seq)
  names(r) <- d$id2
  gc()
  
  a <- bind_rows(parLapply(cluster, mixAndChunkSeqs(r, chunkSize), function(seqChunk){
         library(dplyr)
         library(ShortRead)
         source(file.path(opt$softwareDir, 'lib.R'))
         alignReads.BLAT(seqChunk, refGenome) %>%
         dplyr::filter(alignmentPercentID >= opt$alignReads_genomeAlignment_minPercentID, tNumInsert  <= 1, qNumInsert <= 1, tBaseInsert <= 2, qBaseInsert <= 2)
      }))
  
  i <- match(a$qName, d$id2)
  a$idList <- d[i,]$idList
  a <- tidyr::unnest(a, idList)
  a$qName <- a$idList 
  a$idList <- NULL
  stopCluster(cluster)
  a
}


# Align the leader sequence trimmed anchor reads to the genome.
#
# Tip: adjust opt$alignReads_alignmentChunkSize so that the number of unique reads / chunk size equals number of CPUs.
# Memory limitations may require less cpus or smaller chunk sizes.  2500 chunks across 25 CPUs use ~ 100GB of memory for CanFam3.
anchorReadAlignments <- 
  bind_rows(lapply(split(readSampleMap, readSampleMap$refGenome.id), function(x){
                refGenome <- file.path(opt$softwareDir, 'data', 'blatDBs', paste0(x$refGenome.id[1], '.2bit'))
                align(anchorReads[names(anchorReads) %in% x$id], refGenome, opt$alignReads_CPUs, opt$alignReads_alignmentChunkSize)
               }))

saveRDS(anchorReadAlignments, file.path(opt$outputDir, opt$alignReads_outputDir, 'anchorReadAlignments.save1.rds'))
anchorReadAlignments <- select(anchorReadAlignments, qName, strand, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, queryPercentID, tAlignmentWidth, queryWidth, alignmentPercentID, percentQueryCoverage)

# Select anchor reads where the ends align to the genome.
i <- (anchorReadAlignments$qSize - anchorReadAlignments$qEnd) <= opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned
anchorReadAlignments <- anchorReadAlignments[i,]


if(opt$alignReads_performSecondGenomeAlignment){
   # Some reads with have qStarts of zero which means they start on the first NT.
   # This may be after recognizable sequences were removed if mappings were performed earlier. 
   # Aligning without removing mappings (or even when we do) may get us close to the true boundary but non-native 
   # NTs may of shift the alignemnt. Here we remove the extra unaligned 5' NTs after the first alignment 
   # and align again which may get us closer to the true boundaries.
   #
   # Second alignments
   # --------------------|~~~~~~~~~~~~~~~~
   # --|~~~~~~~~~~~~~~~~~~~~~~~~~
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
}




# We need to reconstruct the leader sequencs by determining what we removed before aligning and adding on additional
# NTs for those alignments which do not start at qStart = 1. Not All reads were trimmed going in because some mappings may of failed.

anchorReadAlignments <- left_join(anchorReadAlignments, tibble(id = names(anchorReads_preLeaderTrim), readSeq = as.character(anchorReads_preLeaderTrim)), by = c('qName' = 'id'))

if('alignReads_mapLeaderSequencesOutputFile' %in% names(opt) & file.exists(file.path(opt$outputDir, opt$alignReads_mapLeaderSequencesOutputFile))){
  message('reconstructing leader sequences')
  anchorReadAlignments <- left_join(anchorReadAlignments, select(m, id, leaderMapping.qEnd), by = c('qName' = 'id')) 

  a <- anchorReadAlignments[is.na(anchorReadAlignments$leaderMapping.qEnd),]
  b <- anchorReadAlignments[! is.na(anchorReadAlignments$leaderMapping.qEnd),]
  
  a$leaderSeq <- substr(a$readSeq, 1, a$qStart)
  b$leaderSeq <- substr(b$readSeq, 1, (b$leaderMapping.qEnd + b$qStart))
  anchorReadAlignments <- bind_rows(a, b)
  anchorReadAlignments$leaderMapping.qEnd <- NULL
} else {
  anchorReadAlignments$leaderSeq <- substr(anchorReadAlignments$readSeq, 1, anchorReadAlignments$qStart)
}
  
anchorReadAlignments$readSeq <-NULL


# Read in adrift reads and limit reads to found in the adrift reads alignment
adriftReads <- shortRead2DNAstringSet(Reduce('append', lapply(f[grepl('adriftReads', f)], function(r){
                 readFasta(file.path(opt$inputFastaDir, r))
               })))


# Limit adrift reads to those ids in the anchor reads.
adriftReads <- adriftReads[names(adriftReads) %in% anchorReadAlignments$qName]


# Limit readSampleMap to those reads found in the adrift reads previously limited to anchor
# read ids that aligned well to the genome.
readSampleMap <- subset(readSampleMap, id %in% names(adriftReads)) %>% 
                 left_join(select(samples, uniqueSample, adriftRead.linker.seq), by = c('sample' = 'uniqueSample')) %>%
                 mutate(adriftRead.linker.length = nchar(adriftRead.linker.seq))


# Limit anchor read alignments to those with leader sequencs >= 10 NT because we need to create adrift over-read adapters.
anchorReadAlignments <- anchorReadAlignments[nchar(anchorReadAlignments$leaderSeq) >= opt$alignReads_genomeAlignment_anchorReadMinStart,]


# Retrieve the 10NTs before anchor read alignments begin which will be used for adrift read over-reading trimming. 
# Need to map by both read id and query start because a read can align with multiple start positions
adriftReadsAdapters <- tibble(id = anchorReadAlignments$qName, 
                              adapter = as.character(reverseComplement(DNAStringSet(substr(anchorReadAlignments$leaderSeq, nchar(anchorReadAlignments$leaderSeq)-9, nchar(anchorReadAlignments$leaderSeq))))))

# An anchor read may align to multiples positions in the genome each with a different start position (qStart).
# We need to trim all posibilities which will reduce the length of the adrift genomic letters in some cases. 

dupReadIds <- adriftReadsAdapters$id[duplicated(adriftReadsAdapters$id)]
adriftReadsAdapters.dups <- subset(adriftReadsAdapters, id %in% dupReadIds)

o <- bind_rows(lapply(split(adriftReadsAdapters.dups, adriftReadsAdapters.dups$id), function(x){
       tibble(id = x$id[1], adapter = paste0(sort(unique(x$adapter)), collapse = ','))
     }))

adriftReadsAdapters <- bind_rows(subset(adriftReadsAdapters, ! id %in% dupReadIds), o)



# Trim adrift read over-reading.
cluster <- makeCluster(opt$alignReads_CPUs)
clusterExport(cluster, c('opt', 'samples', 'adriftReads', 'parseCutadaptLog', 'tmpFile'))


adriftReads <- Reduce('append', parLapply(cluster, split(adriftReadsAdapters, adriftReadsAdapters$adapter), function(x){
#adriftReads <- Reduce('append', lapply(split(adriftReadsAdapters, adriftReadsAdapters$adapter), function(x){
  f <- tmpFile()
  adapter <- paste0(paste0(' --no-indels  -a ', unlist(strsplit(x$adapter[1], ',')), ' '), collapse = ' ')
  
  ShortRead::writeFasta(adriftReads[names(adriftReads) %in% x$id], file = file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReads', f))
  
  system(paste0(opt$command_cutadapt, ' -f fasta  -e 0.15 ', adapter, ' --overlap 2 ',
               ### '--info-file=', file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReadsLogs', paste0(x$adapter[1], '.log')), ' ',
                file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReads', f), ' > ', file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReads', paste0(f, '.trimmed'))), ignore.stderr = TRUE)

  ### parseCutadaptLog(file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReadsLogs', paste0(x$adapter[1], '.log')))
  Biostrings::readDNAStringSet(file.path(opt$outputDir, opt$alignReads_outputDir, 'trimmedAdriftReads', paste0(f, '.trimmed')))
}))

stopCluster(cluster)

# Remove post-trimmed adrift reads that fall below the min. length threshold.
adriftReads <- adriftReads[width(adriftReads) >= opt$alignReads_minAdriftReadLengthPostTrim]


# Update readSampleMap because it drives the alignment function.
readSampleMap <- subset(readSampleMap, id %in% names(adriftReads))


# Align adrift reads to the reference genome.
adriftReadAlignments <- 
  bind_rows(lapply(split(readSampleMap, readSampleMap$refGenome.id), function(x){
    refGenome <- file.path(opt$softwareDir, 'data', 'blatDBs', paste0(x$refGenome.id[1], '.2bit'))
    align(adriftReads[names(adriftReads) %in% x$id], refGenome, opt$alignReads_CPUs, opt$alignReads_alignmentChunkSize)
  }))

saveRDS(adriftReadAlignments, file.path(opt$outputDir, opt$alignReads_outputDir, 'adriftReadAlignments1.rds'))

adriftReadAlignments <- select(adriftReadAlignments, qName, strand, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, queryPercentID, tAlignmentWidth, queryWidth, alignmentPercentID, percentQueryCoverage)


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
