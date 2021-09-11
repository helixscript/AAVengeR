library(ShortRead)
library(dplyr)
library(parallel)
library(readr)
library(yaml)
options(stringsAsFactors = FALSE)

opt <- read_yaml('config.yml')

source(file.path(opt$softwareDir, 'lib.R'))

opt$inputFastaDir <- file.path(opt$outputDir, 'uniqueFasta')
opt$outputDir <- file.path(opt$outputDir, 'alignedReads')
dir.create(opt$outputDir)

samples <- read_tsv(opt$sampleConfigFile, col_types=cols())
if(nrow(samples) == 0) stop('Error - no lines of information was read from the sample configuration file.')
if(! 'subject' %in% names(samples))   samples$subject <- 'subject'
if(! 'replicate' %in% names(samples)) samples$replicate <- 1
if(any(grepl('~|\\|', paste(samples$subject, samples$sample, samples$replicate)))) stop('Error -- tildas (~) are reserved characters and can not be used in the subject, sample, or replicate sample configuration columns.')

if(! all(file.exists(file.path(opt$softwareDir, 'data', 'blatDBs', paste0(unique(samples$refGenome.id), '.2bit'))))) stop('All reference genomes are not available.')

samples$uniqueSample <- paste0(samples$subject, '~', samples$sample, '~', samples$replicate)
if(any(duplicated(samples$uniqueSample))) stop('Error -- all subject, sample, replicate id combinations are not unique.')

f <- list.files(opt$inputFastaDir, full.names = FALSE) 

dir.create(file.path(opt$outputDir, 'trimmedAnchorReads'))
dir.create(file.path(opt$outputDir, 'trimmedAnchorReadsLogs'))
dir.create(file.path(opt$outputDir, 'trimmedAdriftReads'))
dir.create(file.path(opt$outputDir, 'trimmedAdriftReadsLogs'))


cluster <- makeCluster(opt$alignReads_CPUs)
clusterExport(cluster, c('opt', 'samples'))

# Read in each anchor read file and trim the static linker reads may read into.
anchorReads <- Reduce('append', parLapply(cluster, f[grepl('anchorReads', f)], function(r){
#anchorReads <- Reduce('append', lapply(f[grepl('anchorReads', f)], function(r){
  source(file.path(opt$softwareDir, 'lib.R'))
  d <- subset(samples, uniqueSample == unlist(strsplit(r, '\\.'))[1])
  
  # Use the last 10 NTs of the common linker sequence as the adapter to be trimmed.
  anchorReadOverReadSeq <- substr(d$adriftRead.linker.seq, nchar(d$adriftRead.linker.seq) - 10, nchar(d$adriftRead.linker.seq))
  anchorReadOverReadSeq <- Biostrings::DNAString(anchorReadOverReadSeq)
  anchorReadOverReadSeq <- as.character(Biostrings::reverseComplement(anchorReadOverReadSeq))
  
  system(paste0(opt$command_cutadapt, ' -f fasta  -e 0.15 -a ', anchorReadOverReadSeq, ' --overlap 2 ',
                '--info-file=', file.path(opt$outputDir, 'trimmedAnchorReadsLogs', paste0(r, '.log')), ' ',
                file.path(opt$inputFastaDir, r), ' > ', file.path(opt$outputDir, 'trimmedAnchorReads', r)), ignore.stderr = TRUE)

   parseCutadaptLog(file.path(opt$outputDir, 'trimmedAnchorReadsLogs', paste0(r, '.log')))

   Biostrings::readDNAStringSet(file.path(opt$outputDir, 'trimmedAnchorReads', r))
}))

# Remove reads below the minimum anchor read length.
anchorReads <- anchorReads[width(anchorReads) >= opt$alignReads_minAnchorReadLengthPostTrim]


# Create a mapping of read ids to samples and reference genome.
readSampleMap <- bind_rows(lapply(f[grepl('anchorReads', f)], function(r){
  data.frame(sample = unlist(strsplit(r, '\\.'))[1], id = readFasta(opt$inputFastaDir, r)@id) 
})) %>% left_join(select(samples, uniqueSample, refGenome.id), by = c('sample' = 'uniqueSample'))



# Align the anchor reads to their reference genome(s).
align <- function(r, refGenome){
  clusterExport(cluster, 'refGenome', envir = environment())
  
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
  
  # Create a splitting vector for the reads which distribute varied read lengths across cores.
  r <- r[order(width(r))]
  i <- rep(1:opt$alignReads_CPUs, length(r))[1:length(r)]

  a <- bind_rows(parLapply(cluster, split(r, i), function(seqChunk){
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
  a
}

anchorReadAlignments <- 
  bind_rows(lapply(split(readSampleMap, readSampleMap$refGenome.id), function(x){
                refGenome <- file.path(opt$softwareDir, 'data', 'blatDBs', paste0(x$refGenome.id[1], '.2bit'))
                align(anchorReads[names(anchorReads) %in% x$id], refGenome)
               }))


# Select anchor reads where the ends align to the genome.
i <- (anchorReadAlignments$qSize - anchorReadAlignments$qEnd) <= opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned
anchorReadAlignments <- anchorReadAlignments[i,]


# Select anchor reads with alignments which do not start immediately since we expect some AAV remnant before the alignment.
i <- anchorReadAlignments$qStart > opt$alignReads_genomeAlignment_anchorReadMinStart
anchorReadAlignments <- anchorReadAlignments[i,]


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


# Remove leading linker sequences.
adriftReads <- Reduce('append', lapply(split(readSampleMap, readSampleMap$adriftRead.linker.length), function(x){
                 i <- x$adriftRead.linker.length[1] + 1
                 s <- adriftReads[names(adriftReads) %in% x$id]
                 s <- s[width(s) > i]
                 s <- subseq(s, i)
                 s[width(s) >= opt$alignReads_minAdriftReadLengthPostTrim]
               }))

# Update readSampleMap
readSampleMap <- subset(readSampleMap, id %in% names(adriftReads))


# Limit anchor read alignments to adrift read alignments ready to be aligned.
anchorReadAlignments <- subset(anchorReadAlignments, qName %in% names(adriftReads))


# Retrieve the 10NTs before anchor read alignments begin which will be used for adrift read over-reading trimming. 
adriftReadsAdapters <- bind_rows(lapply(split(anchorReadAlignments, anchorReadAlignments$qStart), function(x){
                         o <- anchorReads[names(anchorReads) %in% x$qName]
                         o <- o[match(x$qName, names(o))]
                         
                         tibble(id = names(o), 
                                leaderSeq = as.character(subseq(o, 1, x$qStart[1]-1)),
                                adapter = as.character(reverseComplement(subseq(o, x$qStart[1]-10, x$qStart[1]-1))))
                       }))


# Trim adrift read over-reading.
clusterExport(cluster, c('adriftReads', 'parseCutadaptLog', 'tmpFile'))
adriftReads <- Reduce('append', parLapply(cluster, split(adriftReadsAdapters, adriftReadsAdapters$adapter), function(x){
#adriftReads <- Reduce('append', lapply(split(adriftReadsAdapters, adriftReadsAdapters$adapter), function(x){
  f <- tmpFile()
  ShortRead::writeFasta(adriftReads[names(adriftReads) %in% x$id], file = file.path(opt$outputDir, 'trimmedAdriftReads', f))
  
  system(paste0(opt$command_cutadapt, ' -f fasta  -e 0.15 -a ', x$adapter[1], ' --overlap 2 ',
                '--info-file=', file.path(opt$outputDir, 'trimmedAdriftReadsLogs', paste0(x$adapter[1], '.log')), ' ',
                file.path(opt$outputDir, 'trimmedAdriftReads', f), ' > ', file.path(opt$outputDir, 'trimmedAdriftReads', paste0(f, '.trimmed'))), ignore.stderr = TRUE)

  parseCutadaptLog(file.path(opt$outputDir, 'trimmedAdriftReadsLogs', paste0(x$adapter[1], '.log')))
  Biostrings::readDNAStringSet(file.path(opt$outputDir, 'trimmedAdriftReads', paste0(f, '.trimmed')))
}))


# Remove post-trimmed adrift reads that fall below the min. length threshold.
adriftReads <- adriftReads[width(adriftReads) >= opt$alignReads_minAdriftReadLengthPostTrim]


# Update readSampleMap because it drives the alignment function.
readSampleMap <- subset(readSampleMap, id %in% names(adriftReads))


# Align adrift reads to the reference genome.
adriftReadAlignments <- 
  bind_rows(lapply(split(readSampleMap, readSampleMap$refGenome.id), function(x){
    refGenome <- file.path(opt$softwareDir, 'data', 'blatDBs', paste0(x$refGenome.id[1], '.2bit'))
    align(adriftReads[names(adriftReads) %in% x$id], refGenome)
  }))


# Select adrift reads where the ends align to the genome.
i <- (adriftReadAlignments$qSize - adriftReadAlignments$qEnd) <= opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned
adriftReadAlignments <- adriftReadAlignments[i,]


# Select adrift reads with alignments that start at the beginning of reads.
i <- adriftReadAlignments$qStart <= opt$alignReads_genomeAlignment_adriftReadMaxStart
adriftReadAlignments <- adriftReadAlignments[i,]


i <- base::intersect(anchorReadAlignments$qName, adriftReadAlignments$qName)
anchorReadAlignments <- anchorReadAlignments[anchorReadAlignments$qName %in% i,]
adriftReadAlignments <- adriftReadAlignments[adriftReadAlignments$qName %in% i,]

anchorReadAlignments <- left_join(anchorReadAlignments, select(adriftReadsAdapters, id, leaderSeq), by = c('qName' = 'id'))

anchorReadAlignments <- select(anchorReadAlignments, qName, strand, qSize, qStart, qEnd, tName, 
                               tStart, tEnd, queryPercentID, tAlignmentWidth, queryWidth, 
                               alignmentPercentID, percentQueryCoverage, leaderSeq)

adriftReadAlignments <- select(adriftReadAlignments, qName, strand, qSize, qStart, qEnd, tName, 
                              tStart, tEnd, queryPercentID, tAlignmentWidth, queryWidth, 
                              alignmentPercentID, percentQueryCoverage)

anchorReadAlignments <- left_join(anchorReadAlignments, select(readSampleMap, sample, id), by = c('qName' = 'id'))
adriftReadAlignments <- left_join(adriftReadAlignments, select(readSampleMap, sample, id), by = c('qName' = 'id'))


saveRDS(anchorReadAlignments, file.path(opt$outputDir, 'anchorReadAlignments.rds'))  
saveRDS(adriftReadAlignments, file.path(opt$outputDir, 'adriftReadAlignments.rds')) 



