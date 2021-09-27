library(dplyr)
library(parallel)

opt <- yaml::read_yaml('config.yml')
source(file.path(opt$softwareDir, 'lib.R'))

dups <- tibble()
if('buildFragments_duplicateReadFile' %in% names(opt) & file.exists(file.path(opt$outputDir, opt$buildFragments_duplicateReadFile))){
  dups <- select(readr::read_tsv(file.path(opt$outputDir, opt$buildFragments_duplicateReadFile)), id, n)
  dups <- dups[dups$n > 0,]
}

anchorReadAlignments <- readRDS(file.path(opt$outputDir, opt$buildFragments_anchorReadsAlignmentFile))
adriftReadAlignments <- readRDS(file.path(opt$outputDir, opt$buildFragments_adriftReadsAlignmentFile))

dir.create(file.path(opt$outputDir, opt$buildFragments_outputDir))
dir.create(file.path(opt$outputDir, opt$buildFragments_outputDir, 'multiHitReads'))


anchorReadAlignments <- select(anchorReadAlignments, sample, qName, tName, strand, tStart, tEnd, leaderSeq)
adriftReadAlignments <- select(adriftReadAlignments, sample, qName, tName, strand, tStart, tEnd)

a <- table(anchorReadAlignments$qName)
b <- table(adriftReadAlignments$qName)
x <- unique(c(names(a[a > opt$buildFragments_maxAlignmentsPerRead]), names(b[b > opt$buildFragments_maxAlignmentsPerRead])))

write(x, file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'highAlignmentReads'))

anchorReadAlignments <- subset(anchorReadAlignments, ! qName %in% x)
adriftReadAlignments <- subset(adriftReadAlignments, ! qName %in% x)

names(anchorReadAlignments) <- paste0(names(anchorReadAlignments), '.anchorReads')
names(adriftReadAlignments) <- paste0(names(adriftReadAlignments), '.adriftReads')

ids <- unique(anchorReadAlignments$qName.anchorReads)
id_groups <- split(ids, dplyr::ntile(1:length(ids), round(length(ids)/100)))

cluster <- makeCluster(opt$buildFragments_CPUs)
clusterExport(cluster, c('opt', 'dups', 'anchorReadAlignments', 'adriftReadAlignments'))

frags <- bind_rows(parLapply(cluster, id_groups, function(id_group){
#frags <- bind_rows(lapply(id_groups, function(id_group){
  library(dplyr)
  source(file.path(opt$softwareDir, 'lib.R'))
  
  a <- subset(anchorReadAlignments, qName.anchorReads %in% id_group)
  b <- subset(adriftReadAlignments, qName.adriftReads %in% id_group)
  
  if(nrow(a) == 0 | nrow(b) == 0) return(data.frame())
  
  frags <- left_join(a, b, by = c('qName.anchorReads' = 'qName.adriftReads')) %>% tidyr::drop_na()
  
  # Remove combinations not found on the same chromosome. 
  i <- which(frags$tName.anchorReads != frags$tName.adriftReads)
  if(length(i) > 0) frags <- frags[-i,]
  if(nrow(frags) == 0) return(data.frame())
  
  # Remove combinations which have the same strand since fragment reads are expected to have opposite strands.
  i <- which(frags$strand.anchorReads == frags$strand.adriftReads)
  if(length(i) > 0) frags <- frags[-i,]
  if(nrow(frags) == 0) return(data.frame())
  
  # Determine the start and end of fragments based on their alignment strands
  # and perform some sanity tests then filter on fragment size. 
  
  frags <- mutate(frags, 
                  fragStart  = ifelse(strand.anchorReads == '+', tStart.anchorReads + 1, tStart.adriftReads + 1),
                  fragEnd    = ifelse(strand.anchorReads == '+', tEnd.adriftReads + 1,   tEnd.anchorReads + 1),
                  strand     = ifelse(strand.anchorReads == '+', '+', '-'),
                  chromosome = tName.anchorReads,  
                  fragTest  = ifelse(strand.anchorReads == '+', tStart.anchorReads < tEnd.adriftReads, tStart.adriftReads < tEnd.anchorReads),  
                  fragWidth = (fragEnd - fragStart) + 1) %>%
           filter(fragTest == TRUE, 
                  fragWidth <= opt$buildFragments_maxFragLength,
                  fragWidth >= opt$buildFragments_minFragLength) %>%
                  mutate(uniqueSample = sample.anchorReads, readID = qName.anchorReads) %>%
                  select(uniqueSample, readID, chromosome, strand, fragStart, fragEnd, leaderSeq.anchorReads)
  
  # Reads pairs building more than fragment are multihits.
  # A single fragment can map to multiple reads but a single read can only map to one fragment.
  dupReads <- unique(frags$readID[duplicated(frags$readID)])
  
  
  if(length(dupReads) > 0) write(dupReads, file = file.path(opt$outputDir, opt$buildFragments_outputDir, 'multiHitReads', tmpFile()))
  
  frags <- bind_rows(lapply(split(frags, paste(frags$uniqueSample, frags$strand, frags$fragStart, frags$fragEnd)), function(a){
             a <- left_join(a, dups, by = c('readID' = 'id'))
                         
             # Remove multi-hit reads.
             a <- subset(a, ! readID %in% dupReads)
             if(nrow(a) == 0) return(tibble())
             
             a$reads <-sum(a$n, na.rm =TRUE) + nrow(a)
             r <- representativeSeq(a$leaderSeq.anchorReads)
             
             # Exclude reads where the leaderSeq is not similiar to the representative sequence. 
             i <- stringdist::stringdist(r[[2]], a$leaderSeq.anchorReads) / nchar(r[[2]]) <= opt$buildFragments_maxLeaderSeqDiffScore
             if(all(! i)) return(data.frame())
             
             a <- a[i,]
             a$repLeaderSeq <- r[[2]]
             a[1, c(1,3:7,9,10)]
           }))

  frags
}))


# Record multihit alignments for later analyses.
m <- unique(unlist(lapply(list.files(file.path(opt$outputDir, opt$buildFragments_outputDir, 'multiHitReads'), full.names = TRUE), readLines)))
readr::write_tsv(subset(anchorReadAlignments, qName.anchorReads %in% m), file.path(opt$outputDir, opt$buildFragments_outputDir, 'multiHitAnchorReadAlignments.tsv'))
readr::write_tsv(subset(adriftReadAlignments, qName.adriftReads %in% m), file.path(opt$outputDir, opt$buildFragments_outputDir, 'multiHitAdriftReadAlignments.tsv'))

# Clean up.
unlink(file.path(opt$outputDir, opt$buildFragments_outputDir, 'multiHitReads'), recursive = TRUE)

saveRDS(frags, file.path(opt$outputDir, opt$buildFragments_outputDir, opt$buildFragments_outputFile))

