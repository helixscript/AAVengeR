library(yaml)
library(dplyr)

opt <- read_yaml('config.yml')

anchorReadAlignments <- readRDS(file.path(opt$outputDir, 'alignedReads', 'anchorReadAlignments.rds'))
adriftReadAlignments <- readRDS(file.path(opt$outputDir, 'alignedReads', 'adriftReadAlignments.rds'))

names(anchorReadAlignments) <- paste0(names(anchorReadAlignments), '.anchorReads')
names(adriftReadAlignments) <- paste0(names(adriftReadAlignments), '.adriftReads')


# Here we split the alignments into groups of <= 1000 read ids to prevent the left_join from 
# exceeding R's internal table size limit.

ids <- unique(anchorReadAlignments$qName.anchorReads)
id_groups <- split(ids, dplyr::ntile(1:length(ids), round(length(ids)/1000)))

cluster <- makeCluster(opt$buildFragments_CPUs)

frags <- bind_rows(parLapply(cluster, id_groups, function(id_group){
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
 
  mutate(frags, 
         fragStart = ifelse(strand.anchorReads == '+', tStart.anchorReads + 1, tStart.adriftReads + 1),
         fragEnd   = ifelse(strand.anchorReads == '+', tEnd.adriftReads + 1,   tEnd.anchorReads + 1),
         fragTest  = ifelse(strand.anchorReads == '+', tStart.anchorReads < tEnd.adriftReads, tStart.adriftReads < tEnd.anchorReads),  
         fragWidth = (fragEnd - fragStart) + 1) %>%
    filter(fragTest == TRUE, 
           fragWidth <= opt$buildFragments_maxFragLength,
           fragWidth >= opt$buildFragments_minFragLength) %>%
    mutate(uniqueSample = sample.anchorReads, readID = qName.anchorReads) %>%
    select(-sample.anchorReads, -qName.anchorReads)
}))




