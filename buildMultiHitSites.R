library(dplyr)
library(parallel)
library(igraph)

opt <- yaml::read_yaml('config.yml')
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$buildMultiHitSites_outputDir))

if(! file.exists(file.path(opt$outputDir, opt$buildMultiHitSites_inputReadsFile))) stop('Error - input reads data object not found.')
if(! file.exists(file.path(opt$outputDir, opt$buildMultiHitSites_inputFragsFile))) stop('Error - input frags data object not found.')

samples <- loadSamples()

multiHitFrags <- readRDS(file.path(opt$outputDir, opt$buildMultiHitSites_inputFragsFile))
multiHitReads <- readRDS(file.path(opt$outputDir, opt$buildMultiHitSites_inputReadsFile))

multiHitFrags$width <- multiHitFrags$fragEnd - multiHitFrags$fragStart + 1

multiHitSites <- group_by(multiHitFrags, subject, sample, posid) %>%
                 summarise(reads = sum(reads), estAbund = n_distinct(width)) %>%
                 ungroup()

multiHitSites <- bind_rows(lapply(split(multiHitSites, paste(multiHitSites$subject, multiHitSites$sample)), function(x){
  nodes <- data.frame(name = x$posid, estAbund = x$estAbund, reads = x$reads)
  r <- subset(multiHitReads, subject == x$subject[1] & sample == x$sample[1] & posids %in% nodes$name)
  
  edges <- bind_rows(lapply(split(r, r$readID), function(y){
    o <- RcppAlgos::comboGeneral(unique(y$posids), 2)
    data.frame(from = o[,1], to = o[,2])
  }))
  
  g <- igraph::simplify(graph_from_data_frame(edges, directed=FALSE, vertices=nodes))
  o <- tibble(subject = x$subject[1], sample = x$sample[1], clusters = lapply(igraph::decompose(g), function(x) igraph::V(x)$name))
  
  s <- bind_rows(lapply(split(o, 1:nrow(o)), function(y){
         a <- subset(nodes, name %in% unlist(y$clusters))
         y$reads <- sum(a$reads)
         y$estAbund <- max(a$estAbund)
         y
       }))
  
  s$clusterNum <- 1:nrow(s)
  s$numSites <- unlist(lapply(s$clusters, function(o) length(unlist(o))))
  s$posid <- paste0('chrMulti*', s$clusterNum, ',', s$numSites)
  s$chromosome <- 'chrMulti'
  s$strand <- '*'
  s$position <- NA
  s$fragmentsRemoved <- NA
  s$repLeaderSeq <- NA
  s$flags <- paste0(unique(subset(samples, subject == x$subject[1] & sample == x$sample[1])$flags), collapse = ',')
  dplyr::select(s, subject, sample, chromosome, strand, position, posid, estAbund, reads, fragmentsRemoved, repLeaderSeq, flags)
}))
