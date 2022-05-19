library(dplyr)
library(parallel)
library(igraph)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')

opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$buildMultiHitSites_outputDir))

if(! file.exists(file.path(opt$outputDir, opt$buildMultiHitSites_inputReadsFile))) stop('Error - input reads data object not found.')
if(! file.exists(file.path(opt$outputDir, opt$buildMultiHitSites_inputFragsFile))) stop('Error - input frags data object not found.')

samples <- loadSamples()

multiHitFrags <- readRDS(file.path(opt$outputDir, opt$buildMultiHitSites_inputFragsFile))
multiHitReads <- readRDS(file.path(opt$outputDir, opt$buildMultiHitSites_inputReadsFile))

if(nrow(multiHitFrags) > 0 & nrow(multiHitReads) > 0){
  multiHitReads <- unpackUniqueSampleID(tidyr::unnest(multiHitReads, posids))

  multiHitFrags$width <- multiHitFrags$fragEnd - multiHitFrags$fragStart + 1
  multiHitFrags$posid <- paste0(multiHitFrags$chromosome, multiHitFrags$strand, ifelse(multiHitFrags$strand == '+', multiHitFrags$fragStart, multiHitFrags$fragEnd))

  multiHitSites <- group_by(multiHitFrags, trial, subject, sample, posid) %>%
                   summarise(reads = sum(reads), estAbund = n_distinct(width)) %>%
                   ungroup()

  multiHitSites <- bind_rows(lapply(split(multiHitSites, paste(multiHitSites$subject, multiHitSites$sample)), function(x){

  # Define posid network nodes.
  nodes <- data.frame(name = x$posid, estAbund = x$estAbund, reads = x$reads)
  
  # Retrieve reads IDs for posid nodes. 
  r <- subset(multiHitReads, subject == x$subject[1] & sample == x$sample[1] & posids %in% nodes$name)
    # For each read id, create different permutations of posids connected by the read id.
    edges <- bind_rows(lapply(split(r, r$readIDlist), function(y){
      o <- RcppAlgos::comboGeneral(unique(y$posids), 2)
      data.frame(from = o[,1], to = o[,2])
    }))
    
    # Build simple graph with no loops or multiple edges.
    g <- igraph::simplify(graph_from_data_frame(edges, directed=FALSE, vertices=nodes))
  
    # Separate out individual graphs.
    o <- tibble(trial = x$trial[1], subject = x$subject[1], sample = x$sample[1], clusters = lapply(igraph::decompose(g), function(x) igraph::V(x)$name))
  
    s <- bind_rows(lapply(split(o, 1:nrow(o)), function(y){
           a <- subset(nodes, name %in% unlist(y$clusters))
           y$reads <- sum(a$reads)
           y$estAbund <- max(a$estAbund)
           y
         }))
  
    s$clusterNum <- 1:nrow(s)
    s$numSites <- unlist(lapply(s$clusters, function(o) length(unlist(o))))
    s$posid <- paste0('clust', s$clusterNum, '.', s$numSites)
    s$chromosome <- 'Mult'
    s$strand <- '*'
    s$position <- NA
    s$fragmentsRemoved <- NA
    s$repLeaderSeq <- NA
    s$flags <- paste0(unique(subset(samples, subject == x$subject[1] & sample == x$sample[1])$flags), collapse = ',')
    s$replicate <- NA
    dplyr::select(s, trial, subject, sample, replicate, chromosome, strand, position, posid, estAbund, reads, fragmentsRemoved, repLeaderSeq, flags)
  }))
} else {
  multiHitSites <- tibble()
}

saveRDS(multiHitSites, file.path(opt$outputDir, opt$buildMultiHitSites_outputDir,  opt$buildMultiHitSites_outputFile))

q(save = 'no', status = 0, runLast = FALSE) 

