library(dplyr)
library(lubridate)
library(parallel)
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(igraph)
library(multidplyr)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$buildStdFragments_outputDir))

write(c(paste(now(), 'Reading in fragment file(s).')), file = file.path(opt$outputDir, 'log'), append = TRUE)

samples <- loadSamples()



if('databaseGroup' %in% names(opt)){
  library(RMariaDB)
  library(DBI)
  con <- dbConnect(RMariaDB::MariaDB(), group = opt$databaseGroup)
  
  write(c(paste(now(), '  Reading in fragment data from database.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
  # Retrieve all fragments for unique trial / subject pairs.
  # This may retrieve fragments for samples outside of the current sample list but 
  # this is intended since we need to standardize sites across all samples.
  o <- distinct(dplyr::select(samples, trial, subject))
  frags <- bind_rows(lapply(split(o, 1:nrow(o)), function(x){
    r <- dbGetQuery(con, paste0("select * from fragments where subject='", x$subject, "' and trial='", x$trial, "'"))
    r2 <- tibble()
    if(nrow(r) > 0){
      r2 <- bind_rows(lapply(split(r, 1:nrow(r)), function(x2){
         d <- unserialize(x2$data[[1]])
         d$refGenome <- x2$refGenome
         d
      }))
    }
    r2
  }))
  
  dbDisconnect(con)
} else {
  opt$buildStdFragments_inputFiles <- gsub('\\s*,\\s*', ',', opt$buildStdFragments_inputFiles)

  frags <- bind_rows(lapply(unlist(strsplit(opt$buildStdFragments_inputFiles, ',')), function(x){
             write(c(paste(now(), '  reading in local data file ', x)), file = file.path(opt$outputDir, 'log'), append = TRUE)
             readRDS(file.path(opt$outputDir, x))
           }))
}

# Create fragment ids.
frags <- unpackUniqueSampleID(frags)
frags$fragID <- paste0(frags$trial, ':', frags$subject, ':', frags$sample, ':', frags$replicate, ':', frags$chromosome, ':', 
                       frags$strand, ':', frags$fragStart, ':', frags$fragEnd, ':', frags$randomLinkerSeq.adriftReads)


# Create a tibble that can be turned into a GRange object which then can be used to standardize positions.
# Updated positions can be joined and updated via fragID.


cl <- new_cluster(30)
cluster_library(cl, "dplyr")
f <- frags %>% group_by(fragID) %>% partition(cl) %>% 
     summarise(seqnames = chromosome[1], start = fragStart[1], end = fragEnd[1], strand = strand[1], 
            reads = n_distinct(readID) + sum(n), fragID = fragID[1], s = paste(trial[1], subject[1])) %>%
     collect()


cluster <- makeCluster(opt$buildStdFragments_CPUs)
clusterExport(cluster, 'opt')

# Standardize integration positions.
g <- GenomicRanges::makeGRangesFromDataFrame(f, keep.extra.columns = TRUE)
g <- unlist(GenomicRanges::GRangesList(parallel::parLapply(cluster, split(g, g$s), function(x){
       library(dplyr)
       library(GenomicRanges)
       source(file.path(opt$softwareDir, 'lib.R'))
  
       x$intSiteRefined <- FALSE
  
       out <- tryCatch({
         o <- standardize_sites(x, counts.col = 'reads', sata.gap = 5)
         o$intSiteRefined <- TRUE
         o
       },
       error=function(cond) {
         x
       },
       warning=function(cond) {
         o
       })
  
      return(out)
})))

g <- data.frame(g)

frags$orgFragStart <- frags$fragStart
frags$orgFragEnd <- frags$fragEnd

frags <- left_join(frags, select(g, start, end, fragID), by = 'fragID')
rm(g)

# Assign the new integration positions.
frags$fragStart <- frags$start;  frags$start <- NULL
frags$fragEnd <- frags$end;      frags$end <- NULL
 


frags$posid <- paste0(frags$chromosome, frags$strand, ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd))
frags$fragID2 <- paste0(frags$trial, ':', frags$subject, ':', frags$sample, ':', frags$replicate, ':', frags$chromosome, ':', 
                       frags$strand, ':', frags$fragStart, ':', frags$fragEnd)




# Correct for instances where a read maps to more than fragment but all fragments have the sample integration position.
# These are instances of fuzzy break points and here we select the shortest fragments.

cat(NULL, file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'adjustedVariableBreakFragReads'), append = FALSE)
o <- split(frags, frags$readID)
frags <- bind_rows(lapply(o, function(x){
  if(n_distinct(x$fragID2) > 1 & n_distinct(x$posid) == 1){
    if(x$strand[1] == '+'){
      i <- which(x$fragEnd == min(x$fragEnd))[1]
      write(x[-i,]$fragID2, file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'adjustedVariableBreakFragReads'), append = TRUE)
      x <- x[i,]
    } else {
      i <- which(x$fragStart == max(x$fragStart))[1]
      write(x[-i,]$fragID2, file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'adjustedVariableBreakFragReads'), append = TRUE)
      x <- x[i,]
    }
  }
  x
}))
o <- unique(readLines(file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'adjustedVariableBreakFragReads')))
write(o, file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'adjustedVariableBreakFragReads'), append = FALSE)



# Identify reads which map to more than position id and define these as multi-hit reads.
o <- group_by(frags, readID) %>% partition(cl) %>%
     summarise(n = n_distinct(fragID2)) %>%
     collect() %>%
     dplyr::filter(n > 1)

multiHitFrags <- tibble()

if(nrow(o) > 0){
  
  # If there are reads associated with more than one fragment, we look at the multiple positions
  # associated with each read. If one, and only one of those positions were uniquely resolved then 
  # read is returned to frags.
  multiHitFrags <- subset(frags, readID %in% o$readID)  # Reads which map to more than one fragment.
  frags <- subset(frags, ! readID %in% o$readID)        # Reads which map to a single intSite positions.
  
  # There may be multiHit reads where one of their position ids is the same as those in the frags
  # where the reads mapped uniquely. If only one of the possible alignment positions in multiHitFrags 
  # is in frags (unique alignments), then move those reads to frags.
  
  if(opt$buildStdFragments_salvageMultiHits_within == 'sample'){
    multiHitFrags$s <- paste(multiHitFrags$trial, multiHitFrags$subject, multiHitFrags$sample)
  } else {
    multiHitFrags$s <- paste(multiHitFrags$trial, multiHitFrags$subject)
  }
  
  
  
  # Here we cycle through each sample or subject, create a list of uniquely resolved sites
  # and flag reads to be returned to frags if only one site in the list of potential sites is in the list uniquely resolved sites.
  
  multiHitFrags$returnToFrags <- FALSE
  
  multiHitFrags <- bind_rows(lapply(split(multiHitFrags, multiHitFrags$s), function(x){
                     uniqueSites <- unique(dplyr::filter(frags, trial == x$trial[1], subject == x$subject[1], sample %in% x$sample)$posid)
    
                     bind_rows(lapply(split(x, x$readID), function(x2){
                       if(sum(unique(x2$posid) %in% uniqueSites) == 1) x2[x2$posid %in% uniqueSites,]$returnToFrags <- TRUE
                       x2
                     }))
                   }))
  
  # Return reads to frags if needed.
  
  if(any(multiHitFrags$returnToFrags == TRUE)){
    # Isolate frag reads that should be returned.
    m <- subset(multiHitFrags, returnToFrags == TRUE)
    
    # We correct for instances where a read supports for multiple breaks but a single position
    # because this may of been missed earlier since reads may of muspported multiple positions in the first pass.
    m <- bind_rows(lapply(split(m, m$readID), function(x){
      if(n_distinct(x$fragID2) > 1 & n_distinct(x$posid) == 1){
        if(x$strand[1] == '+'){
          i <- which(x$fragEnd == min(x$fragEnd))[1]
          write(x[-i,]$fragID2, file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'adjustedVariableBreakFragReads'), append = TRUE)
          x <- x[i,]
        } else {
          i <- which(x$fragStart == max(x$fragStart))[1]
          write(x[-i,]$fragID2, file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'adjustedVariableBreakFragReads'), append = TRUE)
          x <- x[i,]
        }
      }
      x
    }))
    o <- unique(readLines(file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'adjustedVariableBreakFragReads')))
    write(o, file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'adjustedVariableBreakFragReads'), append = FALSE)
    
    write(paste0(sprintf("%.2f%%", (nrow(m) / nrow(multiHitFrags))*100), ' of multihit reads recovered because one potential site was in the unabiguous site list.'), file = file.path(opt$outputDir, 'log'), append = TRUE)
    
    # Remove salvaged reads from multi-hit reads.
    multiHitFrags <- subset(multiHitFrags, ! readID %in% m$readID) %>% dplyr::select(-s, -returnToFrags)
    
    # Return salvaged reads to frags.
    frags <- bind_rows(frags, dplyr::select(m, -s, -returnToFrags))
  }
} 

multiHitClusters <- tibble()

if(nrow(multiHitFrags) > 0 & opt$buildStdFragments_createMultiHitClusters){
  
  multiHitFrags$fragWidth = multiHitFrags$fragEnd - multiHitFrags$fragStart + 1
  
  # For each read, create a from -> to data frame and capture the width of the read. 
  multiHitNet_replicates <- bind_rows(lapply(split(multiHitFrags, multiHitFrags$readID), function(x){
    if(n_distinct(x$posid) == 1) return(tibble()) # Cases of break point only variation.
    node_pairs <- RcppAlgos::comboGeneral(unique(x$posid), 2)
    tibble(trial = x[1,]$trial, subject = x[1,]$subject, sample = x[1,]$sample, 
           replicate = x[1,]$replicate, from = node_pairs[,1], to = node_pairs[,2], width = x[1,]$fragWidth, readID = x[1,]$readID)
  }))
  
  # Build networks for each sample.
  multiHitClusters <- bind_rows(lapply(split(multiHitNet_replicates, paste(multiHitNet_replicates$trial, multiHitNet_replicates$subject, multiHitNet_replicates$sample)), function(x){
    
    # Create a local copy of read-level multiHitFrags data frame specific to this sample.
    multiHitFrags <- subset(multiHitFrags, trial = x$trial[1], subject = x$trial[1], sample = x$sample[1])
    
    # Build a graph with the posid nodes and read edges.
    g <- igraph::simplify(graph_from_data_frame(dplyr::select(x, from, to), directed=FALSE, vertices=tibble(name = unique(c(x$from, x$to)))))
    
    # Separate out individual graphs.
    o <- tibble(trial = x$trial[1], subject = x$subject[1], sample = x$sample[1], clusters = lapply(igraph::decompose(g), function(x) igraph::V(x)$name))
    
    # Create cluster ids.
    o$clusterID <- paste0('MHC.', 1:nrow(o))
    
    bind_rows(lapply(split(o, o$clusterID), function(y){
      
      # Subset the larger multiHitFrag read data frame to focus on nodes in this graph.
      a <- subset(x, from %in% unlist(y$clusters) | to %in% unlist(y$clusters))
      
      # Determine node specific values.
      b <- bind_rows(lapply(unique(c(a$to, a$from)), function(b){
        a2 <- subset(a, to == b | from == b)
        tibble(node = b[1], reads = n_distinct(a2$readID), breaks = n_distinct(a2$width))
      }))
      
      y$nodes <- n_distinct(unlist(y$clusters))
      y$reads <- n_distinct(a$readID)
      y$abund <- n_distinct(a$width)
      y$maxNodeReads <- max(b$reads)
      y$minNodeReads <- min(b$reads)
      y$avgNodeReads <- mean(b$reads)
      y$nodeMostReads <- paste0(b[b$reads == max(b$reads),]$node, collapse = ',')
      y$nodeMostBreaks <- paste0(b[b$breaks == max(b$breaks),]$node, collapse = ',')
      y
    }))
  }))
  
  saveRDS(multiHitClusters, file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitClusters.rds'))
}


# Frags are still read level.
frags <- as_tibble(frags)
fragsRemoved <- tibble()

# Now we examine each integration and compare leader sequences to the most common at each position
# Then merge reads into fragments with using fragStart, fragEnd, and random id.

# (?) Here we are looking at all fragment reads assoicated with a posid and determining the representative 
# leader sequence. It would be better search for the ... chr16+59923338

o <- split(frags, paste(frags$trial, frags$subject, frags$sample, frags$replicate, frags$posid, 
                        frags$fragStart, frags$fragEnd, frags$randomLinkerSeq.adriftReads))

f <- bind_rows(lapply(o, function(x){
  
  r <- representativeSeq(x$leaderSeq.anchorReads)
  
  # Exclude reads where the leaderSeq is not similar to the representative sequence.
  i <- stringdist::stringdist(r[[2]], x$leaderSeq.anchorReads) / nchar(r[[2]]) <= opt$buildStdFragments_maxLeaderSeqDiffScore
  if(all(! i)) return(dplyr::tibble())
  
  # if(sum(!i) > 0){
  #   u <- tibble(trial = x$trial[1], subject = x$subject[1], sample = x$sample[1], replicate = x$replicate[1], 
  #               posid = x$posid[1], fragmentsRemoved = sum(!i), percentFragments = (sum(!i)/nrow(x))*100)
  #   fragsRemoved <- bind_rows(fragsRemoved, u)
  # }
  
  x <- x[i,]
  
  # Set the repLeaderSeq for this event.
  x$repLeaderSeq <- r[[2]]
  
  readList <- x$readID
  x$reads <- n_distinct(readList) + sum(x$n)
  x <- x[1,]
  
  if(x$reads < opt$buildStdFragments_minReadsPerFrag) return(tibble())
  
  x$readIDs <- as.list(list(readList))
  x
}))

# Clear out the tmp/ directory.
invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

# Check to see if UMIs are duplicated within samples.
# If they are then select the fragment record with the most reads.
# This will clean up some chatter around the non-standardized break points. 
# Duplicate UMIs across samples were removed earlier.

f <- bind_rows(lapply(split(f, paste(f$trial, f$subject, f$sample)), function(x){
       t <- table(x$randomLinkerSeq.adriftReads) 
       t <- names(t[which(t > 1)])
     
       if(length(t) >= 0){
         a <- subset(x, ! randomLinkerSeq.adriftReads %in% t)
         b <- subset(x, randomLinkerSeq.adriftReads %in% t)
      
         b2 <- bind_rows(lapply(t, function(i){
                 o <- dplyr::arrange(subset(b, randomLinkerSeq.adriftReads == i), desc(reads))
                 o[1,]
               }))
       
         x <- bind_rows(a, b2)
       }
    
       x
     }))


f <- select(f, -uniqueSample, -n, -fragID, -fragID2)
saveRDS(f, file.path(opt$outputDir, opt$buildStdFragments_outputDir, opt$buildStdFragments_outputFile))


q(save = 'no', status = 0, runLast = FALSE) 
