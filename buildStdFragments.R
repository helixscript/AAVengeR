library(dplyr)
library(lubridate)
library(parallel)
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(igraph)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$buildStdFragments_outputDir))

dir.create(file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'randomIDexcludedReads'))

write(c(paste(now(), 'Reading in fragment file(s).')), file = file.path(opt$outputDir, 'log'), append = TRUE)

# Load fragment data from database or local file depending on configuration file.
# Fragments are read in on the read level and contain a variable (n) which is the number 
# of identical read pairs that were removed in prepReads. (n) can be added to fragment
# read counts to account for all reads suporting fragments.


if('databaseGroup' %in% names(opt)){
  library(RMariaDB)
  library(DBI)
  
  samples <- loadSamples()
  
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


# Create fragment ids. (Very slow...)
frags <- unpackUniqueSampleID(frags)


# Look at ITR/LTR remnants on the subject level, order by UMI counts, and assign
# identifiers to each so they can be standardized separately.
frags <- bind_rows(lapply(split(frags, paste(frags$trial, frags$subject)), function(x){
  
  d <- group_by(x, leaderSeq.anchorReads) %>% 
       summarise(reads = n(), UMIs = n_distinct(randomLinkerSeq.adriftReads)) %>%
       arrange(desc(UMIs)) %>%
       mutate(leaderSeqGroup = NA)
  
  g <- 1
  
  invisible(lapply(1:nrow(d), function(i){
    if(! is.na(d[i,]$leaderSeqGroup)) return()
    
    maxEditDist <- round(nchar(as.character(d[i,]$leaderSeq.anchorReads))/opt$buildStdFragments_categorize_anchorReadRemnants_stepSize) + 1
    
    d[i,]$leaderSeqGroup <<- paste0(x$trial[1], '~', x$subject[1], '~', g)
    
    o <- d[is.na(d$leaderSeqGroup),]
    if(nrow(o) == 0) return()
    
    dist <- stringdist::stringdist(d[i,]$leaderSeq.anchorReads, o$leaderSeq.anchorReads)
    
    k <- which(dist <= maxEditDist)
    
    if(length(k) > 0){
      k2 <- which(d$leaderSeq.anchorReads %in% unique(o[k,]$leaderSeq.anchorReads))
      d[k2,]$leaderSeqGroup <<- paste0(x$trial[1], '~', x$subject[1], '~', g)
    }
    
    g <<- g + 1
  }))
  
  left_join(x, select(distinct(d), leaderSeq.anchorReads, leaderSeqGroup), by = 'leaderSeq.anchorReads')
}))


# Define fragment ids -- everything that makes a fragment unique.
frags$fragID <- paste0(frags$trial,     ':', frags$subject,    ':', frags$sample, ':',  
                       frags$replicate, ':', frags$chromosome, ':', frags$strand, ':', 
                       frags$fragStart, ':', frags$fragEnd,    ':', frags$leaderSeqGroup, ':',
                       frags$randomLinkerSeq.adriftReads)


# Create a tibble that can be turned into a GRange object which then can be used to standardize positions.
# frags is on the read level, here we create f which tallies read counts for each fragment and we 
# standardize within ITR/LTR groupings across subjects. This will prevent closely spaced events 
# with different ITR/LTR remnants from being merged. Updated positions can be joined and updated via fragID.

f <- group_by(frags, fragID) %>% 
     summarise(seqnames = chromosome[1], start = fragStart[1], end = fragEnd[1], strand = strand[1], 
               reads = n_distinct(readID) + sum(n), fragID = fragID[1], s = leaderSeqGroup[1]) %>%
     ungroup()


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


# Record original fragment positions.
frags$orgFragStart <- frags$fragStart
frags$orgFragEnd   <- frags$fragEnd

# Join the standardized start/stop positions to the fragments data frame.
frags <- left_join(frags, select(g, start, end, fragID), by = 'fragID')
rm(g)


# Assign the new integration positions.
frags$fragStart <- frags$start
frags$fragEnd   <- frags$end
frags <- select(frags, -start, -end)

# Create position ids including the leader sequence grouping ids.
frags$posid <- paste0(frags$chromosome, frags$strand, 
                      ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd),
                      '.', stringr::str_extract(frags$leaderSeqGroup, '\\d+$'))

# Create position ids including the leader sequence grouping ids.
frags$fragID <- paste0(frags$trial, ':', frags$subject, ':', frags$sample, ':', frags$replicate, ':',
                       frags$chromosome, ':', frags$strand, ':', frags$fragStart, ':', frags$fragEnd, ':', 
                       frags$leaderSeqGroup, ':', frags$randomLinkerSeq.adriftReads)


# Rather than standardizing across subjects, here we standardize breaks within replicate / random id groups
# and parallelize across unique samples.
f <- group_by(frags, fragID) %>%
     summarise(seqnames = chromosome[1], start = fragStart[1], end = fragEnd[1], strand = strand[1], 
               reads = n_distinct(readID) + sum(n), fragID = fragID[1],
               n = n(),
               s1 = uniqueSample[1],
               s2 = paste0(uniqueSample[1], ':', posid[1], ':', randomLinkerSeq.adriftReads[1])) 


# Separate fragment records into those with and those without multiple reads
# since we do not need to standardize fragments with a single read.

z <- table(f$s2)
f1 <- subset(f, ! s2 %in% names(z[z > 1]))
f2 <- subset(f, s2 %in% names(z[z > 1]))


if(nrow(f2) > 0){
  
  message(opt$buildStdFragments_CPUs, ' cores.')
  
  o <- split(f2, f2$s2)
  i <- ntile(1:length(o), opt$buildStdFragments_CPUs)
  f2 <- bind_rows(lapply(1:length(o), function(x){
        a <- o[[x]]
        a$s3 <- i[x]
        a
  }))
  
  f2 <- bind_rows(parLapply(cluster, split(f2, f2$s3), function(x){
          library(dplyr) 
          library(GenomicRanges)
          source(file.path(opt$softwareDir, 'lib.R'))
      
          g <- GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)

          g <- unlist(GenomicRanges::GRangesList(lapply(split(g, g$s2), function(x){
            
                 x$breakPointRefined <- FALSE
  
                 out <- tryCatch({
                           o <- refine_breakpoints(x, counts.col = 'reads')
                           o$breakPointRefined <- TRUE
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
  
          data.frame(g)  
        }))
}

stopCluster(cluster)

f <- bind_rows(f1, f2)
rm(f1, f2)

# Record original fragment positions.
frags$orgFragStart2 <- frags$fragStart
frags$orgFragEnd2   <- frags$fragEnd

# Join the standardized start/stop positions to the fragments data frame.
frags <- left_join(frags, select(f, start, end, fragID), by = 'fragID')

# Assign the new break positions.
frags$fragStart <- frags$start
frags$fragEnd   <- frags$end
frags <- select(frags, -start, -end)


# Update position and fragment IDs to reflect standardization.
frags$posid <- paste0(frags$chromosome, frags$strand, 
                      ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd),
                      '.', stringr::str_extract(frags$leaderSeqGroup, '\\d+$'))

frags$fragID <- paste0(frags$trial,          ':', frags$subject, ':', frags$sample,    ':', frags$replicate, ':',
                       frags$chromosome,     ':', frags$strand,  ':', frags$fragStart, ':', frags$fragEnd,   ':', 
                       frags$leaderSeqGroup, ':', frags$randomLinkerSeq.adriftReads)


# Create an additional fragment ID without the random seq id.
frags$fragID2 <- paste0(frags$trial,     ':', frags$subject,    ':', frags$sample, ':',  
                        frags$replicate, ':', frags$chromosome, ':', frags$strand, ':', 
                        frags$fragStart, ':', frags$fragEnd,    ':', frags$leaderSeqGroup)


# Correct for instances where a read maps to more than fragment but all fragments have the sample integration position.
# These are instances of fuzzy break points and here we select the shortest fragments lengths.

message('f9')

o <- split(frags, frags$readID)
frags <- bind_rows(lapply(o, function(x){
  if(n_distinct(x$fragID2) > 1 & n_distinct(x$posid) == 1){
    if(x$strand[1] == '+'){
      i <- which(x$fragEnd == min(x$fragEnd))[1]
      x <- x[i,]
    } else {
      i <- which(x$fragStart == max(x$fragStart))[1]
      x <- x[i,]
    }
  }
  x
}))


# Identify reads which map to more than position id and define these as multi-hit reads.
o <- group_by(frags, readID) %>% 
     summarise(n = n_distinct(fragID2)) %>%
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
  
  # The posid leader sequence identifiers are calculated on the subject level 
  # since posid ids are standardized across patients not samples.
  
  multiHitFrags$s <- paste(multiHitFrags$trial, multiHitFrags$subject)
  
  # Here we cycle through each sample or subject, create a list of uniquely resolved sites
  # and flag reads to be returned to frags if only one site in the list of potential sites is in the list uniquely resolved sites.
  
  multiHitFrags$returnToFrags <- FALSE
  
  multiHitFrags <- bind_rows(lapply(split(multiHitFrags, multiHitFrags$s), function(x){
                     uniqueSites <- unique(dplyr::filter(frags, trial == x$trial[1], subject == x$subject[1])$posid)
    
                     bind_rows(lapply(split(x, x$readID), function(x2){
                       if(sum(unique(x2$posid) %in% uniqueSites) == 1) x2[x2$posid %in% uniqueSites,]$returnToFrags <- TRUE
                       x2
                     }))
                   }))
  
  # Return reads to frags if needed.
  
  if(any(multiHitFrags$returnToFrags == TRUE)){
    # Isolate frag reads that should be returned.
    m <- subset(multiHitFrags, returnToFrags == TRUE)
    
    # We correct for instances where a read supports for multiple breaks but a single position (again)
    # because reads may of been missed earlier since reads may of muspported multiple positions in the first pass.
    m <- bind_rows(lapply(split(m, m$readID), function(x){
      if(n_distinct(x$fragID2) > 1 & n_distinct(x$posid) == 1){
        if(x$strand[1] == '+'){
          i <- which(x$fragEnd == min(x$fragEnd))[1]
          x <- x[i,]
        } else {
          i <- which(x$fragStart == max(x$fragStart))[1]
          x <- x[i,]
        }
      }
      x
    }))
    
    write(paste0(sprintf("%.2f%%", (nrow(m) / nrow(multiHitFrags))*100), ' of multihit reads recovered because one potential site was in the unabiguous site list.'), file = file.path(opt$outputDir, 'log'), append = TRUE)
    
    # Remove salvaged reads from multi-hit reads.
    multiHitFrags <- subset(multiHitFrags, ! readID %in% m$readID) %>% dplyr::select(-s, -returnToFrags)
    
    # Return salvaged reads to frags.
    frags <- bind_rows(frags, dplyr::select(m, -s, -returnToFrags))
  }
} 

multiHitClusters <- tibble()

save.image(file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitClustersStart.RData'))

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

# Frags are still read level. frags is set classes as a data frame because
# tibbles refuse to store single character vectors as lists needed for read ID lists
# in fragment records.
frags <- data.frame(frags)
fragsRemoved <- tibble()


# Here we split the read level data by boundary corrected fragment ids, 
# call representativeSeq() to perform an alignment of the ITR/LTR remnant sequences
# to find the representative leader sequence and determine fragment read counts. 
# The reads per fragment filter is applied here since it will effect the following 
# clean up of UMIs found across multiple fragment records.

o <- split(frags, frags$fragID)

t <- length(o)
i <- 1

f <- bind_rows(lapply(o, function(x){
       message(i, '/', t); i <<- i + 1
       
       totalReads <- n_distinct(x$readID) + sum(x$n)
  
       if(totalReads < opt$buildStdFragments_minReadsPerFrag) return(tibble())
  
       r <- representativeSeq(x$leaderSeq.anchorReads)
  
       # Prevent repetitive calls to muscle in representativeSeq() from causing 
       # a system level error with too many open connections.
       if(i %% 100 == 0) closeAllConnections()
       
       x$repLeaderSeq <- r[[2]]
  
       leaderSeqs <- x$repLeaderSeq
       readList <- x$readID
  
       x <- x[1,]
       x$reads <- totalReads
       x$maxLeaderSeqDist <- NA # Too slow. Uses all cores. max(stringdist::stringdistmatrix(leaderSeqs))
       x$readIDs <- list(readList)
       x$leaderSeqs <- list(leaderSeqs)
  
       x
}))


# Clear out the tmp/ directory.
invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

# Check to see if UMIs are duplicated within samples.
# We previously checked to make sure they are not spread across samples.
#
# If they are then select the fragment record with the most reads.
# For each fragment retained, we record the percentage of total reads the   
# selected fragment with respect to all reads within the sample with the same random id.

f$randomLinkerSeq.dispute <- FALSE
f$randomLinkerSeq.selectedReadMajority <- NA


# Within each sample, look for instances of duplicate UMIs which typically result
# from PCR cross over events and select the UMI based fragment with the most reads. 
# We record the percentage iof UMI reads attributed to the most read fragment for 
# down stream filtering.

# o$fragWidth <- o$fragEnd - o$fragStart + 1
# select(o, randomLinkerSeq.adriftReads, uniqueSample, posid, fragStart, fragEnd, reads, repLeaderSeq) 

f2 <- bind_rows(lapply(split(f, paste(f$trial, f$subject, f$sample)), function(x){
       t <- table(x$randomLinkerSeq.adriftReads) 
       t <- names(t[which(t > 1)])
       
       if(length(t) > 0){
         a <- subset(x, ! randomLinkerSeq.adriftReads %in% t)
         b <- subset(x, randomLinkerSeq.adriftReads %in% t)
         
         b$randomLinkerSeq.dispute <- TRUE
    
         b2 <- bind_rows(lapply(t, function(i){
                 o <- dplyr::arrange(subset(b, randomLinkerSeq.adriftReads == i), desc(reads))
                 o$percentReads <- (o$reads / sum(o$reads))*100

                 if(o[1,]$percentReads < opt$buildStdFragments_UMI_conflict_minPercentReads){
                   o$msg <- 'All fragments rejected'
                   o$percentReads <- sprintf("%.2f%%", o$percentReads)
                   readr::write_tsv(dplyr::select(o, trial, subject, sample, replicate, readID, randomLinkerSeq.adriftReads, posid, percentReads, msg), 
                                    file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'randomIDexcludedReads', paste0('duplicate_UMIs_withinSample~', x$trial[1], '~', x$subject[1], '~', x$sample[1], '.tsv')))
                   return(tibble())
                 } else {
                   o$msg <- 'Top fragment retained'
                   o$percentReads <- sprintf("%.2f%%", o$percentReads)
                   readr::write_tsv(dplyr::select(o[2:nrow(o),], trial, subject, sample, replicate, readID, randomLinkerSeq.adriftReads, posid, percentReads, msg), 
                                    file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'randomIDexcludedReads', paste0('duplicate_UMIs_withinSample~', x$trial[1], '~', x$subject[1], '~', x$sample[1], '.tsv')))
                 }
                 o[1,]
               }))
       
         x <- bind_rows(a, b2)
       }
    
       x
     }))

f2 <- select(f2, -uniqueSample, -n, -fragID, -fragID2, -readID, -leaderSeq.anchorReads)
saveRDS(f2, file.path(opt$outputDir, opt$buildStdFragments_outputDir, opt$buildStdFragments_outputFile))

q(save = 'no', status = 0, runLast = FALSE) 
