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
setMissingOptions()

dir.create(file.path(opt$outputDir, opt$buildStdFragments_outputDir))
dir.create(file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'tmp'))

# Load fragment data from database or local file depending on configuration file.
# Fragments are read in on the read level and contain a variable (n) which is the number 
# of identical read pairs that were removed in prepReads. (n) can be added to fragment
# read counts to account for all reads supporting fragments.

write(c(paste(now(), '   Reading in fragment file(s).')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = FALSE)

# Create a database connection if requested in the configuration file.
if('databaseGroup' %in% names(opt)){
  library(RMariaDB)

  dbConn <- tryCatch({
    dbConnect(RMariaDB::MariaDB(), group = opt$databaseGroup)
  },
  error=function(cond) {
    write(c(paste(now(), '   Error - could not connect to the database.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE) 
  })
}

# Check for configuration errors.
# ------------------------------------------------------------------------------
if('buildStdFragments_autoPullTrialSamples' %in% names(opt) & ! 'databaseGroup' %in% names(opt)){
  if(opt$buildStdFragments_autoPullTrialSamples){
    write(c(paste(now(), 'Error -- the databaseGroup option must be provided with the buildStdFragments_autoPullTrialSamples option.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE)
  }
}

if('buildStdFragments_trialSubjectList' %in% names(opt) & ! 'databaseGroup' %in% names(opt)){
  write(c(paste(now(), 'Error -- the databaseGroup option must be provided with the buildStdFragments_trialSubjectList option.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE)
}



# Read in fragment data.
# ------------------------------------------------------------------------------

if('buildStdFragments_inputFile' %in% names(opt)){
  # There is an incoming fragment rds file and here we make a table of unique trial / subject pairings that 
  # will be used to pull in other samples from the database. 
  frags <- distinct(readRDS(file.path(opt$outputDir, opt$buildStdFragments_inputFile)))
  
  if(opt$buildStdFragments_autoPullTrialSamples){
    trialSubjects <- unpackUniqueSampleID(tibble(uniqueSample = unique(frags$uniqueSample))) %>% dplyr::select(-uniqueSample, -replicate, -sample) %>% distinct()

    dbFrags <- bind_rows(lapply(split(trialSubjects, 1:nrow(trialSubjects)), function(x) pullDBsubjectFrags(dbConn, x$trial, x$subject)))
    
    dbFragsToAdd <- dbFrags[! dbFrags$uniqueSample %in% frags$uniqueSample,]
    
    if(nrow(dbFragsToAdd) > 0){
      write(c(paste(now(), paste0('   Adding ', nrow(dbFragsToAdd), ' fragment read records to incoming fragment data.'))), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)
      frags <- distinct(bind_rows(frags, dbFragsToAdd))
    }
  }
} else if('buildStdFragments_trialSubjectList' %in% names(opt)){
  # Format:  trial;subject|trial;subject, eg. 'Sabatino;pM50|Sabatino;pLinus'
  
  frags <- distinct(bind_rows(lapply(unlist(base::strsplit(opt$buildStdFragments_trialSubjectList, '\\|')), function(x){
             d <- unlist(base::strsplit(x, ';'))
             pullDBsubjectFrags(dbConn, d[1], d[2])
           })))
} else {
  write(c(paste(now(), 'Error -- neither buildStdFragments_inputFile or buildStdFragments_trialSubjectList options were provided')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE)
}

if('databaseGroup' %in% names(opt)) RMariaDB::dbDisconnect(dbConn)

# Make sure fragments were retrieved.
if(nrow(frags) == 0){
  write(c(paste(now(), 'Error -- no fragments were loaded or retrieved.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE)
}

incomingSamples <- unique(sub('~\\d+$', '', frags$uniqueSample))

# Set aside meta data till the end to save memory and add data back at the end 
# and unpack the uniqueSample ids.
# --------------------------------------------------------------------------------------------------------------
sampleMetaData <- distinct(dplyr::select(frags, uniqueSample, refGenome, vectorFastaFile, seqRunID, flags))

frags <- rbindlist(lapply(split(frags, frags$uniqueSample), function(x){
            o <-  unlist(stringr::str_split(x$uniqueSample[1], '~'))
            x$trial = o[1]; x$subject = o[2]; x$sample = o[3]; x$replicate = o[4]
            x
         })) %>% dplyr::select(-refGenome, -vectorFastaFile, -seqRunID, -flags)
   


# Categorize ITR/LTR remnant sequences.
# -------------------------------------------------------------------------------
write(c(paste(now(), '   Categorizing leader sequences.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)

frags$posid <- paste0(frags$chromosome, frags$strand, ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd))


frags$fragID <- paste0(frags$trial,     ':', frags$subject,    ':', frags$sample, ':',  
                       frags$replicate, ':', frags$chromosome, ':', frags$strand, ':', 
                       frags$fragStart, ':', frags$fragEnd,    ':', frags$randomLinkerSeq)


cluster <- makeCluster(opt$buildStdFragments_CPUs)
clusterExport(cluster, 'opt')

frags <- bind_rows(lapply(split(frags, paste(frags$trial, frags$subject)), function(x){

           bind_rows(parLapply(cluster, split(x, x$chromosome), function(x2){
             library(dplyr)
             processed <- vector()
             
             bind_rows(lapply(split(x2, x2$posid), function(x3){      
                
               # Find neighboring site fragments.
               if(x3$strand[1] == '+'){
                 o <- subset(x2, chromosome == x3$chromosome[1] & 
                                 strand == '+' & 
                                 fragStart >= x3$fragStart[1] - opt$buildStdFragments_leaderSeqGroupingDist & 
                                 fragStart <= x3$fragStart[1] + opt$buildStdFragments_leaderSeqGroupingDist)
               } else {
                 o <- subset(x2, chromosome == x2$chromosome[1] & 
                                 strand == '-' & 
                                 fragEnd >= x3$fragEnd[1] - opt$buildStdFragments_leaderSeqGroupingDist & 
                                 fragEnd <= x3$fragEnd[1] + opt$buildStdFragments_leaderSeqGroupingDist)
               }
           
               o <- filter(o, ! fragID %in% processed)
               if(nrow(o) == 0) return(tibble())
           
               processed <<- c(processed, unique(o$fragID))
           
               if(n_distinct(o$posid) > 1){
             
                 # There are multiple non-standardized position ids within the current one.
                 # Order their leader sequences by abundance and reads.
      
                 d <- mutate(o, widths = abs(fragEnd - fragStart)) %>%
                      group_by(leaderSeq) %>% 
                      summarise(reads = n_distinct(readID), estAbund = n_distinct(widths)) %>%
                      arrange(desc(estAbund), desc(reads)) 
             
                 d$leaderSeqGroup <- 'none'
                 g <- 1

                 invisible(lapply(1:nrow(d), function(i){
                   if(! d[i,]$leaderSeqGroup == 'none') return()

                   maxEditDist <- ceiling(nchar(as.character(d[i,]$leaderSeq))/opt$buildStdFragments_categorize_anchorReadRemnants_stepSize)

                   d[i,]$leaderSeqGroup <<- paste0(x$trial[1], '~', x$subject[1], '~', g)

                   # Retrieve all other leader sequences for which a leader seq group has not been assigned.
                   otherFrags <- d[d$leaderSeqGroup == 'none',]
                   if(nrow(otherFrags) == 0) return()

                   # Calculate the edit distances between current leader seq and other leader seqs.
                   dist <- stringdist::stringdist(d[i,]$leaderSeq, otherFrags$leaderSeq, nthread = 2)

                   k <- which(dist <= maxEditDist)

                   if(length(k) > 0){
                     k2 <- which(d$leaderSeq %in% unique(otherFrags[k,]$leaderSeq))
                     d[k2,]$leaderSeqGroup <<- paste0(x$trial[1], '~', x$subject[1], '~', g)
                   }

                   g <<- g + 1
                 }))
             
                 o <- left_join(o, select(d, leaderSeq, leaderSeqGroup), by = 'leaderSeq')
             } else {
               o$leaderSeqGroup <- paste0(x$trial[1], '~', x$subject[1], '~', 1)
             }
           
             o
             
             }))
         }))
}))


# Define fragment ids including ITR/LTR grouping ids and random ids.
frags$fragID <- paste0(frags$trial,     ':', frags$subject,    ':', frags$sample, ':',  
                       frags$replicate, ':', frags$chromosome, ':', frags$strand, ':', 
                       frags$fragStart, ':', frags$fragEnd,    ':', frags$leaderSeqGroup, ':',
                       frags$randomLinkerSeq)



# Create a tibble that can be turned into a GRange object which can be used to standardize positions.
# frags is on the read level, here we create f which tallies read counts for each fragment and we 
# standardize within subjects + ITR/LTR groupings. This will prevent closely spaced events with 
# different ITR/LTR remnants from being merged. Updated positions can be joined and updated via fragID.

f <- group_by(frags, fragID) %>% 
     summarise(seqnames = chromosome[1], start = fragStart[1], end = fragEnd[1], strand = strand[1], 
               reads = n_distinct(readID) + sum(n), fragID = fragID[1], s = leaderSeqGroup[1]) %>%
     ungroup()

g <- GenomicRanges::makeGRangesFromDataFrame(f, keep.extra.columns = TRUE)



# Standardize integration positions across subjects.

write(c(paste(now(), '   Standardizing integration positions within subjects.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)

g <- unlist(GenomicRanges::GRangesList(parallel::parLapply(cluster, split(g, g$s), function(x){
       library(dplyr)
       library(GenomicRanges)
       source(file.path(opt$softwareDir, 'lib.R'))
       source(file.path(opt$softwareDir, 'stdPos.lib.R'))
  
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


# Join updated positions to the input data frame.
frags <- left_join(frags, select(g, start, end, fragID), by = 'fragID')
rm(g)


# Report percent of fragments 
pIntSitePosUpdated <- sprintf("%.2f%%", (sum(unlist(lapply(split(frags, frags$strand), function(x){
       if(x$strand[1] == '+'){
         return(sum(x$fragStart != x$start))
       } else {
         return(sum(x$fragEnd != x$end))
       }
     }))) / nrow(frags))*100)

write(c(paste(now(), '   ', pIntSitePosUpdated, ' of read records updated with corrected intSite positions.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)


# Assign the new integration positions.
frags$fragStart <- frags$start
frags$fragEnd   <- frags$end
frags <- select(frags, -start, -end)


# Update position and fragment ids.
frags$posid <- paste0(frags$chromosome, frags$strand, 
                      ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd),
                      '.', stringr::str_extract(frags$leaderSeqGroup, '\\d+$'))


# Determine which read level fragment records map to multiple position ids.
u <- group_by(frags, readID) %>% mutate(nPosIDs = n_distinct(posid)) %>% ungroup() %>% filter(nPosIDs == 1) %>% pull(readID)
frags_uniqPosIDs <- frags[frags$readID %in% u,]   # (!) May have sonic break position chatter.
frags_multPosIDs <- frags[! frags$readID %in% u,]


# Remove UMIs found in unique frags from multiHit frags.
# UMIs will server as an abundance metrics in multihits.

# Check w/ real data...


# Create quick look-up index.
frags_uniqPosIDs$i <- paste(frags_uniqPosIDs$trial, frags_uniqPosIDs$subject, frags_uniqPosIDs$sample)

# Identify random ids that span more than one integration position.
r <- group_by(frags_uniqPosIDs, trial, subject, sample, randomLinkerSeq) %>% 
     mutate(n = n_distinct(posid)) %>% 
     ungroup() %>%
     filter(n > 1) %>%
     select(trial, subject, sample, randomLinkerSeq) %>%
     mutate(i = paste(trial, subject, sample))

if(nrow(r) > 0){
  # Quick sonicAbundance table based on non-standardized break points.
  f <- mutate(frags_uniqPosIDs, fragWidth = (fragEnd - fragStart) + 1) %>%
       select(trial, subject, sample, posid, fragWidth) %>%
       group_by(trial, subject, sample, posid) %>% 
       summarise(estAbund = n_distinct(fragWidth)) %>% 
       ungroup() %>% 
       mutate(i = paste(trial, subject, sample))
  
  invisible(lapply(split(r, r$randomLinkerSeq), function(x){
          # Retrieve all reads for this duplicated random linker sequence
          o <- subset(frags_uniqPosIDs, i == x$i[1] & randomLinkerSeq == x$randomLinkerSeq[1])
          t <- data.frame(sort(table(o$posid), decreasing = TRUE)) %>%
               left_join(subset(f, i == x$i[1]), by = c('Var1' = 'posid')) %>% 
               arrange(desc(Freq))
          
          # The top posid associated with this random id is read x times greater than the second - assign all to the first.
          if(nrow(t) > 1 & t[1,]$Freq >= t[2,]$Freq * opt$buildStdFragments_randomIDdupReadMult){
            ind <- which(frags_uniqPosIDs$posid != t[1,]$Var1 & frags_uniqPosIDs$randomLinkerSeq == x$randomLinkerSeq[1])
            frags_uniqPosIDs <<- frags_uniqPosIDs[-ind,]
            return()
          }
          
          # The top posid associated with this random id is more than x times abundant than the second - assign all to the first.
          t <- arrange(t, desc(estAbund))
          if(nrow(t) > 1 & t[1,]$estAbund >= t[2,]$estAbund * opt$buildStdFragments_randomIDdupAbundMult){
            ind <- which(frags_uniqPosIDs$posid != t[1,]$Var1 & frags_uniqPosIDs$randomLinkerSeq == x$randomLinkerSeq[1])
            frags_uniqPosIDs <<- frags_uniqPosIDs[-ind,]
            return()
          }
                  
         # We can not reasonably determine the true source of this random identifier, remove reads with this random id.
         ind <- which(frags_uniqPosIDs$randomLinkerSeq == x$randomLinkerSeq[1])
         frags_uniqPosIDs <<- frags_uniqPosIDs[-ind,]
         
         return()
       }))
}

# Determine if any random tags have more than one width.
f <- mutate(frags_uniqPosIDs, fragWidth = fragEnd - fragStart + 1) %>%
     group_by(trial, subject, sample, randomLinkerSeq) %>% 
     summarise(n = n_distinct(fragWidth), i = i[1]) %>% 
     ungroup()


# Use read counts to determine the most likely random fragment boundary.
if(any(f$n > 1)){
  f2 <- subset(f, n > 1)
  
  invisible(lapply(split(f2, f2$randomLinkerSeq), function(x){
    ind <- which(frags_uniqPosIDs$i == x$i & frags_uniqPosIDs$randomLinkerSeq == x$randomLinkerSeq)
    o <- frags_uniqPosIDs[ind,]
    
    if(o$strand[1] == '+'){
      frags_uniqPosIDs[ind,]$fragEnd <<- as.integer(names(sort(table(o$fragEnd), decreasing = TRUE))[1])
    } else {
      frags_uniqPosIDs[ind,]$fragStart <<- as.integer(names(sort(table(o$fragStart), decreasing = TRUE))[1])
    }
  }))
}

# Correct for instances where a read maps to more than fragment but all fragments have the same integration position.
# These are instances of fuzzy break points and here we select the shortest fragments lengths.

frags_uniqPosIDs <- tidyr::unite(frags_uniqPosIDs, fragID, trial, subject, sample, replicate, 
                                 chromosome, strand, fragStart, fragEnd, leaderSeqGroup, sep = ':', remove = FALSE)

z <- frags_uniqPosIDs$readID[duplicated(frags_uniqPosIDs$readID)]

if(length(z) > 0){
  a <- subset(frags_uniqPosIDs, readID %in% z)
  b <- subset(frags_uniqPosIDs, ! readID %in% z)
  
  a2 <- rbindlist(lapply(split(a, a$readID), function(x){
    if(n_distinct(x$fragID) > 1 & n_distinct(x$posid) == 1){
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
  
  frags_uniqPosIDs <- bind_rows(a2, b)
  rm(a, b, a2)
  gc()
}

# Update position and fragment IDs to reflect standardization.
frags_uniqPosIDs$posid <- paste0(frags_uniqPosIDs$chromosome, frags_uniqPosIDs$strand, 
                      ifelse(frags_uniqPosIDs$strand == '+', frags_uniqPosIDs$fragStart, frags_uniqPosIDs$fragEnd),
                      '.', stringr::str_extract(frags_uniqPosIDs$leaderSeqGroup, '\\d+$'))

frags_uniqPosIDs <- tidyr::unite(frags_uniqPosIDs, fragID, trial, subject, sample, replicate, 
                                 chromosome, strand, fragStart, fragEnd, leaderSeqGroup, sep = ':', remove = FALSE)


# Multihits
if(nrow(frags_multPosIDs) > 0){
  
  frags_multPosIDs$s <- paste(frags_multPosIDs$trial, frags_multPosIDs$subject)
  
  # Here we cycle through each sample or subject, create a list of uniquely resolved sites
  # and flag reads to be returned to frags if only one site in the list of potential sites is in the list uniquely resolved sites.
  
  frags_multPosIDs$returnToFrags <- FALSE
  
  frags_multPosIDs <- bind_rows(lapply(split(frags_multPosIDs, frags_multPosIDs$s), function(x){
                        uniqueSites <- unique(dplyr::filter(frags_uniqPosIDs, trial == x$trial[1], subject == x$subject[1])$posid)
    
                        bind_rows(lapply(split(x, x$readID), function(x2){
                          if(sum(unique(x2$posid) %in% uniqueSites) == 1) x2[x2$posid %in% uniqueSites,]$returnToFrags <- TRUE
                          x2
                        }))
                      }))
  
  # Return reads to frags if needed.
  
  if(any(frags_multPosIDs$returnToFrags == TRUE)){
    # Isolate frag reads that should be returned.
    m <- subset(frags_multPosIDs, returnToFrags == TRUE)

    # We correct for instances where a read supports multiple breaks but a single position (again)
    # because reads may of been missed earlier.
    
    m <- bind_rows(lapply(split(m, m$readID), function(x){
           if(n_distinct(x$fragID) > 1 & n_distinct(x$posid) == 1){
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

    msg <- paste0('   ', sprintf("%.2f%%", (nrow(m) / nrow(frags_multPosIDs))*100), ' of multihit reads recovered because one potential site was in the unabiguous site list.')
    write(c(paste(now(), msg)), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)

    # Remove salvaged reads from multi-hit reads.
    frags_multPosIDs <- subset(frags_multPosIDs, ! readID %in% m$readID) %>% dplyr::select(-s, -returnToFrags)

    # Return salvaged reads to frags.
    frags_uniqPosIDs <- bind_rows(frags_uniqPosIDs, dplyr::select(m, -s, -returnToFrags))
  }
} 

multiHitClusters <- tibble()

if(nrow(frags_multPosIDs) > 0 & opt$buildStdFragments_createMultiHitClusters){
  
  write(c(paste(now(), '   Building mulit-hit networks.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)
 
  # Create a data frame with width data needed by multi-hit calculating worker nodes
  # then export the data to the cluster nodes.
  multiHitFragWidths  <- mutate(frags_multPosIDs, width = fragEnd - fragStart + 1) %>%
                         group_by(trial, subject, sample, posid) %>%
                         summarise(reads = list(readID), UMIs = list(randomLinkerSeq)) %>%
                         ungroup() %>% 
                         distinct()
    
  clusterExport(cluster, 'multiHitFragWidths')
  
  # For each read, create a from -> to data frame and capture the width of the read. 
  multiHitNet_replicates <- rbindlist(lapply(split(frags_multPosIDs, frags_multPosIDs$readID), function(x){
    if(n_distinct(x$posid) == 1) return(tibble()) # Cases of break point only variation.
    
    # Create unique to - from permutations.
    node_pairs <- RcppAlgos::comboGeneral(unique(x$posid), 2)
    data.table(trial = x[1,]$trial, subject = x[1,]$subject, sample = x[1,]$sample, 
               replicate = x[1,]$replicate, from = node_pairs[,1], to = node_pairs[,2], readID = x[1,]$readID)
  }))

  # Create trial/subject/sample grouping indices that will be used to create a splitting vector.
  multiHitNet_replicates$n <- group_by(multiHitNet_replicates, trial, subject, sample) %>% group_indices
  
  # Create a second splitting vector for parallelization that will not break up trial/subject/sample groupings.
  d <- tibble(n = 1:max(multiHitNet_replicates$n), n2 = ntile(1:max(n), opt$buildStdFragments_CPUs))
  multiHitNet_replicates <- left_join(multiHitNet_replicates, d, by = 'n')

  multiHitClusters <- rbindlist(parLapply(cluster, split(multiHitNet_replicates, multiHitNet_replicates$n2), function(x){
  #multiHitClusters <- rbindlist(lapply(split(multiHitNet_replicates, multiHitNet_replicates$n2), function(x){
    library(igraph)
    library(dplyr)
    library(data.table)
    
    # Work within trial/subject/sample groupings...
    
    rbindlist(lapply(split(x, x$n), function(x2){
      # Build a graph with the posid nodes and read edges.
      g <- igraph::simplify(graph_from_data_frame(dplyr::select(x2, from, to), directed=FALSE, vertices=data.table(name = unique(c(x2$from, x2$to)))))
    
      # Separate out individual graphs.
      o <- data.table(trial = x2$trial[1], subject = x2$subject[1], sample = x2$sample[1], clusters = lapply(igraph::decompose(g), function(a) igraph::V(a)$name))
    
      # Create cluster ids.
      o$clusterID <- paste0('MHC.', 1:nrow(o))
      
      tidyr::unnest(o, clusters) %>%
      dplyr::left_join(subset(multiHitFragWidths, sample == x2$sample[1]) %>% dplyr::select(posid, reads, UMIs), by = c('clusters' = 'posid')) %>%
      tidyr::unnest(reads) %>%
      tidyr::unnest(UMIs) %>%
      dplyr::group_by(trial, subject, sample, clusterID) %>%
      dplyr::summarise(nodes = n_distinct(clusters), readIDs = list(unique(reads)), reads = n_distinct(reads), UMIs = n_distinct(UMIs), posids = list(unique(clusters))) %>%
      dplyr::ungroup() %>%
      dplyr::relocate(readIDs, .after = posids)
    }))
  }))
   
  rm(frags_multPosIDs, multiHitNet_replicates, multiHitFragWidths) 
  gc()
  
  saveRDS(multiHitClusters, file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitClusters.rds'))
  readr::write_tsv(multiHitClusters, file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitClusters.tsv.gz'))
  
  if('databaseGroup' %in% names(opt)){
    library(RMariaDB)
    
    dbConn <- tryCatch({
      dbConnect(RMariaDB::MariaDB(), group = opt$databaseGroup)
    },
    error=function(cond) {
      write(c(paste(now(), '   Error - could not connect to the database.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)
      if(opt$core_createFauxSiteDoneFiles) core_createFauxSiteDoneFiles()
      q(save = 'no', status = 1, runLast = FALSE) 
    })
    
    o <- unpackUniqueSampleID(sampleMetaData)
    multiHitClusters$i <- paste0(multiHitClusters$trial, '~', multiHitClusters$subject, '~', multiHitClusters$sample)
    o$i <- paste0(o$trial, '~', o$subject, '~', o$sample)
    multiHitClusters <- left_join(multiHitClusters, distinct(select(o, i, refGenome)) , by = 'i')
    
    invisible(lapply(split(multiHitClusters, multiHitClusters$i), function(x){
  
      dbExecute(dbConn, paste0("delete from multihits where trial='", x$trial[1], "' and subject='", x$subject[1],
                            "' and sample='", x$sample[1], "' and refGenome='", x$refGenome[1], "'"))
      
      f <- tmpFile()
      readr::write_tsv(dplyr::select(x, -trial, -subject, -sample, -i, -refGenome), file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'tmp', f))
      system(paste0('xz ', file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'tmp', f)))
      
      fp <- file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'tmp', paste0(f, '.xz'))
      
      tab <- readBin(fp, "raw", n = as.integer(file.info(fp)["size"])+100)
      
      invisible(file.remove(list.files(file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'tmp'), pattern = f, full.names = TRUE)))
      
      r <- dbExecute(dbConn,
                     "insert into multihits values (?, ?, ?, ?, ?, ?)",
                     params = list(x$trial[1], x$subject[1], x$sample[1], x$refGenome[1], list(serialize(tab, NULL)), as.character(lubridate::today())))
      
      if(r == 0){
        write(c(paste(now(), 'Error -- could not upload fragment data for ', x$sample[1], ' to the database.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)
        if(opt$core_createFauxSiteDoneFiles) core_createFauxSiteDoneFiles()
        q(save = 'no', status = 1, runLast = FALSE)
      } else {
        write(c(paste(now(), '   Uploaded multihit data for ', x$sample[1], ' to the database.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)
      }
    }))
    
    dbDisconnect(dbConn)
  }
}

# Frags are still read level. 
# Switch frags to a data frame because tibbles refuse to store single character vectors as lists.

frags <- group_by(data.frame(frags_uniqPosIDs), fragID) %>% mutate(i = n()) %>% ungroup()

a <- subset(frags, i == 1)
b <- subset(frags, i > 1)

if(nrow(b) > 0){
  o <- split(b, b$fragID)

  write(c(paste(now(), '   Bundling fragment reads into fragment records.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)

  b2 <- bind_rows(lapply(o, function(x){
         totalReads <- n_distinct(x$readID) + sum(x$n)
  
         if(totalReads < opt$buildStdFragments_minReadsPerFrag) return(tibble())
  
         # Select the leader sequence with the most reads as the representative leader sequence.
         x$repLeaderSeq <- names(sort(table(x$leaderSeq), decreasing = TRUE))[1]
  
         leaderSeqs <- x$repLeaderSeq
         readList <- x$readID
  
         x <- x[1,]
         x$reads <- totalReads
         x$readIDs <- list(readList)
         x$leaderSeqs <- list(leaderSeqs)
         x
       }))
} else {
  b2 <- tibble()
}

if(nrow(a) > 0){
  # Set values for fragments with single reads.
  a2 <- dplyr::group_by(a, fragID) %>%
        dplyr::mutate(reads = n+1, 
                      repLeaderSeq = leaderSeq[1], 
                      readIDs = list(readID),
                      leaderSeqs = list(leaderSeq[1])) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
} else {
  a2 <- tibble()
} 

f <- bind_rows(a2, b2)

# Clear out the tmp/ directory.
invisible(file.remove(list.files(file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'tmp'), full.names = TRUE)))

if (nrow(f) > 0) f <- left_join(f, sampleMetaData, by = 'uniqueSample')

s <- unique(paste0(f$trial, '~', f$subject, '~', f$sample))
if(any(! incomingSamples %in% s) & opt$core_createFauxSiteDoneFiles) core_createFauxSiteDoneFiles()

saveRDS(select(f, -leaderSeq, -n, -i, -id, -id, -trialSubject), file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'stdFragments.rds'))
readr::write_tsv(select(f, -leaderSeq, -n, -i, -id, -id, -trialSubject), file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'stdFragments.tsv.gz'))

q(save = 'no', status = 0, runLast = FALSE) 
