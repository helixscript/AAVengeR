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

# Load fragment data from database or local file depending on configuration file.
# Fragments are read in on the read level and contain a variable (n) which is the number 
# of identical read pairs that were removed in prepReads. (n) can be added to fragment
# read counts to account for all reads supporting fragments.

write(c(paste(now(), '   Reading in fragment file(s).')), file = file.path(opt$outputDir, 'log'), append = TRUE)

# Create a database connection if requested in the configuration file.
if('databaseGroup' %in% names(opt)){
  library(RMariaDB)
  dbConn <- dbConnect(RMariaDB::MariaDB(), group = opt$databaseGroup)
  
  if(! 'dbConn' %in% ls()){
    write(c(paste(now(), 'Error -- could not connect to the database.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE)
  }
}

# Check for configuration errors.
# ------------------------------------------------------------------------------
if('buildStdFragments_autoPullTrialSamples' %in% names(opt) & ! 'databaseGroup' %in% names(opt)){
  if(opt$buildStdFragments_autoPullTrialSamples){
    write(c(paste(now(), 'Error -- the databaseGroup option must be provided with the buildStdFragments_autoPullTrialSamples option.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE)
  }
}

if('buildStdFragments_trialSubjectList' %in% names(opt) & ! 'databaseGroup' %in% names(opt)){
  write(c(paste(now(), 'Error -- the databaseGroup option must be provided with the buildStdFragments_trialSubjectList option.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE)
}



# Read in fragment data.
# ------------------------------------------------------------------------------
if('buildStdFragments_inputFile' %in% names(opt)){
  # There is an incoming fragment rds file and here we make a table of unique trial / subject pairings that 
  # will be used to pull in other samples from the database. 
  frags <- readRDS(file.path(opt$outputDir, opt$buildStdFragments_inputFile))
  
  if(opt$buildStdFragments_autoPullTrialSamples){
    trialSubjects <- unpackUniqueSampleID(tibble(uniqueSample = unique(frags$uniqueSample))) %>% dplyr::select(-uniqueSample, -replicate, -sample) %>% distinct()

    dbFrags <- bind_rows(lapply(split(trialSubjects, 1:nrow(trialSubjects)), function(x) pullDBsubjectFrags(dbConn, x$trial, x$subject)))
    
    dbFragsToAdd <- dbFrags[! dbFrags$uniqueSample %in% frags$uniqueSample,]
    
    if(nrow(dbFragsToAdd) > 0){
      write(c(paste(now(), paste0('   Adding ', nrow(dbFragsToAdd), ' fragment read records to incoming fragment data.'))), file = file.path(opt$outputDir, 'log'), append = TRUE)
      frags <- bind_rows(frags, dbFragsToAdd)
    }
  }
} else if('buildStdFragments_trialSubjectList' %in% names(opt)){
  # Format:  trial;subject|trial;subject, eg. 'Sabatino;pM50|Sabatino;pLinus'
  
  frags <- bind_rows(lapply(unlist(base::strsplit(opt$buildStdFragments_trialSubjectList, '\\|')), function(x){
         d <- unlist(base::strsplit(x, ';'))
         pullDBsubjectFrags(dbConn, d[1], d[2])
       }))

} else {
  write(c(paste(now(), 'Error -- neither buildStdFragments_inputFile or buildStdFragments_trialSubjectList options were provided')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE)
}

if(nrow(frags) == 0){
  write(c(paste(now(), 'Error -- no fragments were loaded or retrieved.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE)
}

if('databaseGroup' %in% names(opt)) RMariaDB::dbDisconnect(dbConn)


# Figure out why these are not distinct...
# Read ids may be duplicated for multi-hits.
frags <- distinct(frags)


# Set aside meta data till the end to save memory and add data back at the end 
# and unpack the uniqueSample ids.
# --------------------------------------------------------------------------------------------------------------
sampleMetaData <- distinct(dplyr::select(frags, uniqueSample, refGenome, vectorFastaFile, flags))

frags <- rbindlist(lapply(split(frags, frags$uniqueSample), function(x){
            o <-  unlist(stringr::str_split(x$uniqueSample[1], '~'))
            x$trial = o[1]; x$subject = o[2]; x$sample = o[3]; x$replicate = o[4]
            x
         })) %>% dplyr::select(-refGenome, -vectorFastaFile, -flags)
   
# Spark patch
# frags$leaderSeq <- sub('^TNC', 'TCC', frags$leaderSeq)





# Look at ITR/LTR remnants on the subject level, order by unique fragment widths.
# -------------------------------------------------------------------------------

write(c(paste(now(), '   Categorizing leader sequences.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

frags$posid <- paste0(frags$chromosome, frags$strand, ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd))


frags$fragID <- paste0(frags$trial,     ':', frags$subject,    ':', frags$sample, ':',  
                       frags$replicate, ':', frags$chromosome, ':', frags$strand, ':', 
                       frags$fragStart, ':', frags$fragEnd,    ':', frags$randomLinkerSeq)


frags <- bind_rows(lapply(split(frags, paste(frags$trial, frags$subject)), function(x){
           processed <- vector()
         
           x <- bind_rows(lapply(split(x, x$posid), function(x2){
           
           
           # Find neighboring site fragments.
           if(x2$strand[1] == '+'){
              o <- subset(x, chromosome == x2$chromosome[1] & strand == '+' & fragStart >= x2$fragStart[1] - opt$buildStdFragments_leaderSeqGroupingDist & fragStart <= x2$fragStart[1] + opt$buildStdFragments_leaderSeqGroupingDist)
           } else {
              o <- subset(x, chromosome == x2$chromosome[1] & strand == '-' & fragEnd >= x2$fragEnd[1] - opt$buildStdFragments_leaderSeqGroupingDist  & fragEnd <= x2$fragEnd[1] + opt$buildStdFragments_leaderSeqGroupingDist)
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
               dist <- stringdist::stringdist(d[i,]$leaderSeq, otherFrags$leaderSeq)

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
         
         x
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


cluster <- makeCluster(opt$buildStdFragments_CPUs)
clusterExport(cluster, 'opt')
write(c(paste(now(), '   Standardizing integration positions within subjects.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

# (!) Splitting on leader Seqs

g <- GenomicRanges::makeGRangesFromDataFrame(f, keep.extra.columns = TRUE)
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

frags <- left_join(frags, select(g, start, end, fragID), by = 'fragID')
rm(g)

pIntSitePosUpdated <- sprintf("%.2f%%", (sum(unlist(lapply(split(frags, frags$strand), function(x){
       if(x$strand[1] == '+'){
         return(sum(x$fragStart != x$start))
       } else {
         return(sum(x$fragEnd != x$end))
       }
     }))) / nrow(frags))*100)

write(c(paste(now(), '   ', pIntSitePosUpdated, ' of read records updated with corrected intSite positions.')), file = file.path(opt$outputDir, 'log'), append = TRUE)



# Assign the new integration positions.
# ------------------------------------------------------------------------------
frags$fragStart <- frags$start
frags$fragEnd   <- frags$end
frags <- select(frags, -start, -end)



# Update position and fragment ids.
# ------------------------------------------------------------------------------
frags$posid <- paste0(frags$chromosome, frags$strand, 
                      ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd),
                      '.', stringr::str_extract(frags$leaderSeqGroup, '\\d+$'))


# Update fragment ids with standardized positions.
frags$fragID <- paste0(frags$trial, ':', frags$subject, ':', frags$sample, ':', frags$replicate, ':',
                       frags$chromosome, ':', frags$strand, ':', frags$fragStart, ':', frags$fragEnd, ':', 
                       frags$leaderSeqGroup, ':', frags$randomLinkerSeq)




# Rather than standardizing across subjects, here we standardize breaks within 
# replicate / random id groups and parallelize across unique samples.
# ------------------------------------------------------------------------------
f <- group_by(frags, fragID) %>%
     summarise(seqnames = chromosome[1], 
               start = fragStart[1], 
               end = fragEnd[1], 
               strand = strand[1], 
               reads = n_distinct(readID) + sum(n), 
               fragID = fragID[1],
               s = paste0(uniqueSample[1], ':', posid[1], ':', randomLinkerSeq[1])) 



# Separate fragment records into those with a single fragment per UMI site (f1) 
# and multiple UMI fragments per site (f2) since we do not need to standardize 
# sites supported by single fragments.
# ------------------------------------------------------------------------------

z <- table(f$s)
f1 <- subset(f, ! s %in% names(z[z > 1]))
f2 <- subset(f, s %in% names(z[z > 1]))


if(nrow(f2) > 0){
  
  write(c(paste(now(), '   Standardizing sonic break postions within UMI / replicate groupings.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
  # First split data on site/UMI ids
  o <- split(f2, f2$s)
  
  i <- ntile(1:length(o), opt$buildStdFragments_CPUs)
  f2 <- bind_rows(lapply(1:length(o), function(x){
          a <- o[[x]]
          a$i <- i[x]
          a
        }))
  
  f2 <- bind_rows(parLapply(cluster, split(f2, f2$i), function(x){
          library(dplyr) 
          library(GenomicRanges)
          source(file.path(opt$softwareDir, 'lib.R'))
          source(file.path(opt$softwareDir, 'stdPos.lib.R'))
      
          g <- GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)

          g <- unlist(GenomicRanges::GRangesList(lapply(split(g, g$s), function(x){
            
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

f <- bind_rows(f1, f2)
rm(f1, f2)

# Join the standardized start/stop positions to the fragments data frame.
frags <- left_join(frags, select(f, start, end, fragID), by = 'fragID')


pBreakPosUpdated <- sprintf("%.2f%%", (sum(unlist(lapply(split(frags, frags$strand), function(x){
  if(x$strand[1] == '+'){
    return(sum(x$fragEnd != x$end))
  } else {
    return(sum(x$fragStart != x$start))
  }
}))) / nrow(frags))*100)

write(c(paste(now(), '   ', pBreakPosUpdated, ' of read records updated with corrected break positions.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

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
                       frags$leaderSeqGroup, ':', frags$randomLinkerSeq)


# Create an additional fragment ID without the random seq id.
frags$fragID2 <- paste0(frags$trial,      ':', frags$subject, ':', frags$sample,    ':', frags$replicate, ':', 
                        frags$chromosome, ':', frags$strand,  ':', frags$fragStart, ':', frags$fragEnd,    ':', 
                        frags$leaderSeqGroup)


# Correct for instances where a read maps to more than fragment but all fragments have the sample integration position.
# These are instances of fuzzy break points and here we select the shortest fragments lengths.

write(c(paste(now(), '   Correcting fuzzy break points.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

z <- frags$readID[duplicated(frags$readID)]

if(length(z) > 0){
  a <- subset(frags, readID %in% z)
  b <- subset(frags, ! readID %in% z)

  a2 <- rbindlist(lapply(split(a, a$readID), function(x){
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

  frags <- bind_rows(a2, b)
  rm(a, b, a2)
}

gc()


# Identify reads which map to more than position id and define these as multi-hit reads.
# --------------------------------------------------------------------------------------

write(c(paste(now(), '   Identifying multi-hit reads.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

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
  
  write(c(paste(now(), '   Identifying multi-hit reads that can be returned to list of uniquely called sites.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
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
    
    msg <- paste0('   ', sprintf("%.2f%%", (nrow(m) / nrow(multiHitFrags))*100), ' of multihit reads recovered because one potential site was in the unabiguous site list.')
    write(c(paste(now(), msg)), file = file.path(opt$outputDir, 'log'), append = TRUE)
    
    # Remove salvaged reads from multi-hit reads.
    multiHitFrags <- subset(multiHitFrags, ! readID %in% m$readID) %>% dplyr::select(-s, -returnToFrags)
    
    # Return salvaged reads to frags.
    frags <- bind_rows(frags, dplyr::select(m, -s, -returnToFrags))
  }
} 

rm(o)
gc()

multiHitClusters <- tibble()

if(nrow(multiHitFrags) > 0 & opt$buildStdFragments_createMultiHitClusters){
  
  write(c(paste(now(), '   Building mulit-hit networks.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
 
  # Read and width data needed by multi-hit calculating worker nodes.
  multiHitFragWidths  <- mutate(multiHitFrags, width = fragEnd - fragStart + 1) %>%
                         group_by(trial, subject, sample, posid) %>%
                         summarise(widths = list(width), reads = list(readID)) %>%
                         ungroup() %>% 
                         distinct()
    
  clusterExport(cluster, 'multiHitFragWidths')
  
  # For each read, create a from -> to data frame and capture the width of the read. 
  multiHitNet_replicates <- rbindlist(lapply(split(multiHitFrags, multiHitFrags$readID), function(x){
    if(n_distinct(x$posid) == 1) return(tibble()) # Cases of break point only variation.
    node_pairs <- RcppAlgos::comboGeneral(unique(x$posid), 2)
    data.table(trial = x[1,]$trial, subject = x[1,]$subject, sample = x[1,]$sample, 
           replicate = x[1,]$replicate, from = node_pairs[,1], to = node_pairs[,2], readID = x[1,]$readID)
  }))

  # Create trial/subject/sample grouping indices that will be used to create a splitting vector for parallelization.
  multiHitNet_replicates$n <- group_by(multiHitNet_replicates, trial, subject, sample) %>% group_indices
  d <- tibble(n = 1:max(multiHitNet_replicates$n), n2 = ntile(1:max(n), opt$buildStdFragments_CPUs))
  multiHitNet_replicates <- left_join(multiHitNet_replicates, d, by = 'n')

  
  multiHitClusters <- rbindlist(parLapply(cluster, split(multiHitNet_replicates, multiHitNet_replicates$n2), function(x){
  #multiHitClusters <- rbindlist(lapply(split(multiHitNet_replicates, multiHitNet_replicates$n2), function(x){
    library(igraph)
    library(data.table)
    
    rbindlist(lapply(split(x, x$n), function(x2){
      # Build a graph with the posid nodes and read edges.
      g <- igraph::simplify(graph_from_data_frame(dplyr::select(x2, from, to), directed=FALSE, vertices=data.table(name = unique(c(x2$from, x2$to)))))
    
      # Separate out individual graphs.
      o <- data.table(trial = x2$trial[1], subject = x2$subject[1], sample = x2$sample[1], clusters = lapply(igraph::decompose(g), function(a) igraph::V(a)$name))
    
      # Create cluster ids.
      o$clusterID <- paste0('MHC.', 1:nrow(o))
      
      tidyr::unnest(o, clusters) %>%
      dplyr::left_join(subset(multiHitFragWidths, sample == x2$sample[1]) %>% dplyr::select(posid, reads, widths), by = c('clusters' = 'posid')) %>%
      tidyr::unnest(reads) %>%
      tidyr::unnest(widths) %>%
      dplyr::group_by(trial, subject, sample, clusterID) %>%
      dplyr::summarise(nodes = n_distinct(clusters), reads = n_distinct(reads), abund = n_distinct(widths), posids = list(unique(clusters))) %>%
      dplyr::ungroup()
    }))
  }))
   
  rm(multiHitFrags, multiHitNet_replicates, multiHitFragWidths) 
  gc()
  
  saveRDS(multiHitClusters, file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitClusters.rds'))
  
  
  if('databaseGroup' %in% names(opt)){
    library(RMariaDB)
    
    dbConn <- dbConnect(RMariaDB::MariaDB(), group = opt$databaseGroup)
    
    if(! 'dbConn' %in% ls()){
      write(c(paste(now(), 'Error -- could not connect to the database.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
      q(save = 'no', status = 1, runLast = FALSE)
    }
    
    o <- unpackUniqueSampleID(sampleMetaData)
    multiHitClusters$i <- paste0(multiHitClusters$trial, '~', multiHitClusters$subject, '~', multiHitClusters$sample)
    o$i <- paste0(o$trial, '~', o$subject, '~', o$sample)
    multiHitClusters <- left_join(multiHitClusters, distinct(select(o, i, refGenome)) , by = 'i')
    
    invisible(lapply(split(multiHitClusters, multiHitClusters$i), function(x){
  
      dbExecute(dbConn, paste0("delete from multihits where trial='", x$trial[1], "' and subject='", x$subject[1],
                            "' and sample='", x$sample[1], "' and refGenome='", x$refGenome[1], "'"))
      
      f <- tmpFile()
      readr::write_tsv(dplyr::select(x, -trial, -subject, -sample, -i, -refGenome), file.path(opt$outputDir, 'tmp', f))
      system(paste0('xz ', file.path(opt$outputDir, 'tmp', f)))
      
      fp <- file.path(opt$outputDir, 'tmp', paste0(f, '.xz'))
      
      tab <- readBin(fp, "raw", n = as.integer(file.info(fp)["size"])+100)
      
      invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), pattern = f, full.names = TRUE)))
      
      r <- dbExecute(dbConn,
                     "insert into multihits values (?, ?, ?, ?, ?, ?)",
                     params = list(x$trial[1], x$subject[1], x$sample[1], x$refGenome[1], list(serialize(tab, NULL)), as.character(lubridate::today())))
      
      if(r == 0){
        write(c(paste(now(), 'Error -- could not upload fragment data for ', x$sample[1], ' to the database.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
        q(save = 'no', status = 1, runLast = FALSE)
      } else {
        write(c(paste(now(), '   Uploaded multihit data for ', x$sample[1], ' to the database.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
      }
    }))
  }
}

# Frags are still read level. 
# Switch frags to a data frame because tibbles refuse to store single character vectors as lists.

frags <- data.frame(frags)


# Here we split the read level data by boundary corrected fragment ids, 
# ....

frags <- group_by(frags, fragID) %>% mutate(i = n()) %>% ungroup()
a <- subset(frags, i == 1)
b <- subset(frags, i > 1)
o <- split(b, b$fragID)

write(c(paste(now(), '   Bundling fragment reads into fragment records.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

f <- bind_rows(lapply(o, function(x){
  
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


# Set values for fragments with single reads.
a2 <- dplyr::group_by(a, fragID) %>%
      dplyr::mutate(reads = n+1, 
                   repLeaderSeq = leaderSeq[1], 
                   readIDs = list(readID),
                   leaderSeqs = list(leaderSeq[1])) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
            
f <- bind_rows(f, a2)


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

write(c(paste(now(), '   Correcting instances where the same UMI is associated with more than one sample fragment.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

f$i <- NULL

f2 <- bind_rows(lapply(split(f, paste(f$trial, f$subject, f$sample)), function(x){
       t <- table(x$randomLinkerSeq) 
       t <- names(t[which(t > 1)])
       
       if(length(t) > 0){
         a <- subset(x, ! randomLinkerSeq %in% t)
         b <- subset(x, randomLinkerSeq %in% t)
         
         b$randomLinkerSeq.dispute <- TRUE
    
         b2 <- bind_rows(lapply(t, function(i){
                 o <- dplyr::arrange(subset(b, randomLinkerSeq == i), desc(reads))
                 o$percentReads <- (o$reads / sum(o$reads))*100

                 if(o[1,]$percentReads < opt$buildStdFragments_UMI_conflict_minPercentReads){
                   o$msg <- 'All fragments rejected'
                   o$percentReads <- sprintf("%.2f%%", o$percentReads)
                   readr::write_tsv(dplyr::select(o, trial, subject, sample, replicate, readID, randomLinkerSeq, posid, percentReads, msg), 
                                    file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'randomIDexcludedReads', paste0('duplicate_UMIs_withinSample~', x$trial[1], '~', x$subject[1], '~', x$sample[1], '.tsv')))
                   return(tibble())
                 } else {
                   o$msg <- 'Top fragment retained'
                   o$percentReads <- sprintf("%.2f%%", o$percentReads)
                   readr::write_tsv(dplyr::select(o[2:nrow(o),], trial, subject, sample, replicate, readID, randomLinkerSeq, posid, percentReads, msg), 
                                    file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'randomIDexcludedReads', paste0('duplicate_UMIs_withinSample~', x$trial[1], '~', x$subject[1], '~', x$sample[1], '.tsv')))
                 }
                 
                 o[1,]
               }))
       
         x <- bind_rows(a, b2)
       }
    
       x
     }))

f2 <- dplyr::select(f2, -n, -fragID, -fragID2, -readID, -leaderSeq)
f2 <- left_join(f2, sampleMetaData, by = 'uniqueSample') %>% dplyr::select(-uniqueSample)

saveRDS(f2, file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'stdFragments.rds'))

q(save = 'no', status = 0, runLast = FALSE) 
