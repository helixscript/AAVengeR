library(dplyr)
library(lubridate)
library(parallel)
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(igraph)
library(dtplyr)

configFile <- commandArgs(trailingOnly=TRUE)

if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))
setMissingOptions()
setOptimalParameters()

cluster <- makeCluster(opt$buildStdFragments_CPUs)
clusterExport(cluster, 'opt')

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

frags <- tidyr::separate(frags, uniqueSample, c('trial', 'subject', 'sample', 'replicate'), sep = '~', remove = FALSE) %>% 
         dplyr::select(-refGenome, -vectorFastaFile, -seqRunID, -flags)

frags$trialSubject <- paste0(frags$trial, '~', frags$subject)


# Categorize ITR/LTR remnant sequences.
# -------------------------------------------------------------------------------
write(c(paste(now(), '   Categorizing leader sequences.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)

frags$posid <- paste0(frags$chromosome, frags$strand, ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd))


frags$fragID <- paste0(frags$trial,     ':', frags$subject,    ':', frags$sample, ':',  
                       frags$replicate, ':', frags$chromosome, ':', frags$strand, ':', 
                       frags$fragStart, ':', frags$fragEnd,    ':', frags$randomLinkerSeq)





# Find sites that are near one another with GRanges expanded +/- 100 and then 
# categorize LTR / ITR remnants.
#-------------------------------------------------------------------------------

frags <- data.table(frags)
positionExpansion <- opt$buildStdFragments_intSiteStdWindowWidth * 2

frags <- rbindlist(lapply(split(frags, paste(frags$trial, frags$subject)), function(k){
  
           k$start <- ifelse(k$strand == '+', k$fragStart - positionExpansion, k$fragEnd - positionExpansion)
           k$end   <- k$start + positionExpansion
           k$i <- 1:nrow(k)
  
           # Find the position ids with broad overlap.
           g <- makeGRangesFromDataFrame(distinct(select(k, posid, chromosome, strand, start, end, i)), keep.extra.columns = TRUE)
           g <- g[g$i]  # Ensure the same order as k.
  
           gr <- GenomicRanges::reduce(g, with.revmap = TRUE)

           o <- rbindlist(lapply(gr$revmap, function(x){
                  if(length(x) > 1){
                    return(data.table(n = list(x)))
                  } else {
                    return(NULL)
                  }
                }))
  
           k$leaderSeqGroupNum <- 1
  
           if(nrow(o) > 0){
             ke <- k[, c('readID', 'leaderSeq')]
             clusterExport(cluster, c('opt', 'ke', 'CD_HIT_clusters', 'tmpFile'), envir = environment())
    
             u <- rbindlist(parLapply(cluster, split(o, ntile(1:nrow(o), opt$buildStdFragments_CPUs)), function(d){
             #u <- rbindlist(lapply(split(o, ntile(1:nrow(o), opt$buildStdFragments_CPUs)), function(d){
                    library(GenomicRanges)
                    library(data.table)
                    library(Biostrings)
                    library(dplyr)
      
                    rbindlist(lapply(split(d, 1:nrow(d)), function(j){
                      a <- distinct(ke[unlist(j$n)])
      
                      if(n_distinct(a$leaderSeq) == 1){
                        return(data.table(readID = a$readID, leaderSeqGroupNum = 1))
                      } else {
                        d <- DNAStringSet(a$leaderSeq)
                        names(d) <- a$readID
                        clstrs <- CD_HIT_clusters(d, opt$outputDir, opt$buildStdFragments_remnantClusterParams)
        
                        n <- 0
                        m <- rbindlist(lapply(clstrs, function(x){
                               e <- unlist(stringr::str_extract_all(x, '>[^\\.]+'))
                        
                               if(length(e) > 0){
                                 n <<- n + 1
                                 return(data.table(readID = sub('^>', '', e), leaderSeqGroupNum2 = n))
                               } else {
                                 return(data.table())
                               }
                            }))
        
                        # Fail safe, should not be tripped.
                        if(any(duplicated(m$readID))) m <- m[! duplicated(m$readID),]
        
                        r <- left_join(a, distinct(m), by = 'readID')
        
                        i <- is.na(r$leaderSeqGroupNum2)
                        if(any(i)) r[is.na(r), ] <- 0
        
                        r$leaderSeqGroupNum <- r$leaderSeqGroupNum2
        
                        return(data.table(distinct(select(r, readID, leaderSeqGroupNum))))
                      }
                    }))
                  }))
    
             if(nrow(u) > 0){
                k$leaderSeqGroupNum <- NULL
                k <- left_join(k, distinct(u), by = 'readID', relationship = 'many-to-many') %>% dplyr::select(-start, -end, -i)
                k[is.na(k), ] <- 1  # Handle instances where u does not contain a read IDs. 
             }
           }
           
           k
}))

if('start' %in% names(frags)) frags$start <- NULL
if('end' %in% names(frags))   frags$end <- NULL

frags$leaderSeqGroup <- paste0(frags$trial, '~', frags$subject, '~', frags$leaderSeqGroupNum)
frags$leaderSeqGroupNum <- NULL

# Define fragment ids including ITR/LTR grouping ids and random ids.
frags$fragID <- paste0(frags$trial,     ':', frags$subject,    ':', frags$sample, ':',  
                       frags$replicate, ':', frags$chromosome, ':', frags$strand, ':', 
                       frags$fragStart, ':', frags$fragEnd,    ':', frags$leaderSeqGroup, ':',
                       frags$randomLinkerSeq)



# Create a tibble that can be turned into a GRange object which can be used to standardize positions.
# frags is on the read level, here we create f which tallies read counts for each fragment and we 
# standardize within subjects + ITR/LTR groupings. This will prevent closely spaced events with 
# different ITR/LTR remnants from being merged. Updated positions can be joined and updated via fragID.
#------------------------------------------------------------------------------------------------------

f <- group_by(frags, fragID) %>% 
     summarise(seqnames = chromosome[1], start = fragStart[1], end = fragEnd[1], strand = strand[1], 
               reads = n_distinct(readID) + sum(nDuplicateReads), fragID = fragID[1], s = leaderSeqGroup[1]) %>%
     ungroup()

g <- GenomicRanges::makeGRangesFromDataFrame(f, keep.extra.columns = TRUE)



# Standardize integration positions across subjects.
#-------------------------------------------------------------------------------

write(c(paste(now(), '   Standardizing integration positions within subjects.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)

g <- unlist(GenomicRanges::GRangesList(parallel::parLapply(cluster, split(g, g$s), function(x){
       library(dplyr)
       library(GenomicRanges)
       source(file.path(opt$softwareDir, 'lib.R'))
       source(file.path(opt$softwareDir, 'stdPos.lib.R'))
  
       x$intSiteRefined <- FALSE
  
       out <- tryCatch({
         o <- standardize_sites(x, counts.col = 'reads', sata.gap = opt$buildStdFragments_intSiteStdWindowWidth)
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


# Assign the new integration positions.
frags$fragStart <- frags$start
frags$fragEnd   <- frags$end
frags <- select(frags, -start, -end)


# Update position and fragment ids.
frags$posid <- paste0(frags$chromosome, frags$strand, 
                      ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd),
                      '.', stringr::str_extract(frags$leaderSeqGroup, '\\d+$'))

rm(g)



# Standardize break positions within replicates
#-------------------------------------------------------------------------------
frags$s <- paste0(frags$leaderSeqGroup, '~', frags$sample, '~', frags$replicate)

f <- group_by(frags, fragID) %>% 
     summarise(seqnames = chromosome[1], start = fragStart[1], end = fragEnd[1], strand = strand[1], 
              reads = n_distinct(readID) + sum(nDuplicateReads), fragID = fragID[1], s = s[1]) %>%
     ungroup()

g <- GenomicRanges::makeGRangesFromDataFrame(f, keep.extra.columns = TRUE)

g <- unlist(GenomicRanges::GRangesList(parallel::parLapply(cluster, split(g, g$s), function(x){
  library(dplyr)
  library(GenomicRanges)
  source(file.path(opt$softwareDir, 'lib.R'))
  source(file.path(opt$softwareDir, 'stdPos.lib.R'))
  
  x$breakPointsRefined <- FALSE
  
  out <- tryCatch({
    o <- refine_breakpoints(x, counts.col = 'reads', sata.gap = opt$buildStdFragments_breakPointStdWindowWidth)
    o$breakPointsRefined <- TRUE
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
frags$s <- NULL

# Join updated positions to the input data frame.
frags <- left_join(frags, select(g, start, end, fragID), by = 'fragID')


# Assign the new integration positions.
frags$fragStart <- frags$start
frags$fragEnd   <- frags$end
frags <- select(frags, -start, -end)

# Update position and fragment ids.
frags$posid <- paste0(frags$chromosome, frags$strand, 
                      ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd),
                      '.', stringr::str_extract(frags$leaderSeqGroup, '\\d+$'))


# Determine which read level fragment records map to multiple position ids.
#-------------------------------------------------------------------------------
u <- group_by(frags, readID) %>% 
     mutate(nPosIDs = n_distinct(posid)) %>% 
     ungroup() %>% filter(nPosIDs == 1)  %>% 
     pull(readID)

frags_uniqPosIDs <- frags[frags$readID %in% u,] 
frags_multPosIDs <- frags[! frags$readID %in% u,]

if(nrow(frags_uniqPosIDs) == 0){
  write(c(paste(now(), '   Error - No unique position remain after filtering.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)
  if(opt$core_createFauxSiteDoneFiles) core_createFauxSiteDoneFiles()
  q(save = 'no', status = 1, runLast = FALSE) 
}


frags_uniqPosIDs <- tidyr::unite(frags_uniqPosIDs, fragID, trial, subject, sample, replicate, 
                                 chromosome, strand, fragStart, fragEnd, leaderSeqGroup, sep = ':', remove = FALSE)

# Multihits
# This block is here because some reads may be returned to frags_uniqPosIDs.
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


# Ensure that the read level fragments for each sample position id have similar remnant sequences.
# Set leaderSeqGroups to 1 and enumerate if different remnants are found. This is a more rigorous
# clean-up from the original sorting before we standardized. 


buildConsensusSeq <- function(x){
  x$w <- x$fragEnd - x$fragStart
  tab <- group_by(x, leaderSeq) %>% 
         summarise(nWidths = n_distinct(fragEnd - fragStart), nReads = n()) %>% 
         ungroup() %>%
         arrange(desc(nWidths), desc(nReads))
  as.character(tab[1, 'leaderSeq'])
}

leaderSeqGroupIDs <- list()

frags_uniqPosIDs <- bind_rows(lapply(split(frags_uniqPosIDs, paste(frags_uniqPosIDs$trialSubject, frags_uniqPosIDs$sample, frags_uniqPosIDs$posid)), function(x){
  b <- sub('\\.\\d+$', '', x$posid[1])
  
  if(! b %in% names(leaderSeqGroupIDs)){
    leaderSeqGroupIDs[[b]] <<- 1
    x$leaderSeqGroupNum <- 1
  } else {
    leaderSeqGroupIDs[[b]] <<- leaderSeqGroupIDs[[b]] + 1
    x$leaderSeqGroupNum <- leaderSeqGroupIDs[[b]]
  }
  
  if(n_distinct(x$leaderSeq) != 1){

    d <- DNAStringSet(x$leaderSeq)
    names(d) <- x$readID
    clstrs <- CD_HIT_clusters(d, file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'tmp'), opt$buildStdFragments_remnantClusterParams)
    
    if(length(clstrs) == 1 | any(is.na(clstrs))){
      x$leaderSeq <- buildConsensusSeq(x)
    } else {
      n <- leaderSeqGroupIDs[[b]] - 1
      
      m <- rbindlist(lapply(clstrs, function(x){
        e <- unlist(stringr::str_extract_all(x, '>[^\\.]+'))
        if(length(e) > 0){
          n <<- n + 1
          return(data.table(readID = sub('^>', '', e), leaderSeqGroupNum = n))
        } else {
          return(data.table())
        }
      }))
      
      # Fail safe, should not be tripped.
      if(any(duplicated(m$readID))) m <- m[! duplicated(m$readID),]
      
      x$leaderSeqGroupNum <- NULL
      x <- left_join(x, m, by = 'readID')
      
      leaderSeqGroupIDs[[b]] <<- max(x$leaderSeqGroupNum)
      
      x <- bind_rows(lapply(split(x, x$leaderSeqGroupNum), function(x2){
             if(nrow(x2) > 1) x2$leaderSeq <- buildConsensusSeq(x2)
             x2
           }))
    }
  }
  
  x
}))

# Here we rearrange leader sequence groups so that the must abundant groups have the lowest identifiers.
frags_uniqPosIDs <- bind_rows(lapply(split(frags_uniqPosIDs, paste(frags_uniqPosIDs$trialSubject, frags_uniqPosIDs$sample, sub('\\.\\d+$', '', frags_uniqPosIDs$posid))), function(x){
  if(n_distinct(x$leaderSeqGroupNum) == 1){
   x$leaderSeqGroupNum <- 1 
  } else {
    x$w <- x$fragEnd - x$fragStart
    o <- group_by(x, leaderSeqGroupNum) %>% 
         summarise(nWidths = n_distinct(w), nReads = n()) %>% 
         ungroup() %>% 
         arrange(desc(nWidths), desc(nReads)) %>% 
         mutate(leaderSeqGroupNum2 = 1:n())
    
    x <- left_join(x, select(o, -nWidths, -nReads), by = 'leaderSeqGroupNum')
    x$leaderSeqGroupNum <- x$leaderSeqGroupNum2
    x <- select(x, -leaderSeqGroupNum2, -w)
  }
  x
}))


frags_uniqPosIDs$leaderSeqGroup <- paste0(frags_uniqPosIDs$trialSubject, '~', frags_uniqPosIDs$leaderSeqGroupNum)


# Update position and fragment IDs to reflect standardization.
frags_uniqPosIDs$posid <- paste0(frags_uniqPosIDs$chromosome, frags_uniqPosIDs$strand, 
                                 ifelse(frags_uniqPosIDs$strand == '+', frags_uniqPosIDs$fragStart, frags_uniqPosIDs$fragEnd),
                                 '.', frags_uniqPosIDs$leaderSeqGroupNum)

frags_uniqPosIDs$leaderSeqGroupNum <- NULL

# Create quick look-up index.
frags_uniqPosIDs$i <- paste(frags_uniqPosIDs$trial, frags_uniqPosIDs$subject, frags_uniqPosIDs$sample)


if(opt$processAdriftReadLinkerUMIs){

  frags_uniqPosIDs <- as.data.table(frags_uniqPosIDs)
  
  # Identify random ids that span more than one integration position.
  r <- group_by(lazy_dt(frags_uniqPosIDs), trial, subject, sample, randomLinkerSeq) %>% 
       mutate(n = n_distinct(posid)) %>% 
       ungroup() %>%
       filter(n > 1) %>%
       select(trial, subject, sample, randomLinkerSeq) %>%
       mutate(i = paste(trial, subject, sample)) %>%
       as.data.table()

  if(nrow(r) > 0){
    # Quick sonicAbundance table.
    f <- mutate(lazy_dt(frags_uniqPosIDs), fragWidth = (fragEnd - fragStart) + 1) %>%
         select(trial, subject, sample, posid, fragWidth) %>%
         group_by(trial, subject, sample, posid) %>% 
         summarise(estAbund = n_distinct(fragWidth)) %>% 
         ungroup() %>% 
         mutate(i = paste(trial, subject, sample)) %>%
         as.data.table()
  
    total <- n_distinct(r$randomLinkerSeq)
    counter <- 0
    
    invisible(lapply(split(r, r$randomLinkerSeq), function(x){
            message('UMI cleanup part 1: ', counter, '/', total); counter <<- counter + 1

            # Retrieve all reads for this duplicated random linker sequence
            o <- subset(frags_uniqPosIDs, i == x$i[1] & randomLinkerSeq == x$randomLinkerSeq[1])
            
            # The top posid associated with this random id is more than x times abundant than the second - assign all to the first.
            t <- data.frame(sort(table(o$posid), decreasing = TRUE)) %>%
                 left_join(subset(f, i == x$i[1]), by = c('Var1' = 'posid')) %>%
                 arrange(desc(estAbund))
            
            if(nrow(t) > 1 & t[1,]$estAbund >= t[2,]$estAbund * opt$buildStdFragments_randomIDdupAbundMult){
              ind <- which(frags_uniqPosIDs$posid != t[1,]$Var1 & frags_uniqPosIDs$randomLinkerSeq == x$randomLinkerSeq[1])
              frags_uniqPosIDs <<- frags_uniqPosIDs[-ind,]
              return()
            }
            
            t <- arrange(t, desc(Freq))
          
            # The top posid associated with this random id is read x times greater than the second - assign all to the first.
            if(nrow(t) > 1 & t[1,]$Freq >= t[2,]$Freq * opt$buildStdFragments_randomIDdupReadMult){
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
       ungroup() %>%
       as.data.table()
  
  # Use read counts to determine the most likely random fragment boundary.
  if(any(f$n > 1)){
    f2 <- subset(f, n > 1)
    
    total <- n_distinct(f2$randomLinkerSeq)
    counter <- 1
    
    invisible(lapply(split(f2, f2$randomLinkerSeq), function(x){
      message('UMI cleanup part 2: ', counter, '/', total); counter <<- counter + 1
      
      ind <- which(frags_uniqPosIDs$i == x$i & frags_uniqPosIDs$randomLinkerSeq == x$randomLinkerSeq)
      o <- frags_uniqPosIDs[ind,]
      
      if(o$strand[1] == '+'){
        frags_uniqPosIDs[ind,]$fragEnd <<- as.integer(names(sort(table(o$fragEnd), decreasing = TRUE))[1])
      } else {
        frags_uniqPosIDs[ind,]$fragStart <<- as.integer(names(sort(table(o$fragStart), decreasing = TRUE))[1])
      }
    }))
  }
  
  frags_uniqPosIDs <- as.data.frame(frags_uniqPosIDs)
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
  
  saveRDS(multiHitClusters, file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'multiHitClusters.rds'), compress = opt$compressDataFiles)
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
    
    multiHitClusters$i <- paste0(multiHitClusters$trial, '~', multiHitClusters$subject, '~', multiHitClusters$sample)

    o <- select(sampleMetaData, uniqueSample, refGenome) %>% mutate(i = sub('~\\d+$', '', uniqueSample)) %>% select(-uniqueSample) %>% distinct()
    multiHitClusters <- left_join(multiHitClusters, o, by = 'i')
    
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
                     "insert into multihits values (?, ?, ?, ?, ?)",
                     params = list(x$trial[1], x$subject[1], x$sample[1], x$refGenome[1], list(serialize(tab, NULL))))
      
      if(r == 0){
        write(c(paste(now(), 'Error -- could not upload multihit data for ', x$sample[1], ' to the database.')), file = file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log'), append = TRUE)
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
         totalReads <- n_distinct(x$readID) + sum(x$nDuplicateReads)
  
         if(totalReads < opt$buildStdFragments_minReadsPerFrag) return(tibble())
        
         #tab <- sort(table(x$leaderSeq), decreasing = TRUE)
         #x$leaderSeqs <- paste0(mapply(function(n, x){ paste0(n, ',', x, '|') }, names(tab), tab, SIMPLIFY = FALSE), collapse = '')
         
         x$repLeaderSeq <- x$leaderSeq[1]
         
         readList <- x$readID
         x <- x[1,]
         x$reads <- totalReads
         x$readIDs <- list(readList)
    
         x
       }))
} else {
  b2 <- tibble()
}

if(nrow(a) > 0){
  # Set values for fragments with single reads.
  a2 <- dplyr::group_by(a, fragID) %>%
        dplyr::mutate(reads = nDuplicateReads+1, 
                      repLeaderSeq = leaderSeq[1],
                      readIDs = list(readID)) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::filter(reads >= opt$buildStdFragments_minReadsPerFrag)
} else {
  a2 <- tibble()
} 

f <- bind_rows(a2, b2)

# Clear out the tmp/ directory.
invisible(file.remove(list.files(file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'tmp'), full.names = TRUE)))

if (nrow(f) > 0) f <- left_join(f, sampleMetaData, by = 'uniqueSample')

s <- unique(paste0(f$trial, '~', f$subject, '~', f$sample))
if(any(! incomingSamples %in% s) & opt$core_createFauxSiteDoneFiles) core_createFauxSiteDoneFiles()

saveRDS(select(f, -uniqueSample, -readID, -leaderSeq, -nDuplicateReads, -i, -leaderSeqGroup), file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'stdFragments.rds'), compress = opt$compressDataFiles)
readr::write_tsv(select(f, -uniqueSample, -readID, -leaderSeq, -nDuplicateReads, -i, -leaderSeqGroup), file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'stdFragments.tsv.gz'))

q(save = 'no', status = 0, runLast = FALSE) 
