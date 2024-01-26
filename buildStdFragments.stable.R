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
set.seed(1)

cluster <- makeCluster(opt$buildStdFragments_CPUs)
clusterSetRNGStream(cluster, 1)
clusterExport(cluster, 'opt')

dir.create(file.path(opt$outputDir, opt$buildStdFragments_outputDir))
dir.create(file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'tmp'))

opt$defaultLogFile <- file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'log')


# Load fragment data from database or local file depending on configuration file.
# Fragments are read in on the read level and contain a variable (n) which is the number 
# of identical read pairs that were removed in prepReads. (n) can be added to fragment
# read counts to account for all reads supporting fragments.
#----------------------------------------------------------------------------------------

updateLog('Reading in fragment file(s).')

# Create a database connection if requested in the configuration file.
if('databaseGroup' %in% names(opt)){
  library(RMariaDB)

  dbConn <- tryCatch({
    dbConnect(RMariaDB::MariaDB(), group = opt$databaseGroup)
  },
  error=function(cond) {
    updateLog('Error - could not connect to the database.')
    q(save = 'no', status = 1, runLast = FALSE) 
  })
}

# Check for configuration errors.
# ------------------------------------------------------------------------------
if('buildStdFragments_autoPullTrialSamples' %in% names(opt) & ! 'databaseGroup' %in% names(opt)){
  if(opt$buildStdFragments_autoPullTrialSamples){
    updateLog('Error -- the databaseGroup option must be provided with the buildStdFragments_autoPullTrialSamples option.')
    q(save = 'no', status = 1, runLast = FALSE)
  }
}

if('buildStdFragments_trialSubjectList' %in% names(opt) & ! 'databaseGroup' %in% names(opt)){
  updateLog('Error -- the databaseGroup option must be provided with the buildStdFragments_trialSubjectList option.')
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
      updateLog(paste0('   Adding ', nrow(dbFragsToAdd), ' fragment read records to incoming fragment data.'))
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
  updateLog('Error -- neither buildStdFragments_inputFile or buildStdFragments_trialSubjectList options were provided.')
  q(save = 'no', status = 1, runLast = FALSE)
}

if('databaseGroup' %in% names(opt)) RMariaDB::dbDisconnect(dbConn)

# Make sure fragments were retrieved.
if(nrow(frags) == 0){
  updateLog('Error -- no fragments were loaded or retrieved.')
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
updateLog('Categorizing leader sequences.')

frags$posid <- paste0(frags$chromosome, frags$strand, ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd))



# Categorize LTR / ITR remnants.
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

frags <- tidyr::unite(frags, fragID, trial, subject, sample, replicate, chromosome, 
                      strand, fragStart, fragEnd, leaderSeqGroup, randomLinkerSeq, sep = ':', remove = FALSE)

# Create an index upon which to standardize integration positions within subjects.
frags <- tidyr::unite(frags, z, trial, subject, chromosome, 
                      strand, fragStart, fragEnd, leaderSeqGroup, sep = ':', remove = FALSE)

# Create a tibble that can be turned into a GRange object which can be used to standardize positions.
# frags is on the read level, here we create f which tallies read counts for each fragment and we 
# standardize within subjects + ITR/LTR groupings. This will prevent closely spaced events with 
# different ITR/LTR remnants from being merged. Updated positions can be joined and updated via fragID.
#------------------------------------------------------------------------------------------------------

f <- lazy_dt(frags) %>%
     group_by(z) %>%
     summarise(seqnames = chromosome[1], 
                          start    = fragStart[1], 
                          end      = fragEnd[1], 
                          strand   = strand[1], 
                          reads    = n_distinct(readID) + sum(nDuplicateReads), 
                          s        = paste0(trialSubject[1], '~', leaderSeqGroup[1], '~', chromosome[1])) %>%
    ungroup() %>%
    as.data.table()
     
# Create GenomicRanges object for standardization.
g <- GenomicRanges::makeGRangesFromDataFrame(f, keep.extra.columns = TRUE)


# Standardize integration positions within subjects.
#-------------------------------------------------------------------------------

updateLog('Standardizing integration positions within subjects.')

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
frags <- lazy_dt(frags) %>% left_join(distinct(select(g, start, end, z)), by = 'z') %>% as.data.table()


# Assign the new integration positions.
frags$fragStart <- frags$start
frags$fragEnd   <- frags$end
frags <- select(frags, -start, -end)


# Update position and fragment ids.
frags$posid <- paste0(frags$chromosome, frags$strand, 
                      ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd),
                      '.', stringr::str_extract(frags$leaderSeqGroup, '\\d+$'))

frags <- tidyr::unite(frags, fragID, trial, subject, sample, replicate, 
                      chromosome, strand, fragStart, fragEnd, leaderSeqGroup, randomLinkerSeq, sep = ':', remove = FALSE)


# Create an index upon which to standardize integration positions within specific replicate level integration events.
frags <- tidyr::unite(frags, z, trial, subject, chromosome, replicate,
                      strand, fragStart, fragEnd, posid, sep = ':', remove = FALSE)


rm(f, g)


# Standardize break positions within replicates.
#-------------------------------------------------------------------------------

updateLog('Standardizing break positions within replicates.')

f <- lazy_dt(frags) %>%
  group_by(z) %>%
  summarise(seqnames = chromosome[1], 
            start    = fragStart[1], 
            end      = fragEnd[1], 
            strand   = strand[1], 
            posid    = posid[1],
            reads    = n_distinct(readID) + sum(nDuplicateReads), 
            s        = paste0(uniqueSample[1], '~', chromosome[1])) %>%
  ungroup() %>%
  as.data.table()



g <- parallel::parLapply(cluster, split(f, f$s), function(k){
       library(dplyr)
       library(GenomicRanges)
       source(file.path(opt$softwareDir, 'lib.R'))
       source(file.path(opt$softwareDir, 'stdPos.lib.R'))
  
       k <- GenomicRanges::makeGRangesFromDataFrame(k, keep.extra.columns = TRUE)
  
       unlist(GRangesList(lapply(split(k, k$posid), function(x){
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
     })

g <- data.frame(unlist(GRangesList(g)))       


# Join updated positions to the input data frame.
frags <- lazy_dt(frags) %>% left_join(distinct(select(g, start, end, z)), by = 'z') %>% as.data.table()


# Assign the new integration positions.
frags$fragStart <- frags$start
frags$fragEnd   <- frags$end
frags <- select(frags, -start, -end)

updateLog('Updating positions and fragment ids.')

# Update position and fragment ids.
frags$posid <- paste0(frags$chromosome, frags$strand, 
                      ifelse(frags$strand == '+', frags$fragStart, frags$fragEnd),
                      '.', stringr::str_extract(frags$leaderSeqGroup, '\\d+$'))

frags <- tidyr::unite(frags, fragID, trial, subject, sample, replicate, 
                      chromosome, strand, fragStart, fragEnd, leaderSeqGroup, randomLinkerSeq, sep = ':', remove = FALSE)

rm(f, g)


# Determine which read level fragment records map to multiple position ids.
#-------------------------------------------------------------------------------

updateLog('Idenitfying uniquely called postions.')

u <- group_by(frags, readID) %>% 
     mutate(nPosIDs = n_distinct(posid)) %>% 
     ungroup() %>% filter(nPosIDs == 1)  %>% 
     pull(readID)

frags_uniqPosIDs <- frags[frags$readID %in% u,] 
frags_multPosIDs <- frags[! frags$readID %in% u,]



# Do not continue unless we have at least one uniquely called fragment.
if(nrow(frags_uniqPosIDs) == 0){
  updateLog('Error - No unique position remain after filtering.')
 
  if(opt$core_createFauxSiteDoneFiles) core_createFauxSiteDoneFiles()
  q(save = 'no', status = 1, runLast = FALSE) 
}



# Multihits
# This block is here because some reads may be returned to frags_uniqPosIDs.
if(nrow(frags_multPosIDs) > 0){
  
  updateLog('Resolving multi-hit reads.')
  
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
    updateLog(msg)
    
    # Remove salvaged reads from multi-hit reads.
    frags_multPosIDs <- subset(frags_multPosIDs, ! readID %in% m$readID) %>% dplyr::select(-s, -returnToFrags)
    
    # Return salvaged reads to frags.
    frags_uniqPosIDs <- bind_rows(frags_uniqPosIDs, dplyr::select(m, -s, -returnToFrags))
  }
} 



# Correct for instances where a read maps to more than fragment but all fragments have the same integration position.
# These are instances of fuzzy break points and here we select the shortest fragments lengths.

z <- frags_uniqPosIDs$readID[duplicated(frags_uniqPosIDs$readID)]

if(length(z) > 0){
  updateLog('Correcting fuzzy break points.')
  
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

updateLog('Reevaluating leader sequence classification after standardization.')

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
# This is most applicable to AAV work.
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

frags_uniqPosIDs$posid <- paste0(frags_uniqPosIDs$chromosome, frags_uniqPosIDs$strand, 
                                 ifelse(frags_uniqPosIDs$strand == '+', frags_uniqPosIDs$fragStart, frags_uniqPosIDs$fragEnd),
                                 '.', frags_uniqPosIDs$leaderSeqGroupNum)

frags_uniqPosIDs$leaderSeqGroupNum <- NULL


# Create quick look-up index used for UMI correction to keep UMI reassignments within samples.
frags_uniqPosIDs$i <- paste(frags_uniqPosIDs$trial, frags_uniqPosIDs$subject, frags_uniqPosIDs$sample)


if(opt$processAdriftReadLinkerUMIs){
  
  updateLog('Processing UMIs.')
  
  frags_uniqPosIDs <- as.data.table(frags_uniqPosIDs)
  
  # Identify random ids that span more than one integration position within samples.
  r <- group_by(lazy_dt(frags_uniqPosIDs), trial, subject, sample, randomLinkerSeq) %>% 
       mutate(n = n_distinct(posid)) %>% 
       ungroup() %>%
       filter(n > 1) %>%
       select(trial, subject, sample, randomLinkerSeq) %>%
       mutate(i = paste(trial, subject, sample)) %>%
       distinct() %>%
       as.data.table()

  if(nrow(r) > 0){
    # Quick sonicAbundance table.
    f <- mutate(lazy_dt(frags_uniqPosIDs), fragWidth = (fragEnd - fragStart) + 1) %>%
         select(trial, subject, sample, posid, fragWidth) %>%
         group_by(trial, subject, sample, posid) %>% 
         summarise(estAbund = n_distinct(fragWidth)) %>% 
         ungroup() %>% 
         mutate(i = paste(trial, subject, sample)) %>%
         distinct() %>%
         as.data.table()
  
    total <- n_distinct(r$randomLinkerSeq)
    counter <- 0
    
    abund_correction <- 0
    readCount_correction <- 0
    no_correction <- 0
    
    updateLog(paste0('Evaluating ', total, ' UMI sequences.'))
    
    inds <- unlist(lapply(split(r, r$randomLinkerSeq), function(x){
            if(counter %% 10 == 0) message('UMI cleanup part 1: ', counter, '/', total, ' - ', sprintf("%.2f%%", (counter/total)*100), ' - ',
                                           abund_correction, ' abund corrections; ', 
                                           readCount_correction, ' read count corrections; ',
                                           no_correction, ' no correction possible.')
            counter <<- counter + 1

            # Retrieve all reads for this duplicated random linker sequence
            t <- lazy_dt(frags_uniqPosIDs) %>% 
                 dplyr::filter(i == x$i[1] & randomLinkerSeq == x$randomLinkerSeq[1]) %>%
                 dplyr::group_by(posid) %>%
                 dplyr::summarise(count = n()) %>%
                 dplyr::left_join(subset(f, i == x$i[1]), by = 'posid') %>%
                 dplyr::arrange(desc(estAbund)) %>%
                 as.data.table()
            
            # The top posid associated with this random id is more than x times abundant than the second - assign all to the first.
            if(nrow(t) > 1 & t[1,]$estAbund >= t[2,]$estAbund * opt$buildStdFragments_randomIDdupAbundMult){
              ind <- which(frags_uniqPosIDs$posid != t[1,]$posid & frags_uniqPosIDs$randomLinkerSeq == x$randomLinkerSeq[1])
              abund_correction <<- abund_correction + 1
              return(ind)
            }
            
            t <- arrange(t, desc(count))
          
            # The top posid associated with this random id is read x times greater than the second - assign all to the first.
            if(nrow(t) > 1 & t[1,]$count >= t[2,]$count * opt$buildStdFragments_randomIDdupReadMult){
              ind <- which(frags_uniqPosIDs$posid != t[1,]$posid & frags_uniqPosIDs$randomLinkerSeq == x$randomLinkerSeq[1])
              readCount_correction <<- readCount_correction + 1
              return(ind)
            }
          
            # We can not reasonably determine the true source of this random identifier, remove reads with this random id.
            ind <- which(frags_uniqPosIDs$randomLinkerSeq == x$randomLinkerSeq[1])
            no_correction <<- no_correction + 1
         
            return(ind)
         }))
    
    if(length(inds) > 0) frags_uniqPosIDs <- frags_uniqPosIDs[unique(unname(inds)) * -1]
  }
  
  # Determine if any random tags have more than one width.
  # Each random sequence should be associated with a single posid after previous cleanup.
  
  f <- lazy_dt(frags_uniqPosIDs) %>%
       mutate(fragWidth = fragEnd - fragStart + 1) %>%
       group_by(trial, subject, sample, randomLinkerSeq) %>% 
       summarise(n = n_distinct(fragWidth), p = n_distinct(posid), i = i[1]) %>% 
       ungroup() %>%
       as.data.table()
  
  
  
  # Use read counts to determine the most likely random fragment boundary.
  # A bit slow due to inplace row updates.
  if(any(f$n > 1)){
    f2 <- subset(f, n > 1)
    
    total <- n_distinct(f2$randomLinkerSeq)
    counter <- 1
    
    invisible(lapply(split(f2, f2$randomLinkerSeq), function(x){
      message('UMI cleanup part 2: ', counter, '/', total); counter <<- counter + 1
      
      ind <- which(frags_uniqPosIDs$i == x$i & frags_uniqPosIDs$randomLinkerSeq == x$randomLinkerSeq)
      o <- frags_uniqPosIDs[ind]
      
      if(o$strand[1] == '+'){
        frags_uniqPosIDs[ind]$fragEnd <<- as.integer(names(sort(table(o$fragEnd), decreasing = TRUE))[1])
      } else {
        frags_uniqPosIDs[ind]$fragStart <<- as.integer(names(sort(table(o$fragStart), decreasing = TRUE))[1])
      }
    }))
  }
  
  frags_uniqPosIDs <- as.data.frame(frags_uniqPosIDs)
}



# Build multi-hit clusters if requested.
multiHitClusters <- tibble()

if(nrow(frags_multPosIDs) > 0 & opt$buildStdFragments_createMultiHitClusters){
  
  updateLog('Building mulit-hit networks.')
 
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
      updateLog('Error - could not connect to the database.')
     
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
        updateLog(paste0('Error -- could not upload multihit data for ', x$sample[1], ' to the database.'))
        
        if(opt$core_createFauxSiteDoneFiles) core_createFauxSiteDoneFiles()
        q(save = 'no', status = 1, runLast = FALSE)
      } else {
        updateLog(paste0('Uploaded multihit data for ', x$sample[1], ' to the database.'))
      }
    }))
    
    dbDisconnect(dbConn)
  }
}


# Remove randomLinkerSeq from fragIDs.
frags_uniqPosIDs <- tidyr::unite(frags_uniqPosIDs, fragID, trial, subject, sample, replicate, 
                      chromosome, strand, fragStart, fragEnd, leaderSeqGroup, sep = ':', remove = FALSE)


frags <- group_by(data.frame(frags_uniqPosIDs), fragID) %>% mutate(i = n()) %>% ungroup()

a <- subset(frags, i == 1)   
b <- subset(frags, i > 1)    

# Set values for fragments with more than one read.

# If you include the UMIs in the fragID, frags will be split into
# many smaller pieces and when you take the first row, all the UMIs 
# will be returned. If UMIs are not included in the fragID, then an
# arbitrary UMI will be returned for each fragment length.

filterUMIs <- function(x){

  k <- data.frame(table(x)) %>%
       dplyr::mutate(f = (Freq / sum(Freq) * 100)) %>%
       dplyr::filter(f >= opt$buildStdFragments_UMIminPercentReads) %>%
       dplyr::pull(x) %>%
       as.character()
  
  if(length(k) == 0) k <- vector()
  
  k
}



if(nrow(b) > 0){
  o <- split(b, b$fragID)

  updateLog('Bundling fragment reads into fragment records.')
 
  b2 <- bind_rows(lapply(o, function(x){
         totalReads <- n_distinct(x$readID) + sum(x$nDuplicateReads)
  
         if(totalReads < opt$buildStdFragments_minReadsPerFrag) return(tibble())
         
         x$repLeaderSeq <- x$leaderSeq[1]
         
         readList <- x$readID
         
         rUmiList <- unique(x$randomLinkerSeq)
         fUmiList <- filterUMIs(x$randomLinkerSeq)
         
         x <- x[1,]
         
         x$reads <- totalReads
         
         x$readIDs   <- list(readList)
         x$rUMI_list <- list(rUmiList)
         x$fUMI_list <- list(fUmiList)
    
         x
       }))
} else {
  b2 <- tibble()
}


if(nrow(a) > 0){
  a2 <- dplyr::group_by(a, fragID) %>%
        dplyr::mutate(reads = nDuplicateReads+1, 
                      repLeaderSeq = leaderSeq[1],
                      readIDs = list(readID),
                      rUMI_list = list(randomLinkerSeq),
                      fUMI_list = list(filterUMIs(randomLinkerSeq))) %>%
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

saveRDS(select(f, -z, -uniqueSample, -readID, -leaderSeq, -nDuplicateReads, -i, -leaderSeqGroup, -randomLinkerSeq), file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'stdFragments.rds'), compress = opt$compressDataFiles)
readr::write_tsv(select(f, -z, -uniqueSample, -readID, -leaderSeq, -nDuplicateReads, -i, -leaderSeqGroup, -randomLinkerSeq), file.path(opt$outputDir, opt$buildStdFragments_outputDir, 'stdFragments.tsv.gz'))

updateLog('buildStdFragments completed.')

q(save = 'no', status = 0, runLast = FALSE) 
