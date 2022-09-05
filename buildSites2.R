library(dplyr)
library(sonicLength)
library(Biostrings)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$buildSites_outputDir))

# Read in Standardized fragments.
frags <- readRDS(file.path(opt$outputDir, opt$buildSites_inputFile))

if('buildSites_excludeSites' %in% names(opt)){
  p <- unlist(lapply(strsplit(unlist(strsplit(opt$buildSites_excludeSites, '\\|')), ','), function(x){
    if(length(x) != 3) return()
    subset(frags, seqnames == x[1] & start >= x[2] & start <= x[3])$posid
  }))
  
  frags <- subset(frags, ! posid %in% p)
}

samples <- loadSamples()

frags$fragWidth <- frags$fragEnd - frags$fragStart + 1
frags$replicate <- as.integer(frags$replicate)


calcAbunds <- function(x){
    if(n_distinct(x$fragWidth) >= 5){
      estAbund <- tryCatch(
        {
          o <- sonicLength::estAbund(locations = x$posid, lengths = x$fragWidth)
          as.integer(o$theta)
        },
        error=function(cond) {
          return(NA)
        }
      )    
    }else{
      estAbund <- n_distinct(x$fragWidth)
    }
  
  list(estAbund, n_distinct(x$randomLinkerSeq.adriftReads))
}


# Build replicate level sites.
if(! 'repLeaderSeqGroup' %in% names(frags)) frags$repLeaderSeqGroup <- 'X'
o <- split(frags, paste(frags$trial, frags$subject, frags$sample, frags$replicate, frags$repLeaderSeqGroup, frags$posid))
  
counter <- 1
total <- length(o)

sites <- bind_rows(lapply(o, function(x){
  message(counter, ' / ', total); counter <<- counter + 1

  if(opt$buildStdFragments_categorize_anchorRead_remnants) x$posid <- paste0(x$posid, '.', x$repLeaderSeqGroup)
  
  o <- calcAbunds(x); estAbund <- o[[1]]; molCodes <- o[[2]]
  
  if(nrow(x) == 1){
    
   # Just one fragment for this position.
    return(dplyr::mutate(x, fragments = nrow(x), fragmentWidths = n_distinct(x$fragWidth), sonicAbund = estAbund, UMIs = molCodes, fragmentsRemoved = 0) %>%
             dplyr::select(trial, subject, sample, replicate, posid, estAbund, UMIs, reads, fragmentsRemoved, repLeaderSeq))
    
  } else {
    i <- rep(TRUE, nrow(x))
    r <- representativeSeq(x$repLeaderSeq)
    
    # Exclude fragments that have distinctly different leader sequences compared to the consensus sequence.

    if(r[[1]] > opt$buildStdFragments_maxLeaderSeqDiffScore){
      i <- as.vector(stringdist::stringdistmatrix(x$repLeaderSeq, r[[2]]) / nchar(x$repLeaderSeq) <= opt$buildStdFragments_maxLeaderSeqDiffScore)
     
      # Is there a majority of fragments that are similar to the consensus sequences?
      # If there is we salvage those fragments otherwise we discard this site.
      if(sum(i)/nrow(x) >= opt$buildSites_assemblyConflictResolution){
        x <- x[i,]
        r <- representativeSeq(x$repLeaderSeq)
      } else {
        return(tibble())
      }
      
      # Recalculate abundances since fragments were removed.
      o <- calcAbunds(x); estAbund <- o[[1]]; molCodes <- o[[2]]
    }
      
    return(dplyr::mutate(x, fragments = nrow(x), fragmentWidths = n_distinct(x$fragWidth), sonicAbund = estAbund, UMIs = molCodes, reads = sum(reads), repLeaderSeq = r[[2]], fragmentsRemoved = sum(!i)) %>%
           dplyr::select(trial, subject, sample, replicate, posid, fragments, fragmentWidths, sonicAbund, UMIs, reads, fragmentsRemoved, repLeaderSeq) %>%
           dplyr::slice(1))
  }
}))
  

# Create a wide view of the replicate level sites and create NA cells 
# for replicates where specific sites were not found.
  
minReplicate <- min(sites$replicate)
maxReplicate <- max(sites$replicate)

tbl1 <- bind_rows(lapply(split(sites, paste(sites$trial, sites$subject, sites$sample, sites$posid)), function(x){ 
         o <- bind_cols(lapply(minReplicate:maxReplicate, function(r){
           o <- subset(x, replicate == r)
        
           if(nrow(o) == 1){
             t <- tibble(x1 = o$fragments, x2 = o$fragmentWidths, x3 = o$sonicAbund, x4 = o$UMIs, x5 = o$reads, x6 = o$fragmentsRemoved, x7 = o$repLeaderSeq)
           } else if(nrow(o) > 1){
             stop('Row error 1')  
           } else {
             t <- tibble(x1 = NA, x2 = NA, x3 = NA, x4 = NA, x5 = NA, x6 = NA, x7 = NA)
           }

           names(t) <- c(paste0('rep', r, '-fragments'), 
                         paste0('rep', r, '-fragmentWidths'), 
                         paste0('rep', r, '-sonicAbund'),
                         paste0('rep', r, '-UMIs'), 
                         paste0('rep', r, '-reads'), 
                         paste0('rep', r, '-fragsRemoved'),
                         paste0('rep', r, '-repLeaderSeq'))
          t
         }))
         
      bind_cols(tibble(trial = x$trial[1], subject = x$subject[1], sample = x$sample[1], posid = x$posid[1]), o)
}))


tbl2 <- bind_rows(lapply(1:nrow(tbl1), function(x){
  x <- tbl1[x,]

  f <- subset(frags, trial == x$trial & subject == x$subject & sample == x$sample & posid == x$posid)
  
  o <- calcAbunds(f); estAbund <- o[[1]]; molCodes <- o[[2]]
  
  if(nrow(f) == 1){
    
    # Just one fragment for this position.
    k <- tibble(fragments = nrow(f), fragmentWidths = n_distinct(f$fragWidth), sonicAbund = estAbund, UMIs = molCodes, reads = f$reads, fragmentsRemoved = 0, repLeaderSeq = f$repLeaderSeq) 
    
  } else {
    i <- rep(TRUE, nrow(f))
    r <- representativeSeq(f$repLeaderSeq)
    
    # Exclude fragments that have distinctly different leader sequences compared to the consensus sequence.
    
    if(r[[1]] > opt$buildStdFragments_maxLeaderSeqDiffScore){
      i <- as.vector(stringdist::stringdistmatrix(f$repLeaderSeq, r[[2]]) / nchar(f$repLeaderSeq) <= opt$buildStdFragments_maxLeaderSeqDiffScore)
      
      # Is there a majority of fragments that are similar to the consensus sequences?
      # If there is we salvage those fragments otherwise we discard this site.
      if(sum(i)/nrow(f) >= opt$buildSites_assemblyConflictResolution){
        f <- f[i,]
        r <- representativeSeq(f$repLeaderSeq)
      } else {
        return(tibble())
      }
      
      # Recalculate abundances since fragments were removed.
      o <- calcAbunds(f); estAbund <- o[[1]]; molCodes <- o[[2]]
    }
    
    k <- tibble(fragments = nrow(f), fragmentWidths = n_distinct(f$fragWidth), sonicAbund = estAbund, UMIs = molCodes, reads = sum(f$reads), fragmentsRemoved = sum(!i), repLeaderSeq = r[[2]])
  }
  
  k$nRepsObs <- sum(! is.na(unlist(x[, which(grepl('\\-repLeaderSeq', names(x)))])))
  k$repLeaderSeqDists <- paste0(stringdist::stringdist(k$repLeaderSeq, unlist(x[, which(grepl('\\-repLeaderSeq', names(x)))])), collapse = ', ')
  bind_cols(x[,1:4], k, x[,5:length(x)])
}))

saveRDS(tbl2, file.path(opt$outputDir, opt$buildSites_outputDir, opt$buildSites_outputFile))
openxlsx::write.xlsx(arrange(tbl2, desc(UMIs)), file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.xlsx'))

q(save = 'no', status = 0, runLast = FALSE) 
