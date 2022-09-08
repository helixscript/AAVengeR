library(dplyr)
library(lubridate)

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

# Dual detection
#-------------------------------------------------------------------------------
frags$id  <- paste(frags$trial, frags$subject, frags$sample, frags$replicate)
frags$id2 <- paste(frags$trial, frags$subject, frags$sample, frags$replicate, frags$posid)
frags$fragID <- paste0(frags$trial, ':', frags$subject, ':', frags$sample, ':', frags$replicate, ':', frags$posid, ':', frags$repLeaderSeq)

intSiteFlags <- tibble()

removeAltLeaderSeqs <- function(f1){
  i <- rep(TRUE, nrow(f1))
  r <- representativeSeq(f1$repLeaderSeq)
  if(r[[1]] > opt$buildStdFragments_maxLeaderSeqDiffScore){
    i <- as.vector(stringdist::stringdistmatrix(f1$repLeaderSeq, r[[2]]) / nchar(f1$repLeaderSeq) <= opt$buildStdFragments_maxLeaderSeqDiffScore)
    if(sum(i)/nrow(f1) >= opt$buildSites_assemblyConflictResolution){
      f1 <- f1[i,]
    } else {
      return()
    }
  }
  
  f1
}

# chr3+81746632 u3  chr3-81746636 u5 -> chr3-81746636   379 rows in frags.
invisible(lapply(split(samples, paste(samples$trial, samples$subject, samples$sample)), function(x){
  if(all('IN_u5' %in% x$flags & 'IN_u3' %in% x$flags)){
    u3 <- subset(x, flags == 'IN_u3')
    u5 <- subset(x, flags == 'IN_u5')
    
    # Isolate u3 and u5 sample-replicates 
    u3.ids <- paste(u3$trial, u3$subject, u3$sample, u3$replicate)
    u5.ids <- paste(u5$trial, u5$subject, u5$sample, u5$replicate)
    
    # Isolate corresponding u3 and u5 fragments.
    u3.frags <- subset(frags, id %in% u3.ids)
    u5.frags <- subset(frags, id %in% u5.ids)
    
    # If we find a dual detection, determine its orientation and move the frags 
    # from one strand to another. U3 neg means pos ort.
    lapply(unique(u3.frags$posid), function(a){
      
     
      
      o <- unlist(strsplit(a, '[\\+\\-]'))
      strand <- stringr::str_extract(a, '[\\+\\-]')
 
      alts <- paste0(o[1], ifelse(strand == '+', '-', '+'), (as.integer(o[2])-5):(as.integer(o[2])+5))
      z <- subset(u5.frags, posid %in% alts)
      if(nrow(z) > 0){
        # Just need to change the strand and position id of the records without the 
        # the correction orientation strand. Frags will be pointing away from others 
        # but will be tallied correctly. 
        f1 <- subset(u3.frags, posid == a)           # 379 rows
        f2 <- subset(u5.frags, posid %in% z$posid)   # 286 rows
        
        
        # Asuming dual detections have nothing to correct.
        # Add correction in next release...
        
        ### f1 <- removeAltLeaderSeqs(f1)
        ### f2 <- removeAltLeaderSeqs(f2)
        ### if(nrow(f1) == 0 | nrow(f2) == 0) return()
        
        # Assign the same leaderSeq sequence to u3 and u5 can be processed together
        # otherwise the LTR similarity test would fail. Consider reporting u3,u5 seqs for completeness. 
        n <- names(sort(table(c(frags[frags$id2 %in% f1$id2,]$repLeaderSeq, frags[frags$id2 %in% f2$id2,]$repLeaderSeq)), decreasing = TRUE)[1])
        frags[frags$id2 %in% f1$id2,]$repLeaderSeq <<- n
        frags[frags$id2 %in% f2$id2,]$repLeaderSeq <<- n
        
        
        intSiteFlags <<- bind_rows(intSiteFlags, 
                                   tibble(trial = f1$trial[1], subject = f1$subject[1], 
                                          sample = f1$sample[1], u3_posid = a, u5_posid = f2$posid[1],
                                          posid = ifelse(strand == '+', f2$posid[1], a), flag = 'dual detect'))
        
        
        
        if(strand == '+'){
          frags[frags$id2 %in% f1$id2,]$strand <<- '-' # Set the u3 frag strands to negative to reflect orientation.
          frags[frags$id2 %in% f1$id2,]$posid <<- f2$posid[1] # Set the u3 frag posids to the u5 posid which essentially removes the u3 records.
        } else {
          frags[frags$id2 %in% f1$id2,]$strand <<- '+' # set the u3 strand to positive to reflect orientation.
          frags[frags$id2 %in% f2$id2,]$posid <<- a    # Set the u5 frag posids to the u3 posid which essentially removes the u5 records.
        }
      }
    })
  } else {
    message(x$sample)
  }
}))

# Separate out duel detections since their orientation strand has been set.
a <- subset(frags, posid %in% intSiteFlags$posid)
b <- subset(frags, ! posid %in% intSiteFlags$posid)

# Add u3 / u5 flags to fragments.
samples$id <- paste(samples$trial, samples$subject, samples$sample, samples$replicate)
b <- left_join(b, select(samples, id, flags), by = 'id')

# Shift positions to relect duplication caused by integrase.

# Adjustments
# + strand alignment: +2
# - strand alignment: -2
# flip strand

b1 <- subset(b, strand == '+')
b1$posid <- unlist(lapply(strsplit(b1$posid, '\\+'), function(x) paste0(x[1], '+', as.integer(x[2])+2)))

b2 <- subset(b, strand == '-')
b2$posid <- unlist(lapply(strsplit(b2$posid, '\\+'), function(x) paste0(x[1], '+', as.integer(x[2])-2)))

b <- bind_rows(b1, b2)

updatePosIdStrand <- function(x, s){
  o <- unlist(strsplit(x, '[\\+\\-]'))
  paste0(o[1], s, o[2])
}

# Change strand to reflect orientation.  5072
b1 <- subset(b, strand == '+' & grepl('IN_u3', b$flags))
b2 <- subset(b, strand == '-' & grepl('IN_u3', b$flags))
b3 <- subset(b, strand == '+' & grepl('IN_u5', b$flags))
b4 <- subset(b, strand == '-' & grepl('IN_u5', b$flags))

if(nrow(b1) > 0) b1$posid <- sapply(b1$posid, updatePosIdStrand, '-')
if(nrow(b2) > 0) b2$posid <- sapply(b2$posid, updatePosIdStrand, '+')
if(nrow(b3) > 0) b3$posid <- sapply(b3$posid, updatePosIdStrand, '+')
if(nrow(b4) > 0) b4$posid <- sapply(b4$posid, updatePosIdStrand, '-')

b <- bind_rows(b1, b2, b3, b4)

frags <- bind_rows(a, b)

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

# chr3+81746632 u3 -> chr3-81746636 u5

sites <- bind_rows(lapply(o, function(x){
  message(counter, ' / ', total); counter <<- counter + 1

  if(opt$buildStdFragments_categorize_anchorRead_remnants) x$posid <- paste0(x$posid, '.', x$repLeaderSeqGroup)
  
  o <- calcAbunds(x); estAbund <- o[[1]]; molCodes <- o[[2]]
  
  if(nrow(x) == 1){
    
   # Just one fragment for this position.
    return(dplyr::mutate(x, fragments = nrow(x), fragmentWidths = n_distinct(x$fragWidth), sonicAbund = estAbund, UMIs = molCodes, fragmentsRemoved = 0) %>%
             dplyr::select(trial, subject, sample, replicate, posid, estAbund, UMIs, reads, fragmentsRemoved, repLeaderSeq))
    
  } else {
    
    #if(x$posid[1] == 'chr3-81746636') browser()
    
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

tbl2$id <- paste(tbl2$trial, tbl2$subject, tbl2$sample, tbl2$posid)
intSiteFlags$id <- paste(intSiteFlags$trial, intSiteFlags$subject, intSiteFlags$sample, intSiteFlags$posid)

tbl2 <- left_join(tbl2, select(intSiteFlags, id, flag), by = 'id')
tbl2 <- relocate(tbl2, flag, .before = fragments) 

saveRDS(select(tbl2, -id), file.path(opt$outputDir, opt$buildSites_outputDir, opt$buildSites_outputFile))
openxlsx::write.xlsx(arrange(select(tbl2, -id), desc(UMIs)), file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.xlsx'))

q(save = 'no', status = 0, runLast = FALSE) 
