library(dplyr)
library(lubridate)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$buildSites_outputDir))

# Read in Standardized fragments.
write(c(paste(now(), '   Reading standardized fragment data.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
frags <- readRDS(file.path(opt$outputDir, opt$buildSites_inputFile))

samples <- distinct(tibble(trial = frags$trial, subject = frags$subject, sample = frags$sample, replicate = frags$replicate, flags = frags$flags))

# Process fragments as integrase u5/u3 sites if all samples have a u5/u3 sample flag.
if('IN_u3' %in% frags$flags | 'IN_u5' %in% frags$flags){
  
  write(c(paste(now(), '   Processing IN_u3 / IN_u5 flags.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
  frags$id  <- paste(frags$trial, frags$subject, frags$sample, frags$replicate)
  frags$id2 <- paste(frags$trial, frags$subject, frags$sample, frags$replicate, frags$posid)
  
  intSiteFlags <- tibble()
  
  # Dual detection - look for samples with both u3 and u5 entries.
  invisible(lapply(split(samples, paste(samples$trial, samples$subject, samples$sample)), function(x){
    
    if('IN_u5' %in% x$flags & 'IN_u3' %in% x$flags){
      
      # Subset u3 / u5 isolates.
      u3 <- subset(x, flags == 'IN_u3')
      u5 <- subset(x, flags == 'IN_u5')
    
      u3.ids <- paste(u3$trial, u3$subject, u3$sample, u3$replicate)
      u5.ids <- paste(u5$trial, u5$subject, u5$sample, u5$replicate)
    
      # Isolate corresponding u3 and u5 fragments.
      u3.frags <- subset(frags, id %in% u3.ids)
      u5.frags <- subset(frags, id %in% u5.ids)
      
      # Remove leader sequence groupings to search for neighboring sites.
      u3.frags$posid2 <- sub('\\.\\d+$', '', u3.frags$posid)
      u5.frags$posid2 <- sub('\\.\\d+$', '', u5.frags$posid)
    
      
      # If we find a dual detection, determine its orientation and move the frags 
      # from one strand to another. U3 neg means pos ort.
      invisible(lapply(unique(u3.frags$posid), function(u3_posid){
        
        a <- sub('\\.\\d+$', '', u3_posid)
        o <- unlist(strsplit(a, '[\\+\\-]'))
        strand <- stringr::str_extract(a, '[\\+\\-]')
 
        # Create alternative sites +/- 5 this u3 site.
        alts <- paste0(o[1], ifelse(strand == '+', '-', '+'), (as.integer(o[2])-6):(as.integer(o[2])+6))
        
        # Search u5 fragments for alternative sites.
        # z may have multiple rows, one for each associated fragment.
        # More than one position id would single more than one closely spaced site was included (need to patch).
        z <- subset(u5.frags, posid2 %in% alts)
        

        if(nrow(z) > 0){
          # Just need to change the strand and position id of the records without the 
          # the correction orientation strand. Frags will be pointing away from others 
          # but will be tallied correctly. 
          f1 <- subset(u3.frags, posid2 == a)            # u3
          f2 <- subset(u5.frags, posid2 %in% z$posid2)   # u5
        
          # Assign both most common u3 and u5 leader sequences to each proximal site.
          f1.seq <- names(sort(table(frags[frags$id2 %in% f1$id2,]$repLeaderSeq), decreasing = TRUE))[1]
          f2.seq <- names(sort(table(frags[frags$id2 %in% f2$id2,]$repLeaderSeq), decreasing = TRUE))[1]                             
                                       
          frags[frags$id2 %in% f1$id2,]$repLeaderSeq <<- paste0(f1.seq, '/', f2.seq)
          frags[frags$id2 %in% f2$id2,]$repLeaderSeq <<- paste0(f1.seq, '/', f2.seq)
        
          intSiteFlags <<- bind_rows(intSiteFlags, 
                                     tibble(trial = f1$trial[1], subject = f1$subject[1], 
                                            sample = f1$sample[1], u3_posid = u3_posid, u5_posid = f2$posid[1],
                                            posid = f2$posid[1]))
          
          # [u3][R][u5]------------[u3][R][u5]
          
          if(strand == '-'){
            frags[frags$id2 %in% f1$id2,]$strand <<- '+'         # Set the u3 frag strands to positive to reflect correct orientation. U5 posid already '+'.
            frags[frags$id2 %in% f1$id2,]$posid  <<- f2$posid[1] # Set the u3 frag posids to the u5 posid which is positive causing its fragments to merge with U3 fragments.
          } else {
            frags[frags$id2 %in% f1$id2,]$strand <<- '-'         # Set the u3 frag strands to negative to reflect reverse orientation. U5 posid already '-'.
            frags[frags$id2 %in% f1$id2,]$posid  <<- f2$posid[1] # Set the u3 frag posids to the u5 posid which is negative causing its fragments to merge with U3 fragments.
          }
        }
      }))
    }
  }))
  
  frags <- bind_rows(lapply(split(frags, paste(frags$trial, frags$subject, frags$sample)), function(x){

    dualDetect <- tibble()
    if(nrow(intSiteFlags) > 0) dualDetect <- subset(intSiteFlags, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1])
    
    a <- subset(frags, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1] & posid %in% dualDetect$posid)
    b <- subset(frags, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1] & ! posid %in% dualDetect$posid)
  
    if(nrow(a)){
      a1 <- subset(a, strand == '+')
      if(nrow(a1) > 0) a1$posid <- unlist(lapply(strsplit(a1$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '+', as.integer(x[2])+2, '.', x[3])))
  
      a2 <- subset(a, strand == '-')
      if(nrow(a2) > 0) a2$posid <- unlist(lapply(strsplit(a2$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '-', as.integer(x[2])-2, '.', x[3])))
  
      a <- bind_rows(a1, a2)
      a$flags <- 'dual detect'
    }
    
    # Shift positions to reflect duplication caused by integrase.
    if(nrow(b)){
      b1 <- subset(b, strand == '+')
      if(nrow(b1) > 0) b1$posid <- unlist(lapply(strsplit(b1$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '+', as.integer(x[2])+2, '.', x[3])))

      b2 <- subset(b, strand == '-')
      if(nrow(b2) > 0) b2$posid <- unlist(lapply(strsplit(b2$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '-', as.integer(x[2])-2, '.', x[3])))

      b <- bind_rows(b1, b2)

      updatePosIdStrand <- function(x, s){
        o <- unlist(strsplit(x, '[\\+\\-]'))
        paste0(o[1], s, o[2])
      }

      # Change strand to reflect orientation. 
      b1 <- subset(b, strand == '+' & grepl('IN_u3', b$flags))
      b2 <- subset(b, strand == '-' & grepl('IN_u3', b$flags))
      b3 <- subset(b, strand == '+' & grepl('IN_u5', b$flags))
      b4 <- subset(b, strand == '-' & grepl('IN_u5', b$flags))

      if(nrow(b1) > 0) b1$posid <- sapply(b1$posid, updatePosIdStrand, '-')
      if(nrow(b2) > 0) b2$posid <- sapply(b2$posid, updatePosIdStrand, '+')
      if(nrow(b3) > 0) b3$posid <- sapply(b3$posid, updatePosIdStrand, '+')
      if(nrow(b4) > 0) b4$posid <- sapply(b4$posid, updatePosIdStrand, '-')

      b <- bind_rows(b1, b2, b3, b4)
    }
      
    bind_rows(a, b)
  }))
}

frags$fragWidth <- frags$fragEnd - frags$fragStart + 1
frags$replicate <- as.integer(frags$replicate)


write(c(paste(now(), '   Building replicate level integration sites.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

o <- split(frags, paste(frags$trial, frags$subject, frags$sample, frags$replicate, frags$posid))
counter <- 1
total <- length(o)

sites <- bind_rows(lapply(o, function(x){
  ### if(counter == 4797) browser()
  message(counter, ' / ', total); counter <<- counter + 1
  
  if(nrow(x) == 1){
    return(dplyr::mutate(x, fragments = n_distinct(x$randomLinkerSeq.adriftReads), fragmentWidths = n_distinct(x$fragWidth), maxLeaderSeqDist = 0) %>%
             dplyr::select(trial, subject, sample, replicate, refGenome, posid, flags, fragments, fragmentWidths, reads, maxLeaderSeqDist, repLeaderSeq))
  } else {
    
    s <- unlist(x$leaderSeqs)
    if(length(s) > opt$buildStdFragments_representativeSeqCalc_maxReads){
      set.seed(1)
      s <- sample(s, opt$buildStdFragments_representativeSeqCalc_maxReads)
    }
    
    r <- representativeSeq(s)
    
    return(dplyr::mutate(x, fragments = n_distinct(x$randomLinkerSeq.adriftReads), fragmentWidths = n_distinct(x$fragWidth), reads = sum(reads), repLeaderSeq = r[[2]], 
                         maxLeaderSeqDist = max(stringdist::stringdistmatrix(s))) %>%
           dplyr::select(trial, subject, sample, replicate, refGenome, posid, flags, fragments, fragmentWidths, reads, maxLeaderSeqDist, repLeaderSeq) %>%
           dplyr::slice(1))
  }
}))
  

# Create a wide view of the replicate level sites and create NA cells 
# for replicates where specific sites were not found.
  
minReplicate <- min(sites$replicate)
maxReplicate <- max(sites$replicate)

write(c(paste(now(), '   Creating a wide view of the replicate level integration site data.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

tbl1 <- bind_rows(lapply(split(sites, paste(sites$trial, sites$subject, sites$sample, sites$posid)), function(x){ 
         o <- bind_cols(lapply(minReplicate:maxReplicate, function(r){
           o <- subset(x, replicate == r)
        
           if(nrow(o) == 1){
             t <- tibble(x1 = o$fragments, x2 = o$fragmentWidths, x3 = o$reads, x4 = o$repLeaderSeq)
           } else if(nrow(o) > 1){
             stop('Row error 1')  
           } else {
             t <- tibble(x1 = NA, x2 = NA, x3 = NA, x4 = NA)
           }

           names(t) <- c(paste0('rep', r, '-fragments'), 
                         paste0('rep', r, '-fragmentWidths'), 
                         paste0('rep', r, '-reads'), 
                         paste0('rep', r, '-repLeaderSeq'))
          t
         }))
         
      bind_cols(tibble(trial = x$trial[1], subject = x$subject[1], sample = x$sample[1], refGenome = x$refGenome[1], posid = x$posid[1], flags = x$flags[1]), o)
}))


write(c(paste(now(), '   Buiding sample level integration sites.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

tbl2 <- bind_rows(lapply(1:nrow(tbl1), function(x){
  x <- tbl1[x,]

  f <- subset(frags, trial == x$trial & subject == x$subject & sample == x$sample & posid == x$posid)
  
  if(nrow(f) == 1){
    k <- tibble(fragments = n_distinct(f$randomLinkerSeq.adriftReads), fragmentWidths = n_distinct(f$fragWidth), reads = f$reads, 
                maxLeaderSeqDist = f$maxLeaderSeqDist, repLeaderSeq = f$repLeaderSeq) 
    
  } else {
    s <- unlist(f$leaderSeqs)
    if(length(s) > opt$buildStdFragments_representativeSeqCalc_maxReads){
      set.seed(1)
      s <- sample(s, opt$buildStdFragments_representativeSeqCalc_maxReads)
    }
    
    r <- representativeSeq(s)
    
    k <- tibble(fragments = n_distinct(f$randomLinkerSeq.adriftReads), fragmentWidths = n_distinct(f$fragWidth), reads = sum(f$reads), 
                maxLeaderSeqDist = max(stringdist::stringdistmatrix(s)), repLeaderSeq = r[[2]])
  }
  
  k$nRepsObs <- sum(! is.na(unlist(x[, which(grepl('\\-repLeaderSeq', names(x)))])))
  bind_cols(x[,1:5], k, x[,6:length(x)])
}))

if(any(is.infinite(tbl2$maxLeaderSeqDist))) tbl2[is.infinite(tbl2$maxLeaderSeqDist),]$maxLeaderSeqDist <- NA

saveRDS(tbl2, file.path(opt$outputDir, opt$buildSites_outputDir, opt$buildSites_outputFile))
openxlsx::write.xlsx(arrange(tbl2, desc(fragments)), file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.xlsx'))

q(save = 'no', status = 0, runLast = FALSE) 
