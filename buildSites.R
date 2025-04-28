#!/usr/bin/Rscript

# AAVengeR/buildSites.R
# John K. Everett, Ph.D.
# 
# This script accepts input from the buildStdFragments and assembles fragments 
# into integration events. Dual detections can be called if incoming fragments
# have U3 and U5 flags set.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(RMariaDB))
suppressPackageStartupMessages(library(Biostrings))

# Read in the configuration file and perform basic sanity checks.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) stop('Expected at least one command line argument')
source(file.path(yaml::read_yaml(args[1])$softwareDir, 'lib.R'))
opt <- startModule(args)

createOuputDir()
dir.create(file.path(opt$outputDir, opt$buildSites_outputDir), showWarnings = FALSE)
dir.create(file.path(opt$outputDir, opt$buildSites_outputDir, 'tmp'), showWarnings = FALSE)

# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, opt$buildSites_outputDir, 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)
write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)

quitOnErorr <- function(msg){
  if(opt$core_createFauxSiteDoneFiles) core_createFauxSiteDoneFiles()
  updateLog(msg)
  updateLog(paste0('See log for more details: ', opt$defaultLogFile))
  updateMasterLog()
  q(save = 'no', status = 1, runLast = FALSE) 
}

# Read in Standardized fragments.
updateLog('Reading standardized fragment data.')

frags <- readRDS(file.path(opt$outputDir, opt$buildSites_inputFile)) %>% dplyr::select(-readIDs) 


# Patch for chromosomes that contain dots which conflict with leaderSeq identifiers. 
#--------------------------
if(any(grepl('\\.', frags$chromosome))){
  o <- strsplit(frags$posid, '[\\+\\-]')
  s <- stringr::str_extract(frags$posid, '[\\+\\-]')
  frags$chromosome <- gsub('\\.', '&', unlist(lapply(o, '[', 1)))
  frags$posid <- paste0(gsub('\\.', '&', unlist(lapply(o, '[', 1))), s, unlist(lapply(o, '[', 2)))
}


# Create an identifier for counting unique widths for dual detections.
frags$fragWidths <- paste0(frags$strand, ':', abs(frags$fragEnd - frags$fragStart) + 1)

frags$fragWidths <- paste0(frags$strand, ':', abs(frags$fragEnd - frags$fragStart))

incomingSamples <- unique(paste0(frags$trial, '~', frags$subject, '~', frags$sample))

samples <- distinct(tibble(trial = frags$trial, subject = frags$subject, sample = frags$sample, replicate = frags$replicate, flags = frags$flags))

# Process fragments as integrase u5/u3 sites if all samples have a u5/u3 sample flag.
if('IN_u3' %in% frags$flags | 'IN_u5' %in% frags$flags){
  
  updateLog('Processing IN_u3 / IN_u5 flags.')
  
  if(opt$buildSites_enableDualDetection){
    updateLog('Processing dual detections.')
    
    processed_fragments <- frags[1,]
    processed_fragments[1,] <- NA
    
    # Dual detection - look for samples with both u3 and u5 entries.
    invisible(lapply(split(samples, paste(samples$trial, samples$subject, samples$sample)), function(x){
    
      if('IN_u5' %in% x$flags & 'IN_u3' %in% x$flags){
        # Isolate U3 and U5 fragments for this sample.
        u3.frags <- subset(frags, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1] & flags == 'IN_u3')
        u5.frags <- subset(frags, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1] & flags == 'IN_u5')
      
        # Remove leader sequence groupings to search for neighboring sites.
        u3.frags$posid2 <- sub('\\.\\d+$', '', u3.frags$posid)
        u5.frags$posid2 <- sub('\\.\\d+$', '', u5.frags$posid)
    
        # Cycle through u3 position ids.
        counter <- 1
        invisible(lapply(unique(u3.frags$posid), function(u3_posid){
          
          # Create alternative sites +/- a couple NT this u3 site.
          a <- sub('\\.\\d+$', '', u3_posid)
          o <- unlist(strsplit(a, '[\\+\\-]'))
          strand <- stringr::str_extract(a, '[\\+\\-]')
          alts <- paste0(o[1], ifelse(strand == '+', '-', '+'), (as.integer(o[2]) - opt$buildSites_dualDetectWidth):(as.integer(o[2]) + opt$buildSites_dualDetectWidth))
          
          # Search u5 fragments for alternative sites.
          # z may have multiple rows, one for each associated fragment.
          # More than one position id would single more than one closely spaced site was included (need to patch).
          z <- subset(u5.frags, posid2 %in% alts)
        
          # Are there fragments associated with possible U5 sites?
          if(nrow(z) > 0){
          
            # Retrieve fragments for both sites.
            f1 <- subset(u3.frags, posid2 == a & ! fragID %in% processed_fragments$fragID)        # u3
            f2 <- subset(u5.frags, posid2 %in% alts & ! fragID %in% processed_fragments$fragID)   # u5 - assuming only one alt site was found.
          
            if(nrow(f1) == 0 | nrow(f2) == 0) return()
       
            updateLog(paste0(counter, '. Processing U3 posid ', a, ' as a dual detection with ', nrow(f2), ' U5 fragments.'))
            counter <<- counter + 1
          
            # Records processed u5 fragments 
            i <- which(frags$fragID %in% c(f1$fragID, f2$fragID))
            processed_fragments <<- bind_rows(processed_fragments, frags[i,])
          
            i <- which(frags$fragID %in% c(f1$fragID, f2$fragID) & frags$strand == '+')
            frags[i,]$posid <<- unlist(lapply(strsplit(frags[i,]$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '+', as.integer(x[2]) + opt$buildSites_integraseCorrectionDist, '.', x[3])))
          
            i <- which(frags$fragID %in% c(f1$fragID, f2$fragID) & frags$strand == '-')
            frags[i,]$posid <<- unlist(lapply(strsplit(frags[i,]$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '-', as.integer(x[2]) - opt$buildSites_integraseCorrectionDist, '.', x[3])))
          
            i <- which(frags$fragID %in% c(f1$fragID, f2$fragID))
            frags[i,]$repLeaderSeq <<- paste0(names(sort(table(f1$repLeaderSeq), decreasing = TRUE))[1], '/', names(sort(table(f2$repLeaderSeq), decreasing = TRUE))[1])

            frags[i,]$flags <<- 'dual detect'
            
            chr <- stringr::str_extract(a, 'chr[XY\\d]+')
            pos <- names(sort(table(sub('[\\+\\-]', '', stringr::str_extract(frags[i,]$posid, '[\\+\\-]\\d+'))), decreasing = TRUE))[1]
          
            if(strand == '-'){
              frags[i,]$strand <<- '+'                          # Set the u3 frag strands to positive to reflect correct orientation. U5 posid already '+'.
              frags[i,]$posid  <<- paste0(chr, '+', pos, '.0')  # Set the u3 frag posids to the u5 posid which is positive causing its fragments to merge with U3 fragments.
            } else {
              frags[i,]$strand <<- '-'                          # Set the u3 frag strands to negative to reflect reverse orientation. U5 posid already '-'.
              frags[i,]$posid  <<- paste0(chr, '-', pos, '.0')  # Set the u3 frag posids to the u5 posid which is negative causing its fragments to merge with U3 fragments.
            }
          }
        }))
      }
    }))
  }
  
  frags <- bind_rows(lapply(split(frags, paste(frags$trial, frags$subject, frags$sample)), function(x){

    a <- subset(frags, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1] & flags == 'dual detect')
    b <- subset(frags, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1] & flags != 'dual detect')
    
    if(nrow(b)){
      
      # Shift positions to reflect duplication caused by integrase.
      b1 <- subset(b, strand == '+')
      if(nrow(b1) > 0) b1$posid <- unlist(lapply(strsplit(b1$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '+', as.integer(x[2]) + opt$buildSites_integraseCorrectionDist, '.', x[3])))

      b2 <- subset(b, strand == '-')
      if(nrow(b2) > 0) b2$posid <- unlist(lapply(strsplit(b2$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '-', as.integer(x[2]) - opt$buildSites_integraseCorrectionDist, '.', x[3])))

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

updateLog('Building replicate level integration sites.')

frags$replicate <- as.integer(frags$replicate)

maxRepNum <- max(frags$replicate)

# Dual detentions are assembled on the sample level.
# Move dual detentions to a faux replicate number of 0.
if(any(frags$flags == 'dual detect')) frags[frags$flags == 'dual detect',]$replicate <- 0

minRepNum <- min(frags$replicate)

consensusLeaderSeq <- function(x){
  tab <- dplyr::group_by(x, repLeaderSeq) %>% 
    dplyr::summarise(nWidths = n_distinct(fragWidths), nReads = sum(reads)) %>% 
    dplyr::ungroup() %>%
    dplyr::arrange(desc(nWidths), desc(nReads))
  
  as.character(tab[1, 'repLeaderSeq'])
}


consensusAnchorSeq <- function(x){
  tab <- dplyr::group_by(x, anchorReadSeq) %>% 
         dplyr::summarise(nWidths = n_distinct(fragWidths), nReads = sum(reads)) %>% 
         dplyr::ungroup() %>%
         dplyr::arrange(desc(nWidths), desc(nReads))
  
  as.character(tab[1, 'anchorReadSeq'])
}
  

clusterConsensusSeqs <- function(x){
  if(n_distinct(x$repLeaderSeq) > 1){
    seqs <- DNAStringSet(unique(x$repLeaderSeq))
    names(seqs) <- paste0('s', 1:length(seqs))
    o <- CD_HIT_clusters(seqs, opt$outputDir, opt$buildStdFragments_remnantClusterParams)
    o <- stringr::str_extract_all(o, '>[^\\.]+')
    return(length(o))
  } else {
    return(1)
  }
}


frags <- group_by(frags, trial, subject, sample, posid) %>% 
  mutate(g = cur_group_id()) %>%
  ungroup()


# Determine the best repLeaderSeq for each site.

replicateRepLeaderSeqTable <- buildRepLeaderSeqTable(frags, collapseReplicates = FALSE)
collapsedRepLeaderSeqTable <- buildRepLeaderSeqTable(frags, collapseReplicates = TRUE)

sites <- bind_rows(lapply(split(frags, frags$g), function(x){

  r <- bind_cols(lapply(minRepNum:maxRepNum, function(r){
         b <- tibble(rUMIs = NA, sonicLengths = NA, reads = NA, repLeaderSeq = NA)
         o <- x[x$replicate == r,]
    
         if(nrow(o) >= 1){
           repLeaderSeq <- subset(replicateRepLeaderSeqTable, trial == o$trial[1] & subject == o$subject[1] & sample == o$sample[1] & posid == o$posid[1] & replicate == o$replicate[1])$repLeaderSeq
           
           b$rUMIs <- ifelse(opt$processAdriftReadLinkerUMIs, n_distinct(unlist(o$rUMI_list)), NA)
           b$sonicLengths <- n_distinct(o$fragWidths)
           b$reads <- sum(o$reads)
           b$repLeaderSeq <- repLeaderSeq
         } 
    
         names(b) <- paste0('rep', r, '-', names(b))
         b
       }))
  
  repLeaderSeq <- subset(collapsedRepLeaderSeqTable, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1] & posid == x$posid[1] & replicate == '*')$repLeaderSeq
  
  bind_cols(tibble(trial = x$trial[1], subject = x$subject[1], sample = x$sample[1], 
                   refGenome = x$refGenome[1], 
                   posid = x$posid[1],
                   rUMIs = ifelse(opt$processAdriftReadLinkerUMIs, n_distinct(unlist(x$rUMI_list)), NA),
                   anchorReadClusters = ifelse(any(x$anchorReadCluster), TRUE, NA),
                   sonicLengths = ifelse(opt$buildSites_sumSonicBreaksWithin == 'replicates', 
                                         sum(r[, grepl('sonicLengths', names(r))], na.rm = TRUE),
                                         n_distinct(x$fragWidths)),
                   reads = sum(x$reads),
                   repLeaderSeq = repLeaderSeq, # subset(leaderSeqReps, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1] & posid == x$posid[1])$repLeaderSeq,
                   repLeaderSeqClusters = clusterConsensusSeqs(x),
                   nRepsObs = sum(! is.na(unlist(r[, which(grepl('reads', names(r)))]))),
                   anchorReadConsensus = consensusAnchorSeq(x),
                   flags = x$flags[1],
                   vector = x$vectorFastaFile[1]), r)
})) %>% arrange(desc(sonicLengths))

s <- unique(paste0(sites$trial, '~', sites$subject, '~', sites$sample))
if(any(! incomingSamples %in% s) & opt$core_createFauxSiteDoneFiles) core_createFauxSiteDoneFiles()

# Set nRepsObs to NA for dual detections since these have values of 1 after 
# moving dual detections to rep-0.

sites[sites$flags == 'dual detect',]$nRepsObs <- NA

sites <- group_by(sites, trial, subject, sample) %>%
         mutate(sampleAbund = sum(sonicLengths)) %>%
         ungroup() %>%
         group_by(posid) %>%
         mutate(percentSampleRelAbund = round((sonicLengths/sampleAbund[1]) * 100, 2), .after = 'nRepsObs') %>%
         ungroup() %>%
         select(-sampleAbund)

# Undo chromosome dot patch.
if(any(grepl('\\&', frags$chromosome))) sites$posid <- unname(gsub('&', '.', sites$posid))

# Hide rUMIs unless requested.
if(! opt$processAdriftReadLinkerUMIs) sites <- select(sites, -rUMIs)

# Save outputs.
saveRDS(sites, file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.rds'), compress = opt$compressDataFiles)
openxlsx::write.xlsx(sites, file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.xlsx'))
readr::write_tsv(sites, file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.tsv.gz'))
write(date(), file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.done'))

if(tolower(opt$database_configGroup) != 'none') uploadSitesToDB(sites)

updateLog('buildSites completed.')
updateMasterLog()

q(save = 'no', status = 0, runLast = FALSE) 
