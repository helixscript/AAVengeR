library(dplyr)
library(lubridate)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')

opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))
setMissingOptions()
setOptimalParameters()

dir.create(file.path(opt$outputDir, opt$buildSites_outputDir))
dir.create(file.path(opt$outputDir, opt$buildSites_outputDir, 'tmp'))

opt$defaultLogFile <- file.path(opt$outputDir, opt$buildSites_outputDir, 'log')

# Read in Standardized fragments.
updateLog('Reading standardized fragment data.')
frags <- readRDS(file.path(opt$outputDir, opt$buildSites_inputFile)) %>% dplyr::select(-readIDs) 

# Create an identifier for counting unique widths for dual detections.
frags$fragWidths <- paste0(frags$strand, ':', abs(frags$fragEnd - frags$fragStart) + 1)

incomingSamples <- unique(paste0(frags$trial, '~', frags$subject, '~', frags$sample))

samples <- distinct(tibble(trial = frags$trial, subject = frags$subject, sample = frags$sample, replicate = frags$replicate, flags = frags$flags))

# Process fragments as integrase u5/u3 sites if all samples have a u5/u3 sample flag.
if('IN_u3' %in% frags$flags | 'IN_u5' %in% frags$flags){
  
  updateLog('Processing IN_u3 / IN_u5 flags.')
  
  if(opt$buildSites_enableDualDetection){
    updateLog('Processing dual detections.')
    
    processed_fragments <- tibble()
  
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
            frags[i,]$posid <<- unlist(lapply(strsplit(frags[i,]$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '+', as.integer(x[2])+2, '.', x[3])))
          
            i <- which(frags$fragID %in% c(f1$fragID, f2$fragID) & frags$strand == '-')
            frags[i,]$posid <<- unlist(lapply(strsplit(frags[i,]$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '-', as.integer(x[2])-2, '.', x[3])))
          
            i <- which(frags$fragID %in% c(f1$fragID, f2$fragID))
            frags[i,]$repLeaderSeq <<- paste0(names(sort(table(f1$repLeaderSeq), decreasing = TRUE))[1], '/', names(sort(table(f2$repLeaderSeq), decreasing = TRUE))[1])

            frags[i,]$flags <<- 'dual detect'
          
            chr <- stringr::str_extract(a, 'chr[XY\\d]+')
            pos <- names(sort(table(sub('[\\+\\-]', '', stringr::str_extract(frags[i,]$posid, '[\\+\\-]\\d+'))), decreasing = TRUE))[1]
          
            if(strand == '-'){
              frags[i,]$strand <<- '+'                    # Set the u3 frag strands to positive to reflect correct orientation. U5 posid already '+'.
              frags[i,]$posid  <<- paste0(chr, '+', pos, '.0')  # Set the u3 frag posids to the u5 posid which is positive causing its fragments to merge with U3 fragments.
            } else {
              frags[i,]$strand <<- '-'                    # Set the u3 frag strands to negative to reflect reverse orientation. U5 posid already '-'.
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

updateLog('Building replicate level integration sites.')

frags$replicate <- as.integer(frags$replicate)

# Dual detentions are assembled on the sample level.
# Move dual detentions to a faux replicate number of 0.

frags[frags$flags == 'dual detect',]$replicate <- 0

buildConsensusSeq <- function(x){
  tab <- group_by(x, repLeaderSeq) %>% 
    summarise(nWidths = n_distinct(fragWidths), nReads = sum(reads)) %>% 
    ungroup() %>%
    arrange(desc(nWidths), desc(nReads))
  
  as.character(tab[1, 'repLeaderSeq'])
}

frags <- group_by(frags, trial, subject, sample, posid) %>% 
  mutate(g = cur_group_id()) %>%
  ungroup()

frags$fragWidths <- paste0(frags$strand, ':', abs(frags$fragEnd - frags$fragStart))

sites <- bind_rows(lapply(split(frags, frags$g), function(x){

  r <- bind_cols(lapply(0:max(frags$replicate), function(r){
         b <- tibble(rUMIs = NA, fUMIs = NA, sonicLengths = NA, reads = NA, repLeaderSeq = NA)
         o <- x[x$replicate == r,]
    
         if(nrow(o) > 0){
           b$rUMIs <- ifelse(opt$processAdriftReadLinkerUMIs, n_distinct(unlist(o$rUMI_list)), NA)
           b$fUMIs <- ifelse(opt$processAdriftReadLinkerUMIs, n_distinct(unlist(o$fUMI_list)), NA)
           b$sonicLengths <- n_distinct(o$fragWidths)
           b$reads <- sum(o$reads)
           b$repLeaderSeq <- buildConsensusSeq(o)
         } 
    
         names(b) <- paste0('rep', r, '-', names(b))
         b
       }))
  
  bind_cols(tibble(trial = x$trial[1], subject = x$subject[1], sample = x$sample[1],
                   rUMIs = ifelse(opt$processAdriftReadLinkerUMIs, n_distinct(unlist(x$rUMI_list)), NA),
                   fUMIs = ifelse(opt$processAdriftReadLinkerUMIs, n_distinct(unlist(x$fUMI_list)), NA),
                   sonicLengths = ifelse(opt$buildSites_sumSonicBreaksWithin == 'replicates', 
                                         sum(r[, grepl('sonicLengths', names(r))], na.rm = TRUE),
                                         n_distinct(x$fragWidths)),
                   reads = sum(x$reads),
                   repLeaderSeq = buildConsensusSeq(x),
                   nRepsObs = sum(! is.na(unlist(r[,which(grepl('reads', names(r)))]))),
                   flags = x$flags[1],
                   vector = x$vectorFastaFile[1]), r)
})) %>% arrange(desc(sonicLengths))

s <- unique(paste0(sites$trial, '~', sites$subject, '~', sites$sample))
if(any(! incomingSamples %in% s) & opt$core_createFauxSiteDoneFiles) core_createFauxSiteDoneFiles()

saveRDS(sites, file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.rds'), compress = opt$compressDataFiles)
openxlsx::write.xlsx(sites, file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.xlsx'))
readr::write_tsv(tbl2, file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.tsv.gz'))
write(date(), file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.done'))

if('databaseGroup' %in% names(opt)){
  library(RMariaDB)
  
  dbConn <- tryCatch({
    dbConnect(RMariaDB::MariaDB(), group = opt$databaseGroup)
  },
  error=function(cond) {
    updateLog('Error - could not connect to the database.')
    q(save = 'no', status = 1, runLast = FALSE) 
  })
  
  invisible(lapply(split(sites, paste(sites$trial, sites$subject, sites$sample, sites$replicate, sites$refGenome)), function(x){
    
    dbExecute(dbConn, paste0("delete from sites where trial='", x$trial[1], "' and subject='", x$subject[1],
                          "' and sample='", x$sample[1], "' and refGenome='", x$refGenome[1], "' and replicate='", x$replicate[1], "'"))
    
    o <- unlist(lapply(1:nrow(x), function(n){
           x <- x[n,]  
    
           r <- dbExecute(dbConn,
                          "insert into sites values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                           params = list(x$trial, x$subject, x$sample, x$refGenome,
                                         x$vectorFastaFile, x$flags, x$posid, x$UMIs, x$sonicLengths, x$reads, x$repLeaderSeq, x$replicate))
    }))
    
    if(any(o == 0)){
      updateLog(paste0('Error -- could not upload sites data for ', x$sample[1], ' to the database.'))
      q(save = 'no', status = 1, runLast = FALSE)
    } else {
      updateLog(paste0('Uploaded sites data for ', x$sample[1], ' to the database.'))
     }
  }))
  
  dbDisconnect(dbConn)
}

updateLog('buildSites completed.')

q(save = 'no', status = 0, runLast = FALSE) 
