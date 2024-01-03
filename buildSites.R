library(dplyr)
library(lubridate)
library(parallel)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))
setMissingOptions()
setOptimalParameters()

dir.create(file.path(opt$outputDir, opt$buildSites_outputDir))
dir.create(file.path(opt$outputDir, opt$buildSites_outputDir, 'tmp'))

# Read in Standardized fragments.
write(c(paste(now(), '   Reading standardized fragment data.')), file = file.path(opt$outputDir, opt$buildSites_outputDir, 'log'), append = FALSE)
frags <- readRDS(file.path(opt$outputDir, opt$buildSites_inputFile)) %>% dplyr::select(-readIDs) 

incomingSamples <- unique(paste0(frags$trial, '~', frags$subject, '~', frags$sample))

samples <- distinct(tibble(trial = frags$trial, subject = frags$subject, sample = frags$sample, replicate = frags$replicate, flags = frags$flags))

# Process fragments as integrase u5/u3 sites if all samples have a u5/u3 sample flag.
if('IN_u3' %in% frags$flags | 'IN_u5' %in% frags$flags){
  
  write(c(paste(now(), '   Processing IN_u3 / IN_u5 flags.')), file = file.path(opt$outputDir, opt$buildSites_outputDir, 'log'), append = TRUE)
  
  frags$id <- paste(frags$trial, frags$subject, frags$sample, frags$replicate, frags$posid)
  
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
      invisible(lapply(unique(u3.frags$posid), function(u3_posid){
        
        # Create alternative sites +/- 5 this u3 site.
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
          f1 <- subset(u3.frags, posid2 == a)        # u3
          f2 <- subset(u5.frags, posid2 %in% alts)   # u5 - assuming only one alt site was found.
          
          i <- which(frags$fragID %in% c(f1$fragID, f2$fragID) & frags$strand == '+')
          frags[i,]$posid <<- unlist(lapply(strsplit(frags[i,]$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '+', as.integer(x[2])+2, '.', x[3])))
          
          i <- which(frags$fragID %in% c(f1$fragID, f2$fragID) & frags$strand == '-')
          frags[i,]$posid <<- unlist(lapply(strsplit(frags[i,]$posid, '[\\+\\-\\.]', perl = TRUE), function(x) paste0(x[1], '-', as.integer(x[2])-2, '.', x[3])))
          
          i <- which(frags$fragID %in% c(f1$fragID, f2$fragID))
          frags[i,]$repLeaderSeq <<- paste0(names(sort(table(f1$repLeaderSeq), decreasing = TRUE))[1], '/', names(sort(table(f2$repLeaderSeq), decreasing = TRUE))[1])

          frags[i,]$flags <<- 'dual detect'
          
          chr <- stringr::str_extract(a, 'chr[XY\\d+]')
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

write(c(paste(now(), '   Building replicate level integration sites.')), file = file.path(opt$outputDir, opt$buildSites_outputDir, 'log'), append = TRUE)

frags$fragWidth <- frags$fragEnd - frags$fragStart + 1
frags$replicate <- as.integer(frags$replicate)

# Create a replicate level / posid grouping id then create another grouping id for parallel processing.
frags <- group_by(frags, trial, subject, sample, replicate, posid) %>% 
         mutate(g = cur_group_id()) %>%
         ungroup()

cluster <- makeCluster(opt$buildSites_CPUs)
clusterExport(cluster, c('opt', 'frags'))

#sites <- bind_rows(parLapply(cluster, split(frags, frags$g), function(x){  
sites <- bind_rows(lapply(split(frags, frags$g), function(x){    
           library(dplyr)
           library(Biostrings)
           library(stringdist)
           source(file.path(opt$softwareDir, 'lib.R'))
  
           if(nrow(x) == 1){
              return(dplyr::mutate(x, rUMIs = n_distinct(unlist(x$rUMI_list)),  
                                      fUMIs = n_distinct(unlist(x$fUMI_list)),
                                      sonicLengths = n_distinct(x$fragWidth)) %>%
                     dplyr::select(trial, subject, sample, replicate, refGenome, posid, flags, rUMIs, fUMIs, sonicLengths, reads, repLeaderSeq, vectorFastaFile))
           } else {
              return(dplyr::mutate(x, 
                                   rUMIs = n_distinct(unlist(x$rUMI_list)),
                                   fUMIs = n_distinct(unlist(x$fUMI_list)),
                                   sonicLengths = n_distinct(x$fragWidth), 
                                   reads = sum(reads), 
                                   repLeaderSeq = repLeaderSeq[1]) %>%
                    dplyr::select(trial, subject, sample, replicate, refGenome, posid, flags, rUMIs, fUMIs, sonicLengths, reads, repLeaderSeq, vectorFastaFile) %>%
                    dplyr::slice(1))
           }
         })) 
  
# Create a wide view of the replicate level sites and create NA cells 
# for replicates where specific sites were not found.
  
minReplicate <- min(sites$replicate)
maxReplicate <- max(sites$replicate)

write(c(paste(now(), '   Creating a wide view of the replicate level integration site data.')), file = file.path(opt$outputDir, opt$buildSites_outputDir, 'log'), append = TRUE)

tbl1 <- bind_rows(lapply(split(sites, paste(sites$trial, sites$subject, sites$sample, sites$posid)), function(x){ 
          o <- bind_cols(lapply(minReplicate:maxReplicate, function(r){
          
                 o <- subset(x, replicate == r)
        
                 if(nrow(o) == 1){
                   t <- tibble(x1 = o$rUMIs, x2 = o$fUMIs, x3 = o$sonicLengths, x4 = o$reads, x5 = o$repLeaderSeq)
                 } else if(nrow(o) > 1){
                   stop('Row error 1')  
                 } else {
                   t <- tibble(x1 = NA, x2 = NA, x3 = NA, x4 = NA, x5 = NA)
                 }
           
                 names(t) <- c(paste0('rep', r, '-rUMIs'), 
                               paste0('rep', r, '-fUMIs'),
                               paste0('rep', r, '-sonicLengths'), 
                               paste0('rep', r, '-reads'), 
                               paste0('rep', r, '-repLeaderSeq'))
                t
              }))
         
          bind_cols(tibble(trial = x$trial[1], subject = x$subject[1], sample = x$sample[1], refGenome = x$refGenome[1], posid = x$posid[1], flags = x$flags[1], vectorFastaFile = x$vectorFastaFile[1]), o)
}))


buildConsensusSeq <- function(x){
  tab <- group_by(x, repLeaderSeq) %>% 
         summarise(nWidths = n_distinct(fragWidth), nReads = sum(reads)) %>% 
         ungroup() %>%
         arrange(desc(nWidths), desc(nReads))
  as.character(tab[1, 'repLeaderSeq'])
}

tbl2 <-  bind_rows(lapply(1:nrow(tbl1), function(n){
             x <- tbl1[n,]

             f <- subset(frags, trial == x$trial & subject == x$subject & sample == x$sample & posid == x$posid)

             k <- tibble(rUMIs = n_distinct(unlist(f$rUMI_list)),
                         fUMIs = n_distinct(unlist(f$fUMI_list)), 
                         sonicLengths = n_distinct(f$fragWidth), 
                         reads = sum(f$reads), 
                         repLeaderSeq = buildConsensusSeq(f))
     
             k$nRepsObs <- sum(! is.na(unlist(x[, which(grepl('\\-reads', names(x)))])))
        
             bind_cols(x[,1:5], k, x[,6:length(x)])
          })) %>% arrange(desc(sonicLengths))

s <- unique(paste0(tbl2$trial, '~', tbl2$subject, '~', tbl2$sample))
if(any(! incomingSamples %in% s) & opt$core_createFauxSiteDoneFiles) core_createFauxSiteDoneFiles()

if(! opt$processAdriftReadLinkerUMIs) tbl2$UMIs <- NA

saveRDS(tbl2, file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.rds'), compress = opt$compressDataFiles)
openxlsx::write.xlsx(tbl2, file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.xlsx'))
readr::write_tsv(tbl2, file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.tsv.gz'))
write(date(), file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.done'))

if('databaseGroup' %in% names(opt)){
  library(RMariaDB)
  
  dbConn <- tryCatch({
    dbConnect(RMariaDB::MariaDB(), group = opt$databaseGroup)
  },
  error=function(cond) {
    write(c(paste(now(), '   Error - could not connect to the database.')), file = file.path(opt$outputDir, opt$buildSites_outputDir, 'log'), append = TRUE)
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
      write(c(paste(now(), 'Error -- could not upload sites data for ', x$sample[1], ' to the database.')), file = file.path(opt$outputDir, opt$buildSites_outputDir, 'log'), append = TRUE)
      q(save = 'no', status = 1, runLast = FALSE)
    } else {
      write(c(paste(now(), '   Uploaded sites data for ', x$sample[1], ' to the database.')), file = file.path(opt$outputDir, opt$buildSites_outputDir, 'log'), append = TRUE)
    }
  }))
  
  dbDisconnect(dbConn)
}

q(save = 'no', status = 0, runLast = FALSE) 
