library(dplyr)
library(lubridate)
library(parallel)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

if(! 'core_createFauxFragDoneFiles' %in% names(opt)) opt$core_createFauxFragDoneFiles <- FALSE
if(! 'core_createFauxSiteDoneFiles' %in% names(opt)) opt$core_createFauxSiteDoneFiles <- FALSE

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
  
  #browser()
  
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
      # browser()
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

o <- split(frags, frags$g)
i <- dplyr::ntile(1:length(o), opt$buildSites_CPUs)
f <- mapply(function(x, y){
       x$i <- y
       x
     }, o, i, SIMPLIFY = FALSE) %>% bind_rows()

cluster <- makeCluster(opt$buildSites_CPUs)
clusterExport(cluster, c('opt', 'frags'))

# sites <- bind_rows(parLapply(cluster, split(f, f$i), function(a){
# parallelization here is leading to an issue with represtantiveSeq() 
# not returning a result for some inputs.

sites <- bind_rows(lapply(split(f, f$i), function(a){  
           library(dplyr)
           source(file.path(opt$softwareDir, 'lib.R'))
  
           bind_rows(lapply(split(a, a$g), function(x){

             if(nrow(x) == 1){
               return(dplyr::mutate(x, UMIs = n_distinct(x$randomLinkerSeq), 
                                    sonicLengths = n_distinct(x$fragWidth), 
                                    maxLeaderSeqDist = 0) %>%
                      dplyr::select(trial, subject, sample, replicate, refGenome, posid, flags, UMIs, sonicLengths, reads, maxLeaderSeqDist, repLeaderSeq, vectorFastaFile))
           } else {
             s <- unlist(x$leaderSeqs)

            if(length(s) > opt$buildStdFragments_representativeSeqCalc_maxReads){
              set.seed(1)
              s <- sample(s, opt$buildStdFragments_representativeSeqCalc_maxReads)
            }
    
            r <- representativeSeq(s, tmpDir = file.path(opt$outputDir, opt$buildSites_outputDir, 'tmp'))
    
            return(dplyr::mutate(x, UMIs = n_distinct(x$randomLinkerSeq), 
                                sonicLengths = n_distinct(x$fragWidth), 
                                reads = sum(reads), 
                                repLeaderSeq = r[[2]], 
                                maxLeaderSeqDist = max(stringdist::stringdistmatrix(s))) %>%
                  dplyr::select(trial, subject, sample, replicate, refGenome, posid, flags, UMIs, sonicLengths, reads, maxLeaderSeqDist, repLeaderSeq, vectorFastaFile) %>%
                  dplyr::slice(1))
               }
           }))
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
                   t <- tibble(x1 = o$UMIs, x2 = o$sonicLengths, x3 = o$reads, x4 = o$repLeaderSeq)
                 } else if(nrow(o) > 1){
                   stop('Row error 1')  
                 } else {
                   t <- tibble(x1 = NA, x2 = NA, x3 = NA, x4 = NA)
                 }
           
                 names(t) <- c(paste0('rep', r, '-UMIs'), 
                               paste0('rep', r, '-sonicLengths'), 
                               paste0('rep', r, '-reads'), 
                              paste0('rep', r, '-repLeaderSeq'))
                t
              }))
         
          bind_cols(tibble(trial = x$trial[1], subject = x$subject[1], sample = x$sample[1], refGenome = x$refGenome[1], posid = x$posid[1], flags = x$flags[1], vectorFastaFile = x$vectorFastaFile[1]), o)
}))

write(c(paste(now(), '   Buiding sample level integration sites.')), file = file.path(opt$outputDir, opt$buildSites_outputDir, 'log'), append = TRUE)

tbl1$i <- dplyr::ntile(1:nrow(tbl1), opt$buildSites_CPUs)

tbl2 <- bind_rows(parLapply(cluster, split(tbl1, tbl1$i), function(p){
          library(dplyr)
          source(file.path(opt$softwareDir, 'lib.R'))
  
          bind_rows(lapply(1:nrow(p), function(n){
             x <- p[n,]

             f <- subset(frags, trial == x$trial & subject == x$subject & sample == x$sample & posid == x$posid)
  
             if(nrow(f) == 1){
               k <- tibble(UMIs = n_distinct(f$randomLinkerSeq), 
                           sonicLengths = n_distinct(f$fragWidth), 
                           reads = f$reads, 
                           maxLeaderSeqDist = f$maxLeaderSeqDist, 
                           repLeaderSeq = f$repLeaderSeq) 
    
             } else {
               s <- unlist(f$leaderSeqs)
               
               if(length(s) > opt$buildStdFragments_representativeSeqCalc_maxReads){
                 set.seed(1)
                 s <- sample(s, opt$buildStdFragments_representativeSeqCalc_maxReads)
               }
    
               r <- representativeSeq(s,tmpDir = file.path(opt$outputDir, opt$buildSites_outputDir, 'tmp'))
    
               k <- tibble(UMIs = n_distinct(f$randomLinkerSeq), 
                           sonicLengths = n_distinct(f$fragWidth), 
                           reads = sum(f$reads), 
                           maxLeaderSeqDist = max(stringdist::stringdistmatrix(s)), 
                           repLeaderSeq = r[[2]])
             }
  
             k$nRepsObs <- sum(! is.na(unlist(x[, which(grepl('\\-repLeaderSeq', names(x)))])))
             bind_cols(x[,1:5], k, x[,6:length(x)])
          }))
}))

if(any(is.infinite(tbl2$maxLeaderSeqDist))) tbl2[is.infinite(tbl2$maxLeaderSeqDist),]$maxLeaderSeqDist <- NA
tbl2$i <- NULL

tbl2 <- arrange(tbl2, desc(sonicLengths))

s <- unique(paste0(tbl2$trial, '~', tbl2$subject, '~', tbl2$sample))
if(any(! incomingSamples %in% s) & opt$core_createFauxSiteDoneFiles) core_createFauxSiteDoneFiles()

saveRDS(tbl2, file.path(opt$outputDir, opt$buildSites_outputDir, 'sites.rds'))
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
  
  invisible(lapply(split(tbl2, paste(tbl2$trial, tbl2$subject, tbl2$sample, tbl2$refGenome)), function(x){
    dbExecute(dbConn, paste0("delete from sites where trial='", x$trial[1], "' and subject='", x$subject[1],
                          "' and sample='", x$sample[1], "' and refGenome='", x$refGenome[1], "'"))
    
    f <- tmpFile()
    readr::write_tsv(dplyr::select(x, -trial, -subject, -sample, -refGenome), file.path(opt$outputDir, opt$buildSites_outputDir, 'tmp', f))
    system(paste0('xz ', file.path(opt$outputDir, opt$buildSites_outputDir, 'tmp', f)))
    
    fp <- file.path(opt$outputDir, opt$buildSites_outputDir, 'tmp', paste0(f, '.xz'))
    
    tab <- readBin(fp, "raw", n = as.integer(file.info(fp)["size"])+100)
    
    invisible(file.remove(list.files(file.path(opt$outputDir, opt$buildSites_outputDir, 'tmp'), pattern = f, full.names = TRUE)))
    
    r <- dbExecute(dbConn,
                   "insert into sites values (?, ?, ?, ?, ?, ?, ?, ?)",
                   params = list(x$trial[1], x$subject[1], x$sample[1], x$refGenome[1],
                                 x$vectorFastaFile[1], x$flags[1], list(serialize(tab, NULL)), as.character(lubridate::today())))
    
    if(r == 0){
      write(c(paste(now(), 'Error -- could not upload sites data for ', x$sample[1], ' to the database.')), file = file.path(opt$outputDir, opt$buildSites_outputDir, 'log'), append = TRUE)
      q(save = 'no', status = 1, runLast = FALSE)
    } else {
      write(c(paste(now(), '   Uploaded sites data for ', x$sample[1], ' to the database.')), file = file.path(opt$outputDir, opt$buildSites_outputDir, 'log'), append = TRUE)
    }
  }))
  
  dbDisconnect(dbConn)
}

q(save = 'no', status = 0, runLast = FALSE) 
