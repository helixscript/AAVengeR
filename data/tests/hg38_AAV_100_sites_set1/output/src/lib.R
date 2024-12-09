updateLog <- function(msg, logFile = NULL){
  if(is.null(logFile)) logFile <- opt$defaultLogFile
  msg <- paste0(base::format(Sys.time(), "%m.%d.%Y %l:%M%P"), "\t", msg)
  write(msg, file = logFile, append = TRUE)
}


startModule <- function(args){
  configFile <- args[1]
  if(! file.exists(configFile)) stop('Error - the configuration file does not exists.')
  opt <- yaml::read_yaml(configFile)
  
  if(length(args) == 2){
    invisible(lapply(strsplit(unlist(strsplit(args[2], ',')), ':'), function(x){
      val <- x[2]
      if(! suppressWarnings(is.na(as.numeric(val))) & grepl('\\.', val)) val <- as.numeric(val)
      if(! suppressWarnings(is.na(as.numeric(val))) & ! grepl('\\.', val)) val <- as.integer(val)
      opt[[x[1]]] <<- val
    })) 
  }
  
  opt$configFile <- configFile
  
  optionsSanityCheck(opt)
}


optionsSanityCheck <- function(opt){
  userRequired <- c('mode', 'softwareDir', 'outputDir', 'databaseConfigFile', 'databaseConfigGroup',
                    'sequencingRunID', 'demultiplex_anchorReadsFile', 'demultiplex_adriftReadsFile', 
                    'demultiplex_index1ReadsFile', 'demultiplex_sampleDataFile', 'modules')
  
  if(! all(userRequired %in% names(opt))){
    msg <- paste0('Error, these user required parameters are missing from the configuration file: ', paste0(userRequired[! userRequired %in% names(opt)], collapse = ', '))
    message(msg)
    q(save = 'no', status = 1, runLast = FALSE) 
  }
  
  if(! opt$mode %in% c('integrase', 'AAV', 'transposase', 'manual')){
    message('Error, the mode parameter must be set to integrase, AAV, transposase, or manual.')
    q(save = 'no', status = 1, runLast = FALSE) 
  }
  
  if(tolower(opt$sequencingRunID) == 'none'){
    message("Error - the sequencingRunID paramter is not allowed to be set to 'none'.")
    q(save = 'no', status = 1, runLast = FALSE) 
  }
  
  # Read in the default values and supplement user configs where appropriate.
  defaultConfig <- yaml::read_yaml(file.path(opt$softwareDir, 'defaults.yml'))
  
  allowedToBeAdded <- names(defaultConfig)[! names(defaultConfig) %in% userRequired]
  
  toAdd <- allowedToBeAdded[! allowedToBeAdded %in% names(opt)]
  
  if(length(toAdd) > 0){
    for(x in toAdd){
      opt[[x]] <- defaultConfig[[x]]
    }
  }
  
  # The core_CPU value can not exceed the number of system cores.
  cores <- parallel::detectCores()
  if(opt$core_CPUs > cores) opt$core_CPUs <- cores
  
  # Create list of all options containing 'CPUs' except core_CPUs.
  a <- names(opt)[grepl('_CPUs', names(opt))]
  a <- a[a != 'core_CPUs']
  
  # Do not allow any non-core module exceed the number of system cores.
  # Set their values to core_CPUs value if requested except for cases when the congfig
  # file is written by the core module since it balances the CPUs based on read counts.
  
  if(length(a) > 0){
    for(x in a){
      if(opt[[x]] > cores) opt[[x]] <- cores
      if(opt$core_applyCoreCPUsToAllModules & ! opt$calledFromCore) opt[[x]] <- opt$core_CPUs
    }
  }
  
  # Set mode specific values.
  if(grepl('integrase', opt$mode, ignore.case = TRUE)){
    opt$alignReads_genomeAlignment_anchorRead_maxStartPos <- 3
    opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned <- 5
    opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned <- 5
    opt$prepReads_HMMmatchEnd <- TRUE
    opt$prepReads_HMMmatchTerminalSeq <- 'CA'
    opt$prepReads_limitLeaderSeqsWithQuickAlignFilter <- FALSE
    opt$prepReads_forceAnchorReadStartSeq <- FALSE
  } else if(grepl('AAV', opt$mode, ignore.case = TRUE)){
    opt$demultiplex_quickAlignFilter <- TRUE
    opt$prepReads_limitLeaderSeqsWithQuickAlignFilter <- TRUE
    opt$alignReads_genomeAlignment_anchorRead_maxStartPos <- 300
    opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned <- 5
    opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned <- 10
    opt$prepReads_forceAnchorReadStartSeq <- TRUE # Use provided start sequence if a leader seq model is not returned.
    opt$prepReads_HMMmatchEnd <- FALSE # Triggers leaderSeq extension in alignReads.R.
  } else if(grepl('transposase', opt$mode, ignore.case = TRUE)){
    opt$alignReads_genomeAlignment_anchorRead_maxStartPos <- 3
    opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned <- 5
    opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned <- 5
    opt$prepReads_HMMmatchEnd <- TRUE
    opt$prepReads_HMMmatchTerminalSeq <- 'TA'
  } else if(grepl('manual', opt$mode, ignore.case = TRUE)){
    # Do nothing.
  } else {
    stop('Error -- mode set to an unknown value.')
  }
  
  opt
}


ppNum <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)


createOuputDir <- function(){
  if(! dir.exists(opt$outputDir)){
    invisible(sapply(readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt')), message))
    
    dir.create(file.path(opt$outputDir))
    
    if(! dir.exists(opt$outputDir)){
      message('Error: Can not create the output directory.')
      q(save = 'no', status = 1, runLast = FALSE)
    }
  }
}


previousSampleDatabaseCheck <- function(samples){
  if(tolower(opt$databaseConfigGroup) != 'none'){
      conn <- createDBconnection()
      
      r <- bind_rows(lapply(split(samples, 1:nrow(samples)), function(x){
             o <- dbGetQuery(conn, paste0("select * from samples where trial = '", x$trial, "' and subject = '", x$subject, "' and sample = '", x$sample, "' and replicate = '", x$replicate, "' and refGenome = '", x$refGenome, "'"))
             x <- select(x, trial, subject, sample, replicate, refGenome)
             x$inDatabase <- nrow(o) > 0
             x
           }))
    
      write(pander::pandoc.table.return(r, style = 'simple'), file.path(opt$outputDir, 'databaseSampleCheck.txt'))
    
      dbDisconnect(conn)
      
      if(any(r$inDatabase)){
        message('Error -- one or more sample replicates are already in the database.')
        message('See report located here: ', file.path(opt$outputDir, 'databaseSampleCheck.txt'))
        return(TRUE)
      } else {
        return(FALSE)
      }
  } else {
    return(FALSE)
  }
}


runSam2Psl <- function(samFile, outputFile, cl){
   o <- data.table(line = readLines(samFile))
   o$n <- ntile(1:nrow(o), CPUs)
   o <- split(o, o$n)
   gc()

   sam2psl <- function(x){
     system(paste0(file.path(opt$softwareDir, 'bin', 'sam2psl.py'), ' -i ',  x, ' -o ', paste0(x, '.out')))
     invisible(file.remove(x))
     x <- readLines(paste0(x, '.out'))
     invisible(file.remove(paste0(x, '.out')))
     x
   }

   f <- lapply(o, function(x){
     tmpFile <- paste0('sam2psl.', paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''))
     write(x$line, file.path(opt$outputDir, tmpFile))
     file.path(opt$outputDir, tmpFile)
   })

   write(unlist(parLapply(cl, f, sam2psl)), outputFile)
}


percentSysMemUsed <- function(){
  m <- as.integer(unlist(strsplit(system('free', intern = TRUE)[2], '\\s+'))[2:3])
  (m[2] / m[1]) * 100
}


lowMemoryException <- function(stepDesc = 'none', maxMem = 98){
  if(percentSysMemUsed() > maxMem){
    msg <- paste0('Error - the system has exceeded ', maxMem, '% memory usage -- stopping all processes. Step desc: ', stepDesc)
    updateLog(msg, logFile = file.path(opt$outputDir, 'Error.log'))
    q(save = 'no', status = 1, runLast = FALSE) 
  }
}


waitForMemory <- function(stepDesc = 'none', minMem = 80, maxWaitSecs = 1200, sleepSecs = 10){
  gc()
  t = 0
  
  repeat{
    lowMemoryException(stepDesc = stepDesc)
    
    if(t > maxWaitSecs){
      msg <- paste0('Error - the software has waited for ', maxWaitSecs, 
                    ' seconds for at least ', 100 - minMem, '% of the system memory to be free.',
                    ' Stopping all processes. Step desc: ', stepDesc)
      updateLog(msg, logFile = file.path(opt$outputDir, 'Error.log'))
      q(save = 'no', status = 1, runLast = FALSE) 
    }
    
    m <- percentSysMemUsed()
    if(m <= minMem) break
    
    updateLog(paste0('Waiting for memory to be freed, total memory use currently ', sprintf("%.2f%%", m), 
                     '%, waiting for no more than ', minMem, '% usage - waited for ', t, ' seconds.'))
    
    Sys.sleep(sleepSecs)
    
    t <- t + sleepSecs
  }
}



tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }

# AAVengeR is provided with default configuration files for different types of analyses
# where all the expected options are provided. Removing options from the configuration file
# may break the software in some instances. Some modules, such as the core module, injects 
# non-public options and this function ensures that they are always set even when not needed
# so modules will not throw errors. 


# Helper function that creates the expected done files expected by the core module when modules fail.
# The core module would become stuck in a perpetual loop if the done files did not appear after a failure.

core_createFauxFragDoneFiles <- function(){
  if(! dir.exists(file.path(opt$outputDir, 'buildFragments'))) dir.create(file.path(opt$outputDir,'buildFragments'))
  write(date(), file = file.path(opt$outputDir, 'buildFragments', 'fragments.done'))
  write(date(), file = file.path(opt$outputDir, 'buildFragments', 'xxx'))
}


core_createFauxSiteDoneFiles <- function(){
  if(! dir.exists(file.path(opt$outputDir, 'buildSites'))) dir.create(file.path(opt$outputDir, 'buildSites'))
  write(date(), file = file.path(opt$outputDir,  'buildSites', 'sites.done'))
  write(date(), file = file.path(opt$outputDir,  'buildSites', 'xxx'))
}


# Check the system path for the presence of required system level software.

checkSoftware <- function(){
  s <- c('blat', 'cutadapt', 'hmmbuild', 'hmmsearch', 'blastn', 'makeblastdb', 'python2', 'cd-hit-est')
  
  o <- bind_rows(lapply(s, function(x){
         r <- tryCatch({
                system(paste0('which ', x), intern = TRUE)
              },
              error=function(cond) {
                return(NA)
              },
              warning=function(cond) {
                return(NA)
              })
           tibble(program = x, path = r)
        }))   
  
        o$program <- paste0('                       ', o$program)
        
        write(paste0('Paths of expected software packages:'), file = file.path(opt$outputDir, 'log'), append = TRUE)
        write(pandoc.table.return(o, style = "simple", split.tables = Inf, plain.ascii = TRUE, justify = 'left'), file.path(opt$outputDir, 'log'), append = TRUE)
  
       if(any(is.na(o$path))){
         write(c(paste('Error - one or more of the expected softwares are not in your search path.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
         q(save = 'no', status = 1, runLast = FALSE) 
       } 
}


CD_HIT_clusters <- function(x, dir, params){
  f <- tmpFile()

  writeXStringSet(x, file.path(dir, paste0(f, '.fasta')))
  
  comm <- paste0("cd-hit-est -i ", file.path(dir, paste0(f, '.fasta')), 
                " -o ", file.path(dir, f), " -T ", opt$buildStdFragments_CPUs, 
                "  ", params, ' 1> /dev/null')
  
  m <- system(comm)
  
  if(m != 0) quitOnErorr(paste0('cd-hit-est failed in run dir: ', dir))
  
  r <- paste0(readLines(paste0(file.path(dir, f), '.clstr')), collapse = '')
  invisible(file.remove(list.files(dir, pattern = f, full.names = TRUE)))

  o <- unlist(strsplit(r, '>Cluster'))
  o[2:length(o)]
}


# Last Path Element -- return the last element of a file path delimited by slashes.
lpe <- function(x){
  o <- unlist(strsplit(x, '/'))
  o[length(o)]
}


# Convert a ShortRead object to a BioString DNA object and
# removing trailing lane information from reads ids.
shortRead2DNAstringSet <- function(x){
  r <- x@sread
  names(r) <- sub('\\s+.+$', '', as.character(x@id))
  r
}


# Helper function to wait for a file to appear on the file system.
waitForFile <- function(f, seconds = 1){
  repeat
  {
    Sys.sleep(seconds)
    if(file.exists(f)) break
  }
  return(TRUE)
}


demultiplex <- function(x){
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ShortRead))
 
  source(file.path(opt$softwareDir, 'lib.R'))
  
  chunk <- x$n[1]
  
  f <- file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', x$file[grepl('index1', x$file)])
  if(! file.exists(f)) waitForFile(f)
  index1Reads <- readFastq(f)
  invisible(file.remove(f))
  
  waitForMemory(stepDesc = 'Demultiplex', minMem = opt$system_minMemThreshold, maxWaitSecs = opt$system_minMemWaitTime, sleepSecs = opt$system_minMemSleepTime)
  
  f <- file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', x$file[grepl('anchor', x$file)])
  if(! file.exists(f)) waitForFile(f)
  anchorReads <- readFastq(f)
  invisible(file.remove(f))
  
  waitForMemory(stepDesc = 'Demultiplex', minMem = opt$system_minMemThreshold, maxWaitSecs = opt$system_minMemWaitTime, sleepSecs = opt$system_minMemSleepTime)
  
  f <- file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', x$file[grepl('adrift', x$file)])
  if(! file.exists(f)) waitForFile(f)
  adriftReads <- readFastq(f)
  invisible(file.remove(f))
  
  waitForMemory(stepDesc = 'Demultiplex', minMem = opt$system_minMemThreshold, maxWaitSecs = opt$system_minMemWaitTime, sleepSecs = opt$system_minMemSleepTime)
  
  anchorReads <- trimTailw(anchorReads, opt$demultiplex_qualtrim_events, opt$demultiplex_qualtrim_code, opt$demultiplex_qualtrim_halfWidth)
  anchorReads <- anchorReads[width(anchorReads) >= opt$demultiplex_qualtrim_minLength]
  if(length(anchorReads) == 0) return()
  
  adriftReads <- trimTailw(adriftReads, opt$demultiplex_qualtrim_events, opt$demultiplex_qualtrim_code, opt$demultiplex_qualtrim_halfWidth)
  adriftReads <- adriftReads[width(adriftReads) >= opt$demultiplex_qualtrim_minLength]
  if(length(adriftReads) == 0) return()
  
  index1Reads <- shortRead2DNAstringSet(index1Reads)
  anchorReads <- shortRead2DNAstringSet(anchorReads)
  adriftReads <- shortRead2DNAstringSet(adriftReads)
  
  reads <- syncReads(index1Reads, anchorReads, adriftReads)
  index1Reads <- reads[[1]];  anchorReads  <- reads[[2]]; adriftReads  <- reads[[3]]
  
  if(length(adriftReads) == 0) return()
  
  # Correct Golay encoded barcodes if requested.
  if(opt$demultiplex_correctGolayIndexReads){

    tmpDir <- file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', tmpFile())
    dir.create(tmpDir)
    
    index1Reads.org <- index1Reads
    index1Reads <- golayCorrection(index1Reads, tmpDirPath = tmpDir)
    percentChanged <- (sum(! as.character(index1Reads.org) == as.character(index1Reads)) / length(index1Reads))*100
    
    updateLog(paste0('Golay correction complete. ', sprintf("%.2f", percentChanged), '% of reads updated via Golay correction.'))
    
    rm(index1Reads.org)
    invisible(unlink(tmpDir, recursive = TRUE))
    
    reads <- syncReads(index1Reads, anchorReads, adriftReads)
    index1Reads <- reads[[1]];  anchorReads  <- reads[[2]]; adriftReads  <- reads[[3]]
  }
  
  
  # Loop through samples in sample data file to demultiplex and apply read specific filters.
  invisible(lapply(1:nrow(samples), function(r){
    r <- samples[r,]
      
    v0 <- rep(TRUE, length(anchorReads))
    if('anchorReadStartSeq' %in% names(r)){
      v0 <- vcountPattern(r$anchorReadStartSeq, subseq(anchorReads, 1, nchar(r$anchorReadStartSeq)), max.mismatch = opt$demultiplex_anchorReadStartSeq.maxMisMatch) == 1
    }
    
    # Create barcode demultiplexing vectors.
    v1 <- vcountPattern(r$index1Seq, index1Reads, max.mismatch = opt$demultiplex_index1ReadMaxMismatch) > 0
    
    log.report <- tibble(sample = r$uniqueSample, demultiplexedIndex1Reads = sum(v1))
    
    # Create break read linker barcode demultiplexing vector.
    v2 <- rep(TRUE, length(adriftReads))
    if(opt$demultiplex_useAdriftReadUniqueLinkers){
      testSeq <- substr(r$adriftReadLinkerSeq, r$adriftRead.linkerBarcode.start, r$adriftRead.linkerBarcode.end)
      v2 <- vcountPattern(testSeq, subseq(adriftReads, r$adriftRead.linkerBarcode.start, r$adriftRead.linkerBarcode.end), max.mismatch = opt$demultiplex_adriftReadLinkerBarcodeMaxMismatch) > 0
      log.report$demultiplexedLinkerReads <- sum(v2)
    } else {
      log.report$demultiplexedLinkerReads <- NA
    }
    
    # Test to see if any reads demultiplex to this row of the sample table and then subset reads to this sample.
    i <- Reduce(base::intersect, list(which(v0), which(v1), which(v2)))
    if(length(i) == 0){
      log.report$demultiplexedReads <- 0
    } else {
      reads <- syncReads(index1Reads[i], anchorReads[i], adriftReads[i])
      index1Reads <- reads[[1]];  anchorReads  <- reads[[2]]; adriftReads  <- reads[[3]]
      
      if(length(index1Reads) == 0){
        log.report$demultiplexedReads <- 0
      } else {
        writeFasta(subseq(adriftReads, r$adriftRead.linkerRandomID.start, r$adriftRead.linkerRandomID.end), 
                   file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0(r$uniqueSample, '.', chunk, '.randomAdriftReadIDs')), compress = opt$compressDataFiles)
        
        writeFasta(anchorReads, file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0(r$uniqueSample, '.', chunk, '.anchorReads')), compress = opt$compressDataFiles)
        writeFasta(adriftReads, file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0(r$uniqueSample, '.', chunk, '.adriftReads')), compress = opt$compressDataFiles)
        log.report$demultiplexedReads <- length(index1Reads)
      }
    }
    
    write.table(log.report, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'logs', paste0(r$uniqueSample, '.', chunk, '.logReport')))
  }))
}


blat <- function(y, ref, dir){
  f <- file.path(dir, tmpFile())
  write(paste0('>', y$id, '\n', y$seq), f)
  
  occFlag <- ''
  if(opt$alignReads_blatUseOocFile){
    occFlag <- paste0(' -ooc=', file.path(opt$outputDir, opt$alignReads_outputDir, paste0(opt$alignReads_genomeAlignment_blatTileSize, '.ooc'))) 
  }
  
  system(paste0('blat ', ref, ' ', f, ' ', paste0(f, '.psl'), 
                ' -tileSize=', opt$alignReads_genomeAlignment_blatTileSize, 
                ' -stepSize=', opt$alignReads_genomeAlignment_blatStepSize, 
                ' -repMatch=', opt$alignReads_genomeAlignment_blatRepMatch, 
                ' -out=psl -t=dna -q=dna -minScore=0 -minIdentity=0 -noHead -noTrimA',
                occFlag, ifelse(opt$alignReads_blat_fastMap, ' -fastMap', '')), ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  write('done', paste0(f, '.done'))
}


alignReadEndsToVector <- function(y){
  suppressPackageStartupMessages(library(Biostrings))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(lubridate))
  
  s <- DNAStringSet(y$testSeq)
  names(s) <- y$readID
  
  b <- blastReads(s, tmpDirPath = file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'), dbPath = file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd'))
  
  if(nrow(b) > 0){
    b$alignmentLength <- b$qend - b$qstart + 1
    b <- dplyr::filter(b, pident >= opt$prepReads_vectorAlignmentTest_minPercentID, alignmentLength >= floor(opt$prepReads_vectorAlignmentTestLength * (opt$prepReads_vectorAlignmentTest_minPercentCoverage/100)), gapopen <= 1)
  }
  
  if(nrow(b) > 0){
    b <- left_join(b, data.table(qname = names(s), testSeq = as.character(s)), by = 'qname')
    b$start  <- ifelse(b$sstart > b$send, b$send, b$sstart)
    b$end    <- ifelse(b$sstart > b$send, b$sstart, b$send)
    b$strand <- ifelse(b$sstart > b$send, '-', '+')
    b <- group_by(b, qname) %>% dplyr::top_n(1, wt = bitscore) %>% dplyr::arrange(desc(strand)) %>% dplyr::slice(1) %>% ungroup()
  }
  
  b
}


reconstructDBtable <- function(o, tmpDirPath){
  r <- tibble()
  
  if(nrow(o) > 0){
    r <- bind_rows(lapply(split(o, 1:nrow(o)), function(y){
           f <- tmpFile()
           writeBin(unserialize(y$data[[1]]), file.path(tmpDirPath, paste0(f, '.xz')))
           system(paste0('unxz ', file.path(tmpDirPath, paste0(f, '.xz'))))
           d <- readr::read_tsv(file.path(tmpDirPath, f))
           invisible(file.remove(file.path(tmpDirPath, f)))
           d
         }))
  }
  
  r
}


pullDBsubjectFrags <- function(dbConn, trial, subject, tmpDirPath){
  o <- dbGetQuery(dbConn, paste0("select * from fragments where trial = '", trial, "' and subject = '", subject, "'"))
  reconstructDBtable(o, tmpDirPath)
}


pullDBsubjectSites <- function(dbConn, trial, subject, tmpDirPath){
  o <- dbGetQuery(dbConn, paste0("select * from sites where trial = '", trial, "' and subject = '", subject, "'"))
  reconstructDBtable(o, tmpDirPath)
}


syncReads <- function(...){
  arguments <- list(...)
  n <- Reduce(base::intersect, lapply(arguments, names))
  lapply(arguments, function(x) x[match(n, names(x))])
}


createDBconnection <- function(){
  suppressPackageStartupMessages(library(RMariaDB))
  
  if(! file.exists('~/.my.cnf')) file.copy(opt$databaseConfigFile, '~/.my.cnf')
  if(! file.exists('~/.my.cnf')) quitOnErorr('Error - can not find ~/.my.cnf file.')
  
  tryCatch({
    dbConnect(RMariaDB::MariaDB(), group = opt$databaseConfigGroup)
  },
  error=function(cond) {
    quitOnErorr('Error - could not connect to the database.')
  })
}


uploadSitesToDB <- function(sites){
  suppressPackageStartupMessages(library(RMariaDB))
  
  updateLog('Writing sites to the database.')

  if(! dir.exists(file.path(opt$outputDir, 'tmp'))) dir.create(file.path(opt$outputDir, 'tmp'))
  
  conn <- createDBconnection()
  
  invisible(lapply(split(sites, paste(sites$trial, sites$subject, sites$sample, sites$refGenome)), function(x){
    updateLog(paste0('Uploading sites for: ', sites$trial[1], '/', sites$subject[1], '/', sites$sample[1], '/',sites$refGenome[1]))
    
    dbExecute(conn, paste0("delete from sites where trial='", x$trial[1], "' and subject='", x$subject[1],
                           "' and sample='", x$sample[1], "' and refGenome='", x$refGenome[1], "'"))
    
    f <- tmpFile()
    readr::write_tsv(x, file.path(opt$outputDir, 'tmp', f))
    system(paste0('xz --best ', file.path(opt$outputDir, 'tmp', f)))
    
    fp <- file.path(opt$outputDir, 'tmp', paste0(f, '.xz'))
    
    updateLog(paste0('Reading tar ball as a byte stream (', file.size(fp), ' bytes).'))
    
    tab <- readBin(fp, "raw", n = as.integer(file.info(fp)["size"])+100)
    
    invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), pattern = f, full.names = TRUE)))
    
    updateLog('Pushing byte object to database.')
    
    r <- dbExecute(conn,
                   "insert into sites values (?, ?, ?, ?, ?)",
                   params = list(x$trial[1], x$subject[1], x$sample[1], x$refGenome[1], list(serialize(tab, NULL))))
    if(r == 0){
      quitOnErorr(paset0('Error -- could not upload site data for ', x$trial[1], '~', x$subject[1], '~', x$sample[1], ' to the database.'))
    } else {
      updateLog(paste0('Uploaded sites data for ',  x$trial[1], '~', x$subject[1], '~', x$sample[1], ' to the database.'))
    }
  }))
  
  dbDisconnect(conn)
}



loadSamples <- function(){
  
  samples <- readr::read_tsv(opt$demultiplex_sampleDataFile, col_types = readr::cols(), comment = "#")
  
  if(nrow(samples) == 0){
    updateLog('Error - no lines of information was read from the sample configuration file.')
    q(save = 'no', status = 1, runLast = FALSE) 
  }
  
  if('refGenome' %in% names(samples)){
    if(! all(sapply(unique(file.path(opt$softwareDir, 'data', 'referenceGenomes', 'blat', paste0(samples$refGenome, '.2bit'))), file.exists))){
      updateLog("Error - one or more blat database files could not be found in AAVengeR's data/referenceGenomes directory.")
      q(save = 'no', status = 1, runLast = FALSE) 
    }
  } else {
    samples$refGenome <- NA
  }
  
  if('vectorFastaFile' %in% names(samples)){
    if(! all(sapply(unique(file.path(opt$softwareDir, 'data', 'vectors', samples$vectorFastaFile)), file.exists))){
      updateLog("Error - one or more vector FASTA files could not be found in AAVengeR's data/vectors directory.")
      q(save = 'no', status = 1, runLast = FALSE) 
    }
  } else {
    samples$vectorFastaFile <- NA
  }
  
  if('leaderSeqHMM' %in% names(samples)){
    samples$leaderSeqHMM <- file.path(opt$softwareDir, 'data', 'hmms', samples$leaderSeqHMM)
    
    if(! all(sapply(unique(samples$leaderSeqHMM), file.exists))){
      updateLog("Error - one or more leader sequence HMM files could not be found in AAVengeR's data/hmms directory.")
      q(save = 'no', status = 1, runLast = FALSE) 
    }
  }

  # Create missing columns and set to NA.
  if(! 'adriftRead.linkerBarcode.start' %in% names(samples)){
    samples$adriftRead.linkerBarcode.start <- NA
    samples$adriftRead.linkerBarcode.end   <- NA
  }
  
  if(! 'adriftRead.linkerRandomID.start' %in% names(samples)){
    samples$adriftRead.linkerRandomID.start <- NA
    samples$adriftRead.linkerRandomID.end   <- NA
  }
  
  samples <- bind_rows(lapply(split(samples, 1:nrow(samples)), function(x){
    Ns <- stringr::str_extract(x$adriftReadLinkerSeq, 'N+')
    
    # Detecting Ns in the adrift linker and having adriftRead.linkerBarcode.start 
    # set to NA triggers the determination of positions.
    
    if(! is.na(Ns) & is.na(x$adriftRead.linkerBarcode.start)){  
      o <- stringr::str_locate(x$adriftReadLinkerSeq, Ns)
      x$adriftRead.linkerBarcode.start  <- 1
      x$adriftRead.linkerBarcode.end    <- o[1,1] - 1
    }
    
    # Detecting Ns in the adrift linker and having adriftRead.linkerRandomID.start 
    # set to NA triggers the determination of positions.
    
    if(! is.na(Ns) & is.na(x$adriftRead.linkerRandomID.start)){ 
      o <- stringr::str_locate(x$adriftReadLinkerSeq, Ns)
      x$adriftRead.linkerRandomID.start  <- o[1,1]
      x$adriftRead.linkerRandomID.end    <- o[1,2]
    }
    x
  }))
  
  requiredColumns <- c("trial", "subject", "sample", "replicate",  
                       "adriftReadLinkerSeq", "index1Seq", "refGenome", "vectorFastaFile", "flags")
  
  if(! all(requiredColumns %in% names(samples))){
    missingCols <- paste0(requiredColumns[! requiredColumns %in% names(samples)], collapse = ', ')
    updateLog(paste0('Error - the following columns were missing from the sample data file: ', missingCols))
    q(save = 'no', status = 1, runLast = FALSE)
  }
  
  if(any(grepl('\\s|~|\\||\\.', paste0(samples$trial, samples$subject, samples$sample, samples$replicate)))){
    updateLog("Error -- spaces, tildas (~), pipes (|), and dots (.) are reserved characters and can not be used in the trial, subject, sample, or replicate sample configuration columns.")
    q(save = 'no', status = 1, runLast = FALSE)
  }
  
  if(opt$demultiplex_useAdriftReadUniqueLinkers){
    if(any(is.na(samples$adriftRead.linkerBarcode.start)) | any(is.na(samples$adriftRead.linkerBarcode.end))){
      updateLog('Error - adriftRead.linkerBarcode.start or adriftRead.linkerBarcode.end is set to NA when requesting demultiplex_useAdriftReadUniqueLinkers')
      q(save = 'no', status = 1, runLast = FALSE)
    }
  }
  
  if(any(is.na(samples$adriftRead.linkerRandomID.start)) | any(is.na(samples$adriftRead.linkerRandomID.end))){
      updateLog('Error - adriftRead.linkerRandomID.start or adriftRead.linkerRandomID.end is set to NA.')
      q(save = 'no', status = 1, runLast = FALSE)
  }
  
  samples$uniqueSample <- paste0(samples$trial, '~', samples$subject, '~', samples$sample, '~', samples$replicate)
  
  if(any(duplicated(samples$uniqueSample))){
    updateLog('Error -- There are one ore more duplications of trial~subject~sample~replicate ids in the sample data file.')
    q(save = 'no', status = 1, runLast = FALSE)
  }
  
  samples
}


golayCorrection <- function(x, tmpDirPath = NA){
  suppressPackageStartupMessages(library(ShortRead))
  suppressPackageStartupMessages(library(dplyr))
  
  f <- file.path(tmpDirPath, paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') )
  writeFasta(x, file = f)
  
  system(paste('python2', file.path(opt$softwareDir, 'bin', 'golayCorrection', 'processGolay.py'), f))
  invisible(file.remove(f))
  
  corrected <- readFasta(paste0(f, '.corrected'))
  invisible(file.remove(paste0(f, '.corrected')))
  
  a <- left_join(tibble(id = names(x), seq = as.character(x)),
                 tibble(id = as.character(corrected@id), seq2 = as.character(sread(corrected))), by = 'id')
  a <- rowwise(a) %>% mutate(editDist = as.integer(adist(seq, seq2))) %>% ungroup()
  
  i <- which(a$editDist <= 2)
  a[i,]$seq <- a[i,]$seq2
  
  if(! all(names(x) == a$id)){
    updateLog('Error - There was an ordering error during the Golay correction step.')
    q(save = 'no', status = 1, runLast = FALSE) 
  }
  
  r <- DNAStringSet(a$seq)
  names(r) <- a$id
  r
}


buildRearrangementModel <- function(b, seqsMinAlignmentLength){
  r <- vector()
  counter <- 0
  b <- b[(b$qend - b$qstart + 1) >= seqsMinAlignmentLength,]
  b$qstartBin <- cut(b$qstart, c(seq(0, max(b$qstart), 5), Inf), labels = FALSE)
  b <- arrange(b, qstartBin, desc(bitscore))

  while(nrow(b) != 0){
    counter <- counter + 1
    b <- arrange(b, qstart, evalue)
    x <- b[1,]
    r <- paste0(r, ';', x$qstart, '..', x$qend, '[', x$sstart2, x$strand, x$send2, ']')

    b <- b[b$qstart >= x$qend - 5 & b$qstartBin != x$qstartBin,] 
    b <- arrange(b, qstartBin, desc(bitscore))
    
    if(counter == 1000) break
  } 
  
  sub('^;', '', r)
}


blast2rearangements <- function(b, maxMissingTailNTs, minLocalAlignmentLength){
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(IRanges))
  suppressPackageStartupMessages(library(data.table))
  
  bind_rows(lapply(split(b, b$qname), function(b){
    
    # Sort BLAST results by query start position and evalue (low to high).
    b$strand <- ifelse(b$send < b$sstart, '-', '+')
    
    # Alignment to the negative strand will result in the subject end to come before the start.
    # Switch it back so that they are sequential.
    b$sstart2 <- ifelse(b$sstart > b$send, b$send, b$sstart)
    b$send2   <- ifelse(b$sstart > b$send, b$sstart, b$send)
    
    r <- tibble(qname = b$qname[1], rearrangement = buildRearrangementModel(b, minLocalAlignmentLength))
    
    k <- unlist(stringr::str_extract_all(r$rearrangement, '\\.\\.\\d+'))
    lastRangeEnd <- as.integer(sub('\\.\\.', '', stringr::str_extract(k[length(k)], '\\.\\.\\d+')))
    
    # Add missing internal segments.
    o <- unlist(strsplit(r$rearrangement, ';'))
    
    if(length(o) > 1){
      for(i in c(2:length(o))){
        n1 <- as.integer(sub('\\[', '', stringr::str_extract(o[i-1], '\\d+\\[')))
        n2 <- as.integer(stringr::str_extract(o[i], '^\\d+'))
        if(n2 - n1 >= 10) o[i-1] <- paste0(o[i-1], ';', n1+1, '..', n2-1, '[x]')
      }
      
      r$rearrangement <- paste0(o, collapse = ';')
    }
      
    # Add an [x] segment to the end of models if it falls short of the end of the analyzed sequence.
    if(b$qlen[1] - lastRangeEnd >= maxMissingTailNTs) r$rearrangement <- paste0(r$rearrangement, ';', (lastRangeEnd+1), '..', b$qlen[1], '[x]')

    r
  }))
}


captureHMMleaderSeq <- function(reads, hmm, tmpDirPath = NA){
  outputFile <- file.path(tmpDirPath, tmpFile())
  
  if(opt$prepReads_useDefaultHMMsettings){
    localOpt <- opt # Force local copy opt opt.
    
    f <- sub('\\.hmm$', '.settings', hmm)
    if(! file.exists(f)) stop(paste0('Error -- hmm settings file: ', f, ' does not exists'))
    
    # Update local copy.
    y <- yaml::read_yaml(f)  
    for(i in names(y)){
      localOpt[[i]] <- y[[i]]
    }
    
    # rename local copy to original which will now be local and updated.
    opt <- localOpt
  }
  
  endPos <- ifelse(width(reads) > opt$prepReads_HMMsearchReadEndPos, opt$prepReads_HMMsearchReadEndPos, width(reads))
  writeXStringSet(subseq(reads, opt$prepReads_HMMsearchReadStartPos, endPos), outputFile)
  
  comm <- paste0('hmmsearch --max --tblout ', outputFile, '.tbl --domtblout ', outputFile, '.domTbl ', hmm, ' ', outputFile, ' > ', outputFile, '.hmmsearch')
  system(comm)
    
  o <- readr::read_table(paste0(outputFile, '.domTbl'), col_names = FALSE, col_types = NULL, comment = "#")
  
  if(nrow(o) == 0) return(tibble())
  
  names(o) <- c('targetName', 'targetAcc', 'tlen', 'queryName', 'queryAcc', 'queryLength', 'fullEval', 
                'fullScore', 'fullBias', 'domNum', 'totalDoms', 'dom_c-Eval', 'dom_i-Eval', 'domScore', 
                'domBias', 'hmmStart', 'hmmEnd', 'targetStart', 'targetEnd', 'envStart', 'envEnd', 
                'meanPostProb',  'desc') 
    
  h <- readLines(hmm)
  hmmLength <- as.integer(unlist(strsplit(h[grepl('^LENG', h)], '\\s+'))[2])
  hmmName <- unlist(strsplit(h[grepl('^NAME', h)], '\\s+'))[2]
  
  # Subset HMM results such that alignments start at the start of reads, the end of the HMM
  # which contains the CA is includes and the alignment has a significant alignment scores.
  
  if(opt$prepReads_HMMmatchEnd){
    o <- subset(o, targetStart <= opt$prepReads_HMMmaxStartPos &
                  hmmEnd == hmmLength &
                  fullScore >= as.numeric(opt$prepReads_HMMminFullBitScore))
  } else {
    o <- subset(o, targetStart <= opt$prepReads_HMMmaxStartPos &
                  fullScore >= as.numeric(opt$prepReads_HMMminFullBitScore))
  }
  
  if(nrow(o) == 0) return(tibble())
  
  # Limit reads to those with matching HMM hits.
  reads2 <- reads[names(reads) %in% o$targetName]
  
  rm(reads)
  suppressMessages(gc())
  
  # Arrange reads to match o data frame.
  reads2 <- reads2[match(o$targetName, names(reads2))]
  
  if(! grepl('none', opt$prepReads_HMMmatchTerminalSeq, ignore.case = TRUE)){
    reads2 <- reads2[as.character(subseq(reads2, o$targetEnd-(nchar(opt$prepReads_HMMmatchTerminalSeq)-1), o$targetEnd)) == opt$prepReads_HMMmatchTerminalSeq]
  }
  
  if(length(reads2) == 0) return(tibble())
  
  tibble(id = names(reads2),
         LTRname = hmmName,
         LTRseq = as.character(subseq(reads2, 1, o$targetEnd)))
}


blastReads <- function(reads, wordSize = 5, evalue = 100, tmpDirPath = NA, dbPath = NA){
  suppressPackageStartupMessages(library(Biostrings))
  suppressPackageStartupMessages(library(dplyr))
  
  f <- tmpFile()
  writeXStringSet(reads,  file.path(tmpDirPath, paste0(f, '.fasta')))
  
  comm <- paste0('blastn -dust no -soft_masking false -word_size ', wordSize, ' -evalue ', evalue,' -outfmt 6 -query ',
                 file.path(tmpDirPath, paste0(f, '.fasta')), ' -db ',
                 dbPath, ' -num_threads 1 -out ', file.path(tmpDirPath, paste0(f, '.blast')))
  
  system(comm, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  waitForFile(file.path(tmpDirPath, paste0(f, '.blast')))
  
  if(file.info(file.path(tmpDirPath, paste0(f, '.blast')))$size > 0){
    b <- read.table(file.path(tmpDirPath, paste0(f, '.blast')), sep = '\t', header = FALSE)
    names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
    
    o <- reads[names(reads) %in% b$qname]
    d <- tibble(qname = names(o), qlen = width(o))
    b <- left_join(b, d, by = 'qname')
    return(b)
  } else {
    return(tibble())
  }
}



nearestGene <- function(posids, genes, exons, CPUs = 20){
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(parallel))
  
  o <- base::strsplit(posids, '[\\+\\-]')
  d <- tibble(chromosome = unlist(lapply(o, '[', 1)),
              position = unlist(lapply(o, '[', 2)),
              strand = stringr::str_extract(posids, '[\\+\\-]'))
  d$n <- ntile(1:nrow(d), CPUs)
  d$position <- as.integer(d$position)
  
  cluster <- makeCluster(CPUs)
  clusterSetRNGStream(cluster, 1)
  clusterExport(cluster, c('genes', 'exons'), envir = environment())
  
  r <- bind_rows(parLapply(cluster, split(d, d$n), function(x){
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(GenomicRanges))
    
    x$exon <- NA
    x$nearestGene <- NA
    x$nearestGeneStrand <- NA
    x$nearestGeneDist <- NA
    x$inGene <- FALSE
    x$inExon <- FALSE
    x$beforeNearestGene <- NA
    
    r <- GenomicRanges::makeGRangesFromDataFrame(x, 
                                                 seqnames.field = 'chromosome',
                                                 start.field = 'position',
                                                 end.field = 'position')
    
    o <- suppressWarnings(GenomicRanges::findOverlaps(r, exons, select='all', ignore.strand=TRUE, type='any'))
    
    if(length(queryHits(o)) > 0){
      for(i in unique(queryHits(o))){
        x[i,]$exon <- paste0(unique(exons[subjectHits(o[queryHits(o) == i]),]$name2), collapse = ', ')
        x[i,]$nearestGeneStrand <- paste0(unique(as.character(strand(exons[subjectHits(o[queryHits(o) == i])]))), collapse = ',')
      }
      
      if(any(! is.na(x$exon))){
        x[! is.na(x$exon),]$nearestGene <- x[! is.na(x$exon),]$exon
        x[! is.na(x$exon),]$nearestGeneDist <- 0
      }
    }
    
    x1 <- tibble()
    x2 <- tibble()
    
    if(any(! is.na(x$exon))) x1 <- x[! is.na(x$exon),]  # Exon hits
    if(any(is.na(x$exon)))   x2 <- x[is.na(x$exon),]    # Non-exon hits
    
    if(nrow(x2) > 0){
      r <- GenomicRanges::makeGRangesFromDataFrame(x2, 
                                                   seqnames.field = 'chromosome',
                                                   start.field = 'position',
                                                   end.field = 'position')
    
      o <- suppressWarnings(GenomicRanges::distanceToNearest(r, genes, select='all', ignore.strand=TRUE))
      
      if(length(queryHits(o)) > 0){
        for(i in unique(queryHits(o))){
          x2[i,]$nearestGene <- paste0(unique(genes[subjectHits(o[queryHits(o) == i]),]$name2), collapse = ', ')
          x2[i,]$nearestGeneStrand <- paste0(unique(as.character(strand(genes[subjectHits(o[queryHits(o) == i])]))), collapse = ',')
          x2[i,]$nearestGeneDist <- unique(mcols(o[queryHits(o) == i])$distance) # Always single?
          x2[i,]$beforeNearestGene <- ifelse(x2[i,]$position < min(start(genes[subjectHits(o[queryHits(o) == i])])), TRUE, FALSE)
        }
      }
      
      x2 <- x2[! is.na(x2$nearestGene),]
      
      if(any(x2$nearestGeneDist == 0)) x2[which(x2$nearestGeneDist == 0),]$beforeNearestGene <- NA
    }
    
    if(nrow(x1) > 0  & nrow(x2) > 0)  x <- bind_rows(x1, x2)
    if(nrow(x1) == 0 & nrow(x2) > 0)  x <- x2
    if(nrow(x1) > 0  & nrow(x2) == 0) x <- x1
    
    x
  }))
  
  stopCluster(cluster)
  
  if(any(r$nearestGeneDist == 0)) r[which(r$nearestGeneDist == 0),]$inGene <- TRUE
  if(any(! is.na(r$exon))) r[which(! is.na(r$exon)),]$inExon <- TRUE
  if('n' %in% names(r)) r$n <- NULL
  if('exon' %in% names(r)) r$exon <- NULL
  r
}


parseBLAToutput <- function(f){
  if(! file.exists(f) | file.info(f)$size == 0) return(tibble::tibble())
  b <- readr::read_delim(f, delim = '\t', col_names = FALSE, col_types = 'iiiiiiiicciiiciiiiccc')
  
  x <- read.table(textConnection(system(paste(file.path(opt$softwareDir, 'bin', 'pslScore.pl') ,  f), intern = TRUE)), sep = '\t')
  names(x) <- c('tName', 'tStart', 'tEnd', 'hit', 'pslScore', 'percentIdentity')
  
  names(b) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand',
                'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')
  
  b$queryPercentID <- x$percentIdentity
  b$pslScore <- x$pslScore
  b$tStart <- x$tStart
  b$tEnd <- x$tEnd - 1 # Correct for zero-based half open coordinate system.
  
  dplyr::select(b, qName, matches, strand, qSize, qStart, qEnd, tName, tNumInsert, qNumInsert, tBaseInsert, qBaseInsert, tStart, tEnd, queryPercentID, pslScore)
}


blast2rearangements_worker <-  function(b){
  library(IRanges)
  library(dplyr)
  library(data.table)
  
  data.table::rbindlist(lapply(split(b, b$qname), function(b2){
    
    # Alignment to the negative strand will result in the subject end to come before the start.
    # Switch it back so that they are sequential.
    b2 <- arrange(b2, qstart, evalue)
    b2$sstart2 <- ifelse(b2$sstart > b2$send, b2$send, b2$sstart)
    b2$send2   <- ifelse(b2$sstart > b2$send, b2$sstart, b2$send)
    
    
    ir <- IRanges(start = b2$qstart, end = b2$qend)
    if(length(ir) == 0) return(data.table())
    
    r <- IRanges::reduce(ir, min.gapwidth = opt$prepReads_mapLeaderSeqsMaxGapBetweenAlignments)
    r <- r[start(r) <= opt$prepReads_buildReadMaps_minMapStartPostion]
    if(length(r) == 0) return(data.table())
    
    o <- data.table(qname = b2$qname[1], end = IRanges::end(r), width = width(r))
    o <- o[o$width == max(o$width),]
    o[1,]
  }))
}


intSiteCallerConversions <- function(){
  r <- list()
  r[['vector_CART19.fa']][['vector']]   = 'Bushman_CART19_vector.fasta'
  r[['vector_CART19.fa']][['hmm']]      = 'Bushman_CART19.hmm'
  r[['vector_CART19.fa']][['startSeq']] = 'GAAAATC'
  
  r
}
