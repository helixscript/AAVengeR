

tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }

# AAVegeneR is provided with default configuration files for different types of analyses
# where all the expected options are provided. Removing options from the configuration file
# may break the software in some instances. Some modules, such as the core module, injects 
# non-public options and this function ensures that they are always set even when not needed
# so modules will not throw errors. 

setMissingOptions <- function(){
  if(! 'core_createFauxFragDoneFiles' %in% names(opt)) opt$core_createFauxFragDoneFiles <<- FALSE
  if(! 'core_createFauxSiteDoneFiles' %in% names(opt)) opt$core_createFauxSiteDoneFiles <<- FALSE
}


# Helper function that creates the expected done files expected by the core module when modules fail.
# The core module would become stuck in perpetual loops if the done files did not appear after a failure.

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
  s <- c('blastn', 'blat', 'cutadapt', 'hmmbuild', 'hmmsearch', 'mafft', 'makeblastdb', 'muscle', 'python2')
  
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
        
        write(paste0(now(),         '    Paths of expected software packages:'), file = file.path(opt$outputDir, 'log'), append = TRUE)
        readr::write_tsv(o, file = file.path(opt$outputDir, 'log'), append = TRUE, col_names = FALSE)
  
       if(any(is.na(o$path))){
         write(c(paste(lubridate::now(), 'Error - one or more of the expected softwares are not in your search path.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
         q(save = 'no', status = 1, runLast = FALSE) 
       } 
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


# Pull non-standardized fragments from a database based on trial and subject ids.

pullDBsubjectFrags <- function(dbConn, trial, subject, tmpDirPath){
  
  o <- dbGetQuery(dbConn, paste0("select * from fragments where trial = '", trial, "' and subject = '", subject, "'"))
  r <- tibble()
  
  if(nrow(o) > 0){
    r <- bind_rows(lapply(split(o, 1:nrow(o)), function(y){
      f <- tmpFile()
      writeBin(unserialize(o$data[[1]]), file.path(tmpDirPath, paste0(f, '.xz')))
      system(paste0('unxz ', file.path(tmpDirPath, paste0(f, '.xz'))))
      d <- readr::read_tsv(file.path(tmpDirPath, f))
      invisible(file.remove(file.path(tmpDirPath, f)))
      d
    }))
  }
  
  r
}


# Quality trim ShortRead reads based on parameters in the configuration file 
# and remove reads that fall below a specific length post trimming. 

qualTrimReads <- function(f, chunkSize, label, ouputDir){
  
  # Create a pointer like object to the read file.
  strm <- FastqStreamer(f, n = as.integer(chunkSize))
  
  if(! 'demultiplex_qualtrim_events' %in% names(opt)) opt$demultiplex_qualtrim_events <- 2
  if(! 'demultiplex_qualtrim_halfWidth' %in% names(opt)) opt$demultiplex_qualtrim_halfWidth <- 5
  
  n <- 1
  repeat {
    fq <- yield(strm)
    if(length(fq) == 0) break
    
    fq <- trimTailw(fq, opt$demultiplex_qualtrim_events, opt$demultiplex_qualtrim_code, opt$demultiplex_qualtrim_halfWidth)
    fq <- fq[width(fq) >= opt$demultiplex_qualtrim_minLength]
    
    if(length(fq) > 0) writeFastq(fq, file = file.path(ouputDir, paste0(label, '.', n)), compress = FALSE)
    
    n <- n + 1
  }
}


# Given a set of Biostring objects, subject each object in the set such that 
# they share a common set of read ids.

syncReads <- function(...){
  arguments <- list(...)
  n <- Reduce(base::intersect, lapply(arguments, names))
  lapply(arguments, function(x) x[match(n, names(x))])
}


# Helper function to deconstruct uniqueSampleIDs into separate data frame columns.

unpackUniqueSampleID <- function(d){
  d$trial     <- unlist(lapply(strsplit(d$uniqueSample, '~'), '[[', 1))
  d$subject   <- unlist(lapply(strsplit(d$uniqueSample, '~'), '[[', 2))
  d$sample    <- unlist(lapply(strsplit(d$uniqueSample, '~'), '[[', 3))
  d$replicate <- unlist(lapply(strsplit(d$uniqueSample, '~'), '[[', 4))
  d
}


# Function to remove minor sequence variation from short sequences by identifying
# a abundant sequences and conforming those with minor differences to those sequences.

conformMinorSeqDiffs <- function(x, editDist = 1, abundSeqMinCount = 10, nThreads = 10){
  tab <- table(unname(x))
  if(length(tab[tab >= abundSeqMinCount]) == 0) return(x)
  
  abundant_codes <- names(tab[tab >= abundSeqMinCount])
  nonAbundant_codes <- names(tab[tab < abundSeqMinCount])
  
  # All codes are the abundant code.
  if(length(nonAbundant_codes) == 0) return(x)
  
  conversion_table <- tibble(a = nonAbundant_codes,
                             b = unlist(lapply(nonAbundant_codes, function(x){
                               z <- stringdist::stringdist(x, abundant_codes, nthread = nThreads)
                               i <- which(z <= editDist)
                               if(length(i) == 1) x <- abundant_codes[i]
                               if(length(i) > 1) x <- paste0(rep('N', nchar(x)), collapse = '')
                               x  
                             }))) %>% filter(a != b)
  
  unlist(lapply(x, function(y){
    if(y %in% conversion_table$a) y <- conversion_table[match(y, conversion_table$a),]$b
    y
  }))
}



loadSamples <- function(){
  samples <- readr::read_tsv(opt$demultiplex_sampleDataFile, col_types = readr::cols())

  # browser()
  
  if(nrow(samples) == 0){
    write(c(paste(lubridate::now(), 'Error - no lines of information was read from the sample configuration file.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE) 
  }
  
  if('refGenome' %in% names(samples)){
    if(! all(sapply(unique(file.path(opt$softwareDir, 'data', 'blatDBs', paste0(samples$refGenome, '.2bit'))), file.exists))){
      write(c(paste(now(), "Error - one or more blat database files could not be found in AAVengeR's data/blatDBs directory")), file = file.path(opt$outputDir, 'log'), append = TRUE)
      q(save = 'no', status = 1, runLast = FALSE) 
    }
  } else {
    samples$refGenome <- NA
  }
  
  if('vectorFastaFile' %in% names(samples)){
    if(! all(sapply(unique(file.path(opt$softwareDir, 'data', 'vectors', samples$vectorFastaFile)), file.exists))){
      write(c(paste(now(), "Error - one or more vector FASTA files could not be found in AAVengeR's data/vectors directory")), file = file.path(opt$outputDir, 'log'), append = TRUE)
      q(save = 'no', status = 1, runLast = FALSE) 
    }
  } else {
    samples$vectorFastaFile <- NA
  }
  
  if('leaderSeqHMM' %in% names(samples)){
    samples$leaderSeqHMM <- file.path(opt$softwareDir, 'data', 'hmms', samples$leaderSeqHMM)
    
    if(! all(sapply(unique(samples$leaderSeqHMM), file.exists))){
      write(c(paste(now(), "Error - one or more leader sequence HMM files could not be found in AAVengeR's data/hmms directory")), file = file.path(opt$outputDir, 'log'), append = TRUE)
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
    write(c(paste(lubridate::now(), 'Error - the following columns were missing from the sample data file: '), missingCols), file = file.path(opt$outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE)
  }
  
  if(any(grepl('\\s|~|\\||\\.', paste0(samples$trial, samples$subject, samples$sample, samples$replicate)))){
    write('Error -- spaces, tildas (~), pipes (|), and dots (.) are reserved characters and can not be used in the trial, subject, sample, or replicate sample configuration columns.', 
          file = file.path(opt$outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE)
  }
  
  if(opt$demultiplex_useAdriftReadUniqueLinkers){
    if(any(is.na(samples$adriftRead.linkerBarcode.start)) | any(is.na(samples$adriftRead.linkerBarcode.end))){
      write(c(paste(lubridate::now(), 'Error - adriftRead.linkerBarcode.start or adriftRead.linkerBarcode.end is set to NA when requesting demultiplex_useAdriftReadUniqueLinkers')), file = file.path(opt$outputDir, 'log'), append = TRUE)
      q(save = 'no', status = 1, runLast = FALSE)
    }
  }
  
  if(any(is.na(samples$adriftRead.linkerRandomID.start)) | any(is.na(samples$adriftRead.linkerRandomID.end))){
      write(c(paste(lubridate::now(), 'Error - adriftRead.linkerRandomID.start or adriftRead.linkerRandomID.end is set to NA.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
      q(save = 'no', status = 1, runLast = FALSE)
  }
  
  samples$uniqueSample <- paste0(samples$trial, '~', samples$subject, '~', samples$sample, '~', samples$replicate)
  
  if(any(duplicated(samples$uniqueSample))){
    write('Error -- There are one ore more duplications of trial~subject~sample~replicate ids in the sample data file.', 
          file = file.path(opt$outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE)
  }
  
  samples
}


parseBLAToutput <- function(f){
  if(! file.exists(f) | file.info(f)$size == 0) return(tibble::tibble())
  b <- readr::read_delim(f, delim = '\t', col_names = FALSE, col_types = readr::cols())
  
  names(b) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand',
                'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')
  
  b$tEnd <- b$tEnd - 1
  
  b$queryPercentID       <- (b$matches/b$qSize)*100
  b$tAlignmentWidth      <- (b$tEnd - b$tStart) + 1
  b$queryWidth           <- (b$qEnd - b$qStart) + 1
  b$alignmentPercentID   <- (b$matches/b$tAlignmentWidth)*100
  b$percentQueryCoverage <- (b$queryWidth/b$qSize)*100
  b$qStarts              <- as.character(b$qStarts)
  b$tStarts              <- as.character(b$tStarts)
  b
}


parseCutadaptLog <- function(f){
  l <- readLines(f)
  invisible(file.remove(f))

  if(length(l) > 0){
    tbl <- read.table(text = l, sep = '\t', fill = TRUE)

    if(length(tbl) >= 7){
      tbl <- tbl[,1:7]
      names(tbl) <- c('read.id', 'numErr', 'adapterStart', 'adapterEnd', 'seqBeforeAdapter',
                      'adapterSeq', 'seqAfterAdapter')
    }
    write.table(tbl, sep = '\t', file = f, col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}


representativeSeq <- function(s, percentReads = 95, tmpDirPath = NA){
  if(length(s) == 1 | dplyr::n_distinct(s) == 1) return(list(0, s[1]))
  
  k <- data.frame(table(s))
  s <- unique(s)
  
  if(length(s) > opt$buildStdFragments_representativeSeqCalc_maxReads){
    set.seed(1)
    message('Sampling ', opt$buildStdFragments_representativeSeqCalc_maxReads, ' reads in representativeSeq()')
    s <- sample(s, opt$buildStdFragments_representativeSeqCalc_maxReads)
  }
  
  # Align sequences in order to handle potential indels.
  f <- tmpFile()
  inputFile <- file.path(tmpDirPath, paste0(f, '.fasta'))
  fileConn <- file(inputFile)
  write(paste0('>', paste0('s', 1:length(s)), '\n', s), file = fileConn)
  close(fileConn)
  
  # Align sequences.
  outputFile <- file.path(tmpDirPath, paste0(f, '.representativeSeq.muscle'))
  system(paste('muscle -quiet -maxiters 1 -diags -in ', inputFile, ' -out ', outputFile))
  if(! file.exists(outputFile)) waitForFile(outputFile)

  # Alignments are read in with dashes.
  s <- as.character(Biostrings::readDNAStringSet(outputFile))
  
  # Close the connection created by readDNAStringSet which are sometimes left often.
  g <- grepl(f, showConnections(all = TRUE)[,1])
  if(any(g)) close.connection(getConnection(which(grepl(f, g)) - 1))

  invisible(file.remove(c(inputFile, outputFile)))

  # Create an all vs. all edit distance matrix.
  m <- as.matrix(stringdist::stringdistmatrix(s, nthread = opt$buildStdFragments_representativeSeqCalc_CPUs))

  # Select the sequence with the least amount of difference to other sequences.
  d <- apply(m, 1, sum)
  lowestEditDistSeqs1 <- unique(s[which(d == min(d))])
  lowestEditDistSeqs2 <- gsub('\\-', '', lowestEditDistSeqs1)
  
  if(length(lowestEditDistSeqs1) > 1){
    lowestEditDistSeq <- dplyr::filter(k, s %in% lowestEditDistSeqs2) %>% slice_max(Freq, n = 1, with_ties = FALSE) %>% dplyr::pull(s) %>% as.character()
  } else {
    lowestEditDistSeq <- lowestEditDistSeqs2
  }
  
  # list(max(stringdist::stringdist(lowestEditDistSeqs1, s))/nchar(lowestEditDistSeq), lowestEditDistSeq)
  list(NA, lowestEditDistSeq)
}


standardizationSplitVector <- function(d, v){
  if(v == 'replicate'){
    return(paste(d$trial, d$subject, d$sample, d$replicate))
  } else if (v == 'sample'){
    return(paste(d$trial, d$subject, d$sample))
  } else if( v == 'subject'){
    return(paste(d$trial, d$subject))
  } else {
    write(c(paste(lubridate::now(), 'Error - standardization vector can not be created.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE) 
  }
}


standardizedFragments <- function(frags, opt, cluster){
  g <- GenomicRanges::makeGRangesFromDataFrame(frags, start.field = 'fragStart', end.field = 'fragEnd', keep.extra.columns = TRUE)
  g$s <- standardizationSplitVector(g, opt$buildStdFragments_standardizeSitesBy)
  g$s <- paste(g$s, g$repLeaderSeqGroup)
  
  if(opt$buildStdFragments_standardizeSites){
    g <- unlist(GenomicRanges::GRangesList(parallel::parLapply(cluster, split(g, g$s), function(x){
           library(dplyr)
           library(GenomicRanges)
           source(file.path(opt$softwareDir, 'lib.R'))
      
           x$intSiteRefined <- FALSE
         
           out <- tryCatch({
                             o <- standardize_sites(x, counts.col = 'reads', sata.gap = 5)
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
  } else {
    g$intSiteRefined <- FALSE
  }

  g$s <- standardizationSplitVector(g, opt$buildStdFragments_standardizeBreaksBy)
  g$s <- paste(g$s, g$repLeaderSeqGroup)
  
  if(opt$buildStdFragments_standardizeBreaks){
    g <- unlist(GenomicRanges::GRangesList(parallel::parLapply(cluster, split(g, g$s), function(x){
           library(dplyr) 
           library(GenomicRanges)
           source(file.path(opt$softwareDir, 'lib.R'))

           x$breakPointRefined <- FALSE
           out <- tryCatch({
                              o <- refine_breakpoints(x, counts.col = 'reads')
                              o$breakPointRefined <- TRUE
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
  } else {
    g$breakPointRefined <- FALSE
  }

  g$s <- NULL
  data.frame(g)
}


golayCorrection <- function(x, tmpDirPath = NA){
  library(ShortRead)
  library(dplyr)
  
  f <- file.path(tmpDirPath, paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') )
  writeFasta(x, file = f)
  
  system(paste('python2', file.path(opt$softwareDir, 'bin', 'golayCorrection', 'processGolay.py'), f))
  file.remove(f)
  
  corrected <- readFasta(paste0(f, '.corrected'))
  file.remove(paste0(f, '.corrected'))
  
  a <- left_join(tibble(id = names(x), seq = as.character(x)),
                 tibble(id = as.character(corrected@id), seq2 = as.character(sread(corrected))), by = 'id')
  a <- rowwise(a) %>% mutate(editDist = as.integer(adist(seq, seq2))) %>% ungroup()
  
  i <- which(a$editDist <= 2)
  a[i,]$seq <- a[i,]$seq2
  
  if(! all(names(x) == a$id)){
    write(c(paste(lubridate::now(), 'Error - There was an ordering error during the Golay correction step.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE) 
  }
  
  r <- DNAStringSet(a$seq)
  names(r) <- a$id
  r
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
    

buildRearrangementModel <- function(b, seqsMinAlignmentLength = 15){
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

blast2rearangements <- function(b, maxMissingTailNTs = 5){
  library(dplyr)
  library(IRanges)
  library(data.table)
  
  bind_rows(lapply(split(b, b$qname), function(b){
    
    # Sort BLAST results by query start position and evalue (low to high).
    b$strand <- ifelse(b$send < b$sstart, '-', '+')
    
    # Alignment to the negative strand will result in the subject end to come before the start.
    # Switch it back so that they are sequential.
    b$sstart2 <- ifelse(b$sstart > b$send, b$send, b$sstart)
    b$send2   <- ifelse(b$sstart > b$send, b$sstart, b$send)
    
    r <- tibble(qname = b$qname[1], rearrangement = buildRearrangementModel(b))
    
    # Add missing internal segments.
    o <- unlist(strsplit(r$rearrangement, ';'))
    if(length(o) > 1){
      for(i in c(2:length(o))){
        n1 <- as.integer(sub('\\[', '', stringr::str_extract(o[i-1], '\\d+\\[')))
        n2 <- as.integer(stringr::str_extract(o[i], '^\\d+'))
        if(n2 - n1 >= 10) o[i-1] <- paste0(o[i-1], ';', n1+1, '..', n2-1, '[x]')
      }
      r$rearrangement <- paste0(o, collapse = ';')
      
      o <- unlist(stringr::str_extract_all(r$rearrangement, '\\.\\.\\d+'))
      lastRangeEnd <- as.integer(sub('\\.\\.', '', stringr::str_extract(o[length(o)], '\\.\\.\\d+')))
      browser()
      if(b$qlen[1] - lastRangeEnd >= maxMissingTailNTs){
        
        r$rearrangement <- paste0(r$rearrangement, ';', (lastRangeEnd+1), '..', b$qlen[1], '[x]')
      }
    }
    
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
    
  reads <- reads[width(reads) > opt$prepReads_HMMsearchReadEndPos]
  if(length(reads) == 0) return(tibble())
  
  writeXStringSet(subseq(reads, opt$prepReads_HMMsearchReadStartPos, opt$prepReads_HMMsearchReadEndPos), outputFile)

  comm <- paste0('hmmsearch --max --tblout ', outputFile, '.tbl --domtblout ', outputFile, '.domTbl ', hmm, ' ', outputFile, ' > ', outputFile, '.hmmsearch')
  system(comm)

  r <- readLines(paste0(outputFile, '.domTbl'))
  r <- r[!grepl('^\\s*#', r)]
  r <- strsplit(r, '\\s+')
  
  o <- bind_rows(lapply(r, function(x) data.frame(t(x))))
  
  if(nrow(o) == 0) return(tibble())
  
  names(o) <- c('targetName', 'targetAcc', 'tlen', 'queryName', 'queryAcc', 'queryLength', 'fullEval', 
                'fullScore', 'fullBias', 'domNum', 'totalDoms', 'dom_c-Eval', 'dom_i-Eval', 'domScore', 
                'domBias', 'hmmStart', 'hmmEnd', 'targetStart', 'targetEnd', 'envStart', 'envEnd', 
                'meanPostProb',  'desc') 
  write.table(o, sep = '\t', file = paste0(outputFile, '.domTbl2'), col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  o <- readr::read_delim(paste0(outputFile, '.domTbl2'), '\t', col_types = readr::cols())
  o$fullScore <- as.numeric(o$fullScore)
  o$fullEval  <- as.numeric(o$fullEval)
  
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
  gc()
  
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
  library(Biostrings)
  library(dplyr)

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



blatReads <- function(reads, minIdentity=90, stepSize = 11, tileSize = 11, tmpDirPath = NA){
  library(Biostrings)
  library(dplyr)
  
  f <- tmpFile()
  writeXStringSet(reads,  file.path(tmpDirPath, paste0(f, '.fasta')))
  
  system(paste0(file.path(opt$softwareDir, 'bin', 'blat'), ' ', file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd'), ' ', 
                file.path(tmpDirPath, paste0(f, '.fasta')), ' ',
                file.path(tmpDirPath, paste0(f, '.psl')), ' -minIdentity=', minIdentity, 
                ' -stepSize=', stepSize, ' -tileSize=', tileSize, ' -out=psl -noHead -minScore=0'),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  waitForFile(file.path(tmpDirPath, paste0(f, '.psl')))
  
  if(file.info(file.path(tmpDirPath, paste0(f, '.psl')))$size > 0){
    return(parseBLAToutput(file.path(tmpDirPath, paste0(f, '.psl'))))
  } else {
    return(tibble())
  }
}


nearestGene <- function(posids, genes, exons, CPUs = 20){
  library(dplyr)
  library(parallel)
  
  o <- base::strsplit(posids, '[\\+\\-]')
  d <- tibble(chromosome = unlist(lapply(o, '[', 1)),
              position = unlist(lapply(o, '[', 2)),
              strand = stringr::str_extract(posids, '[\\+\\-]'))
  d$n <- ntile(1:nrow(d), CPUs)
  d$position <- as.integer(d$position)
  
  cluster <- makeCluster(CPUs)
  clusterExport(cluster, c('genes', 'exons'), envir = environment())
  
  r <- bind_rows(parLapply(cluster, split(d, d$n), function(x){
  #r <- bind_rows(lapply(split(d, d$n), function(x){
    library(dplyr)
    library(GenomicRanges)
    
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


determine_RC_I1 <- function(){
  I1 <- as.character(ShortRead::readFastq(opt$demultiplex_index1ReadsFile)@sread)
  
  d <- select(samples, subject, sample, replicate, index1Seq)
  
  d <- dplyr::bind_rows(lapply(split(d, 1:nrow(d)), function(x){
         x$barcodePercent <- sum(I1 %in% x$index1Seq)/length(I1) * 100
         x$barcodePercentRC <- sum(I1 %in% as.character(Biostrings::reverseComplement(Biostrings::DNAString(x$index1Seq))))/length(I1) * 100
         x
       }))
  
  ifelse(sum(d$barcodePercent) > sum(d$barcodePercentRC), FALSE, TRUE)
}





#' Annotate genomic ranges
#'
#' Create a dataframe of sequencing run ids associated with 1 or more patient identifiers.
#'
#' @param d Data frame containing genomic ranges to be added to UCSC track.
#' @param abundCuts Cut points for estimated abundance (estAbund) values. 
#' @param posColors Color codes for binned abundances (positive integrations).
#' @param negColors Color codes for binned abundances (positive integrations).
#' @param title Track title.
#' @param outputFile Track output file.
#' @param visibility Track default visibility (0 - hide, 1 - dense, 2 - full, 3 - pack, and 4 - squish).
#' @param position Deafult track position.
#' @param padSite Number of NTs to pad sites with for increased visibility.
#' @param siteLabel Text to appear next to sites, ie. 'Patient X, chr12+1052325'.
#' 
#' @return Nothing.
#'
#' @export
createIntUCSCTrack <- function(d, abundCuts = c(5,10,50), 
                               posColors = c("#8C9DFF", "#6768E3", "#4234C7", "#1D00AB"),
                               negColors = c("#FF8C8C", "#E35D5D", "#C72E2E", "#AB0000"),
                               title = 'intSites', outputFile = 'track.ucsc', visibility = 1, 
                               position = 'chr7:127471196-127495720', padSite = 0,
                               siteLabel = NA){
  
  # Check function inputs.
  if(length(posColors) != length(negColors)) 
    stop('The pos and neg color vectors are not the same length.')
  
  if(length(abundCuts) != length(posColors) - 1) 
    stop('The number of aundance cut offs must be one less than the number of provided colors.')
  
  if(! all(c('start', 'end', 'strand', 'seqnames', 'estAbund') %in% names(d))) 
    stop("The expected column names 'start', 'end', 'strand', 'seqnames', 'estAbund' were not found.") 
  
  if(is.na(siteLabel) | ! siteLabel %in% names(d)) 
    stop('The siteLabel parameter is not defined or can not be found in your data.')
  
  
  # Cut the abundance data. Abundance bins will be used to look up color codes.
  # We flank the provided cut break points with 0 and Inf in order to bin all values outside of breaks.
  cuts <- cut(d$estAbund, breaks = c(0, abundCuts, Inf), labels = FALSE)
  
  
  # Convert Hex color codes to RGB color codes. 
  # col2rgb() returns a matrix, here we collapse the columns into comma delimited strings.
  #   grDevices::col2rgb(posColors)
  #         [,1] [,2] [,3] [,4]
  #   red    140  103   66   29
  #   green  157  104   52    0
  #   blue   255  227  199  171
  
  posColors <- apply(grDevices::col2rgb(posColors), 2, paste0, collapse = ',')
  negColors <- apply(grDevices::col2rgb(negColors), 2, paste0, collapse = ',')
  
  
  # Create data fields needed for track table.
  d$score <- 0
  d$color <- ifelse(d$strand == '+', posColors[cuts], negColors[cuts])
  
  # Pad the site n NTs to increase visibility.
  if(padSite > 0){
    d$start <- floor(d$start - padSite/2)
    d$end   <- ceiling(d$end + padSite/2)
  }
  
  # Define track header.
  trackHead <- sprintf("track name='%s' description='%s' itemRgb='On' visibility=%s\nbrowser position %s",
                       title, title, visibility, position)
  
  # Write out track table.
  write(trackHead, file = outputFile, append = FALSE)
  write.table(d[, c('seqnames', 'start', 'end', siteLabel, 'score', 'strand', 'start', 'end', 'color')], 
              sep = '\t', col.names = FALSE, row.names = FALSE, file = outputFile, append = TRUE, quote = FALSE)
}

setOptimalParameters <- function(){
  if(grepl('integrase', opt$mode, ignore.case = TRUE)){
    opt$alignReads_genomeAlignment_anchorRead_maxStartPos <<- 3
    opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned <<- 5
    opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned <<- 5
    opt$prepReads_HMMmatchEnd <<- TRUE
    opt$prepReads_HMMmatchTerminalSeq <<- 'CA'
  } else if(grepl('AAV', opt$mode, ignore.case = TRUE)){
    opt$alignReads_genomeAlignment_anchorRead_maxStartPos <<- 300
    opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned <<- 5
    opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned <<- 10
  } else if(grepl('transposase', opt$mode, ignore.case = TRUE)){
    opt$alignReads_genomeAlignment_anchorRead_maxStartPos <<- 3
    opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned <<- 5
    opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned <<- 5
    opt$prepReads_HMMmatchEnd <<- TRUE
    opt$prepReads_HMMmatchTerminalSeq <<- 'TA'
  }
  else if(grepl('manual', opt$mode, ignore.case = TRUE)){
    # No action
  } else {
    stop('Error -- mode set to an unknown value.')
  }

  if(grepl('quick', opt$mode, ignore.case = TRUE)){
    opt$alignReads_genomeAlignment_repMatch <<- 1000
    opt$buildStdFragments_createMultiHitClusters <<- FALSE
  }
}


: TRUE                                                     # [boolean]   (^) Require a match to the end of the HMM profile. 
prepReads_HMMmatchTerminalSeq: CA 