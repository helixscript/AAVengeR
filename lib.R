tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }

mixAndChunkSeqs <- function(reads, n){
  reads <- reads[order(width(reads))]
  a <- round(length(reads) / 2)
  b <- suppressWarnings(unique(c(rbind(c(1:a),  rev(c((a+1):length(reads)))))))
  reads <- reads[b]
  split(reads, ceiling(seq_along(reads)/n))
}


lpe <- function(x){
  o <- unlist(strsplit(x, '/'))
  o[length(o)]
}


shortRead2DNAstringSet <- function(x){
  r <- x@sread
  names(r) <- sub('\\s+.+$', '', as.character(x@id))
  r
}

waitForFile <- function(f, seconds = 1){
  repeat
  {
    Sys.sleep(seconds)
    if(file.exists(f)) break
  }
  return(TRUE)
}


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

syncReads <-function(...){
  arguments <- list(...)

  # Create a list of read IDs common to all read arguments.
  n <- Reduce(base::intersect, lapply(arguments, names))

  z <- lapply(arguments, function(x){
    x <- x[names(x) %in% n]
    ### x[order(names(x))]
    x[match(n, names(x))]
  })
}


unpackUniqueSampleID <- function(d){
  d$trial     <- unlist(lapply(strsplit(d$uniqueSample, '~'), '[[', 1))
  d$subject   <- unlist(lapply(strsplit(d$uniqueSample, '~'), '[[', 2))
  d$sample    <- unlist(lapply(strsplit(d$uniqueSample, '~'), '[[', 3))
  d$replicate <- unlist(lapply(strsplit(d$uniqueSample, '~'), '[[', 4))
  d
}




# Function to remove minor sequence variation from short sequences.
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




parse_cdhitest_output <- function (file) {
  clusters <- readChar(file, file.info(file)$size)
  clusters <- unlist(base::strsplit(clusters, ">Cluster"))
  clusters <- clusters[2:length(clusters)]
  lapply(clusters, function(x) {
    gsub("\\.\\.\\.", "", unlist(lapply(stringr::str_match_all(x,
                                                      ">([^\\s]+)"), function(y) {
                                                        y[, 2]
                                                      })))
  })
}


loadSamples <- function(){
  samples <- readr::read_tsv(opt$sampleConfigFile, col_types = readr::cols())
  
  if(nrow(samples) == 0){
    write(c(paste(lubridate::now(), 'Error - no lines of information was read from the sample configuration file.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE) 
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
    Ns <- stringr::str_extract(x$adriftRead.linker.seq, 'N+')
    
    # Detecting Ns in the adrift linker and having adriftRead.linkerBarcode.start 
    # set to NA triggers the determination of positions.
    
    if(! is.na(Ns) & is.na(x$adriftRead.linkerBarcode.start)){  
      o <- stringr::str_locate(x$adriftRead.linker.seq, Ns)
      x$adriftRead.linkerBarcode.start  <- 1
      x$adriftRead.linkerBarcode.end    <- o[1,1] - 1
    }
    
    # Detecting Ns in the adrift linker and having adriftRead.linkerRandomID.start 
    # set to NA triggers the determination of positions.
    
    if(! is.na(Ns) & is.na(x$adriftRead.linkerRandomID.start)){ 
      o <- stringr::str_locate(x$adriftRead.linker.seq, Ns)
      x$adriftRead.linkerRandomID.start  <- o[1,1]
      x$adriftRead.linkerRandomID.end    <- o[1,2]
    }
    x
  }))
  
  requiredColumns <- c("trial", "subject", "sample", "replicate",  
                       "adriftRead.linker.seq", "index1.seq", "refGenome.id", "vectorFastaFile", "flags")
  
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


representativeSeq <- function(s, percentReads = 95){
  if(length(s) == 1 | dplyr::n_distinct(s) == 1) return(list(0, s[1]))

  k <- data.frame(table(s))
  s <- unique(s)
  
  if(length(s) > opt$buildStdFragments_representativeSeqCalc_maxReads){
    set.seed(1)
    message('Sampling ', opt$buildStdFragments_representativeSeqCalc_maxReads, ' reads in representativeSeq()')
    s <- sample(s, opt$buildStdFragments_representativeSeqCalc_maxReads)
  }
  
  # Align sequences in order to handle potential indels.
  # Here we create and close connection so that we do not reach an open file limit when called in a loop.
  f <- tmpFile()
  inputFile <- file.path(opt$outputDir, 'tmp', paste0(f, '.fasta'))
  fileConn <- file(inputFile)
  write(paste0('>', paste0('s', 1:length(s)), '\n', s), file = fileConn)
  close(fileConn)
  
  outputFile <- file.path(opt$outputDir, 'tmp', paste0(f, '.representativeSeq.muscle'))

  system(paste(opt$command_muscle, '-quiet -maxiters 1 -diags -in ', inputFile, ' -out ', outputFile))
  
  if(! file.exists(outputFile)) waitForFile(outputFile)

  # Alignments are read in with dashes.
  s <- as.character(Biostrings::readDNAStringSet(outputFile))
  
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


golayCorrection <- function(x){
  library(ShortRead)
  library(dplyr)
  
  f <- file.path(opt$outputDir, 'tmp', paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') )
  writeFasta(x, file = f)
  
  system(paste(opt$command_python2, file.path(opt$softwareDir, 'bin', 'golayCorrection', 'processGolay.py'), f))
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
  library(dplyr)
  library(IRanges)
  library(data.table)
  
  rbindlist(lapply(split(b, b$qname), function(b2){
    
    # Sort BLAST results by query start position and evalue (low to high).
    b2 <- arrange(b2, qstart, evalue)
    b2$strand <- ifelse(b2$send < b2$sstart, '-', '+')
    
    # Alignment to the negative strand will result in the subject end to come before the start.
    # Switch it back so that they are sequential.
    b2$sstart2 <- ifelse(b2$sstart > b2$send, b2$send, b2$sstart)
    b2$send2   <- ifelse(b2$sstart > b2$send, b2$sstart, b2$send)
    
    # Create IRanges
    ir <- IRanges(start = b2$qstart, end = b2$qend)
    if(length(ir) == 0) return(data.frame())
    
    # Name the ranges with the binned query positions followed by the actual subject positions.
    names(ir) <- paste0(b2$qstart, '..', b2$qend, '[', b2$sstart2, b2$strand, b2$send2, ']')
    
    o <- ir[1]
    invisible(lapply(split(ir, 1:length(ir)), function(a){
      if(all(! countOverlaps(o, a, minoverlap = 2) > 0)){
        o <<- c(o, a)
      }
    }))
    
    if(length(o) == 0) return(data.frame())
    
    r <- paste0(unique(names(o)), collapse = ';')
    
    data.frame(qname = b2$qname[1], rearrangement = r)
  }))
}




captureLTRseqsLentiHMM <- function(reads, hmm){
  outputFile <- file.path(opt$outputDir, 'tmp', tmpFile())
  
  reads <- reads[width(reads) > opt$prepReads_HMMsearchReadEndPos]
  if(length(reads) == 0) return(tibble())
  
  writeXStringSet(subseq(reads, opt$prepReads_HMMsearchReadStartPos, opt$prepReads_HMMsearchReadEndPos), outputFile)
  # score 5
  comm <- paste0(opt$command.hmmsearch, ' --max --tblout ', outputFile, '.tbl --domtblout ', outputFile, '.domTbl ', hmm, ' ', outputFile, ' > ', outputFile, '.hmmsearch')
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

  o <- subset(o, targetStart <= opt$prepReads_HMMmaxStartPos & 
                 hmmEnd == hmmLength & 
                 fullScore >= as.numeric(opt$prepReads_HMMminFullBitScore))
 
  if(nrow(o) == 0) return(tibble())
  
  reads2 <- reads[names(reads) %in% o$targetName]
  rm(reads)
  gc()
  reads2 <- reads2[match(o$targetName, names(reads2))]
  
  # Make sure all HMMs alignments result in an CA in the target sequences.
  reads2 <- reads2[as.character(subseq(reads2, o$targetEnd-1, o$targetEnd)) == 'CA']
  if(length(reads2) == 0) return(tibble())
  
  tibble(id = names(reads2),
         LTRname = hmmName,
         LTRseq = as.character(subseq(reads2, 1, o$targetEnd)))
}


blastReads <- function(reads, wordSize = 6, evalue = 10){
  library(Biostrings)
  library(dplyr)

  f <- tmpFile()
  writeXStringSet(reads,  file.path(opt$outputDir, 'tmp', paste0(f, '.fasta')))
  
  system(paste0(opt$command_blastn, ' -word_size ', wordSize, ' -evalue ', evalue,' -outfmt 6 -query ',
                file.path(opt$outputDir, 'tmp', paste0(f, '.fasta')), ' -db ',
                file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd'),
                ' -num_threads 1 -out ', file.path(opt$outputDir, 'tmp', paste0(f, '.blast'))),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  waitForFile(file.path(opt$outputDir, 'tmp', paste0(f, '.blast')))
  
  if(file.info(file.path(opt$outputDir, 'tmp', paste0(f, '.blast')))$size > 0){
    b <- read.table(file.path(opt$outputDir, 'tmp', paste0(f, '.blast')), sep = '\t', header = FALSE)
    names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
    return(b)
  } else {
    return(tibble())
  }
}



blatReads <- function(reads, minIdentity=90, stepSize = 11, tileSize = 11){
  library(Biostrings)
  library(dplyr)
  
  f <- tmpFile()
  writeXStringSet(reads,  file.path(opt$outputDir, 'tmp', paste0(f, '.fasta')))
  
  system(paste0(opt$command_blat, ' ', file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd'), ' ', 
                file.path(opt$outputDir, 'tmp', paste0(f, '.fasta')), ' ',
                file.path(opt$outputDir, 'tmp', paste0(f, '.psl')), ' -minIdentity=', minIdentity, 
                ' -stepSize=', stepSize, ' -tileSize=', tileSize, ' -out=psl -noHead -minScore=0'),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  waitForFile(file.path(opt$outputDir, 'tmp', paste0(f, '.psl')))
  
  if(file.info(file.path(opt$outputDir, 'tmp', paste0(f, '.psl')))$size > 0){
    return(parseBLAToutput(file.path(opt$outputDir, 'tmp', paste0(f, '.psl'))))
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

readSamplePlot <- function(reads, n){
  ds <- sample(unique(reads), n)
  dp <- lapply(strsplit(sort(ds), ''), function(x){ tibble(base = x, n = 1:length(x)) })
  dp <- bind_rows(mapply(function(x, n){ x$read <- n; x}, dp, 1:length(dp), SIMPLIFY = FALSE))
  dp$base <- factor(dp$base, levels = c('A', 'T', 'C', 'G', 'N'))
  
  ggplot(dp, aes(n, read, fill = base)) + theme_bw() + geom_tile() +
    scale_fill_manual(values =  c('red', 'green', 'blue', 'gold', 'gray50')) +
    scale_x_continuous(limits = c(0, width(reads[1])), expand = c(0, 0)) +
    scale_y_continuous(label=comma, limits = c(0, n), expand = c(0, 0)) +
    labs(x = 'Position', y = 'Reads') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}


determine_RC_I1 <- function(){
  I1 <- as.character(ShortRead::readFastq(opt$demultiplex_index1ReadsFile)@sread)
  
  d <- select(samples, subject, sample, replicate, index1.seq)
  
  d <- dplyr::bind_rows(lapply(split(d, 1:nrow(d)), function(x){
         x$barcodePercent <- sum(I1 %in% x$index1.seq)/length(I1) * 100
         x$barcodePercentRC <- sum(I1 %in% as.character(Biostrings::reverseComplement(Biostrings::DNAString(x$index1.seq))))/length(I1) * 100
         x
       }))
  
  ifelse(sum(d$barcodePercent) > sum(d$barcodePercentRC), FALSE, TRUE)
}


# The functions below were originally found in the gintools package
# and have been updated to work with more recent Biconductor packages.
# https://github.com/cnobles/gintools
#--------------------------------------------------------------------------------------------------


#' Refine or resolve sonic break points within a dataset
#'
#' \code{refine_breakpoints} returns a GRanges object where the sonic break
#' point positions have been adjusted based on the dataset and counts.
#'
#' @description Given a GRanges object ....
#'
#' @usage
#' refine_breakpoints(input.sites)
#'
#' refine_breakpoints(input.sites, count = "counts", min.gap = 1L, sata.gap = 3L)
#'
#' @param input.sites GRanges object. Ranges of alignments or sites to adjust.
#' @param counts.col character string. Name of column holding range count
#' information (ie. read counts). If not supplied, assume count of 1 for each
#' row in sites.
#' @param min.gap integer minimum gap to consider combining break points.
#' @param sata.gap integer maximum distance to consider combining break points.
#' @param details logical If TRUE, data columns will be appended to the metadata
#' of the GRanges object noting the original and adjusted position for the site.
#' FALSE by default.
#' @param quiet logical during processing, should messages be suppressed which
#' report findings? TRUE by default.
#'
#' @examples
#' gr <- gintools:::generate_test_granges()
#'
#' refine_breakpoints(gr)
#'
#' @author Christopher L. Nobles, Ph.D.
#'
#' @importFrom magrittr %>%
#'
#' @export
#'

refine_breakpoints <- function(input.sites, counts.col = NULL, min.gap = 1L,
                               sata.gap = 3L, details = FALSE, quiet = TRUE){
  
  stopifnot(class(input.sites) == "GRanges")
  sites <- GenomicRanges::granges(input.sites)
  
  if(!quiet){
    message(paste0("Refining ", length(sites), " break points."))
    message(
      "Generating initial graph by connecting positions with 1 nt difference.")
  }
  
  # Retain original order
  sites$ori.order <- seq_along(sites)
  
  # Identify counts or abundance info, assume 1 if not given, but check for
  # a column named counts first and use, otherwise error out if column is not
  # found.
  
  if(!is.null(counts.col)){
    if(counts.col %in% names(GenomicRanges::mcols(input.sites))){
      counts_pos <- grep(counts.col, names(GenomicRanges::mcols(input.sites)))
      sites$counts <- GenomicRanges::mcols(input.sites)[,counts_pos]
    }else{
      stop("Could not identify 'counts' column.")
    }
  }else{
    if(!quiet) message("Assuming abundance of 1 for each row of sites object.")
    sites$counts <- rep(1, length(sites))
  }
  
  # Reduce the genomic locations of break points down to only unique positions,
  # and identify the abundance of the positions
  
  red_sites <- GenomicRanges::reduce(
    GenomicRanges::flank(sites, -1, start = FALSE),
    min.gapwidth = 0L,
    with.revmap = TRUE)
  red_sites$breakpoint.id <- seq_along(red_sites)
  
  # Summarise count data for reduced site object
  red_counts <- GenomicRanges::mcols(sites) %>%
    as.data.frame(row.names = NULL)
  red_counts <- red_counts[unlist(red_sites$revmap),] %>%
    dplyr::mutate(grp = as.numeric(S4Vectors::Rle(
      values = seq_along(red_sites),
      lengths = S4Vectors::elementNROWS(red_sites$revmap)))) %>%
    dplyr::group_by(grp) %>%
    dplyr::summarise(abund = sum(counts)) %>%
    dplyr::ungroup()
  
  red_sites$abund <- red_counts$abund
  
  # Identify which sites are adjacent to each other
  red_hits <- GenomicRanges::as.data.frame(
    GenomicRanges::findOverlaps(
      red_sites, maxgap = min.gap - 1L, drop.self = TRUE
    )
  )
  
  # Organize data and filter for constructing directed graph
  red_hits <- red_hits %>%
    dplyr::mutate(
      hits.q = queryHits,
      hits.s = subjectHits,
      pos.q = GenomicRanges::start(red_sites[hits.q]),
      pos.s = GenomicRanges::start(red_sites[hits.s]),
      abund.q = red_sites[hits.q]$abund,
      abund.s = red_sites[hits.s]$abund,
      strand = as.vector(GenomicRanges::strand(red_sites[hits.q])),
      is.downstream = ifelse(strand == "+", pos.q > pos.s, pos.q < pos.s),
      keep = abund.q > abund.s,
      keep = ifelse(abund.q == abund.s, is.downstream, keep)) %>%
    dplyr::filter(keep)
  
  # Construct initial graph, the graph will be modified during other parts of
  # the function and can be used to identify clusters as well as 'sources'
  # within those clusters as it's a directed graph. Sources are identified and
  # used to identify the original position given a distribution.
  
  g <- igraph::make_empty_graph(n = length(red_sites), directed = TRUE) %>%
    igraph::add_edges(with(red_hits, vzip(hits.q, hits.s)))
  
  red_sites$clus.id <- igraph::clusters(g)$membership
  
  if(!quiet) message(paste0("Initial cluster count: ", igraph::clusters(g)$no))
  
  # Identify satalite positions that should be included in clusters up to the
  # sata.gap max. This portion of the function tries to reach out to up to the
  # sata.gap distance from the boundry of a cluster to see if there are any
  # further "satalite" positions that have been annotated. It does this by
  # iterively increasing the size from 2nt to the sata.gap by 1nt increments.
  
  if(!quiet){
    message(paste0(
      "Connecting satalite positions up to ", sata.gap, " nt apart."))
  }
  
  null <- lapply(2:sata.gap, function(gap){
    g <<- connect_satalite_vertices(red_sites, g, gap, "downstream")
    red_sites$clus.id <<- igraph::clusters(g)$membership
  })
  
  if(!quiet){
    message("Clusters after satalite connecting: ", igraph::clusters(g)$no)
  }
  
  # During these steps of expanding the graph, it's possible that clusters grew
  # larger than they needed to be. This could be identified by seeing multiple
  # sources within a single cluster. This may arrize by reaching out and
  # connecting the edge of one smaller cluster to a larger one. If this occures,
  # the path between the sources will be 'snipped' between the nodes that share
  # the largest distance appart from one another.
  
  g <- break_connecting_source_paths(red_sites, g, "downstream")
  red_sites$clus.id <- igraph::clusters(g)$membership
  
  if(!quiet){
    message(paste0("Final break point cluster count: ", igraph::clusters(g)$no))
  }
  
  # As cluster membership has been determined, the remaining section of the
  # function serves to consolidate the information and correct the original
  # sites object with the adjusted breakpoints.
  
  src_nodes <- sources(g)
  
  clus_data <- data.frame(
    "clus.id" = seq_along(igraph::clusters(g)$csize),
    "chr" = GenomicRanges::seqnames(red_sites[src_nodes]),
    "strand" = GenomicRanges::strand(red_sites[src_nodes]),
    "breakpoint" = GenomicRanges::start(red_sites[src_nodes]),
    "width" = GenomicRanges::width(unlist(range(
      GenomicRanges::split(red_sites, red_sites$clus.id))))
  )
  
  sites <- sites[unlist(IRanges::as.list(red_sites$revmap))]
  sites$clus.id <- as.numeric(S4Vectors::Rle(
    igraph::clusters(g)$membership,
    S4Vectors::elementNROWS(red_sites$revmap)))
  sites$called.bp <- ifelse(
    GenomicRanges::strand(sites) == "+",
    GenomicRanges::end(sites),
    GenomicRanges::start(sites))
  sites$adj.bp <- clus_data[
    match(sites$clus.id, clus_data$clus.id), "breakpoint"]
  
  if(!quiet){
    message(paste0(
      "Cumulative displacement per range: ",
      round(
        sum( abs(sites$called.bp - sites$adj.bp) ) / length(sites),
        digits = 1)))
  }
  
  GenomicRanges::ranges(sites) <- IRanges::IRanges(
    start = ifelse(
      GenomicRanges::strand(sites) == "+",
      GenomicRanges::start(sites),
      sites$adj.bp),
    end = ifelse(
      GenomicRanges::strand(sites) == "+",
      sites$adj.bp,
      GenomicRanges::end(sites))
  )
  
  # Transfer range information from sites to output_sites object
  sites <- sites[order(sites$ori.order)]
  output_sites <- input.sites
  GenomicRanges::ranges(output_sites) <- GenomicRanges::ranges(sites)
  
  if(details){
    GenomicRanges::mcols(output_sites)$called.pos <- GenomicRanges::mcols(
      sites)$called.bp
    GenomicRanges::mcols(output_sites)$adj.pos <- GenomicRanges::mcols(
      sites)$adj.bp
  }
  
  output_sites
}



#' Standardize genomic site positions within a dataset
#'
#' \code{standardize_sites} returns a GRanges object where the site positions
#' have been standardized with sites within the gap distance.
#'
#' @description Given a GRanges object ...
#'
#' @usage
#' standardize_sites(input.sites)
#'
#' standardize_sites(input.sites, counts.col = "counts", min.gap = 1L, sata.gap = 5L)
#'
#' @param input.sites GRanges object of sites to standardize.
#' @param counts.col character string. Name of column holding range count
#' information (ie. read counts). If not supplied, assume count of 1 for each
#' row in sites.
#' @param min.gap integer minimum gap to consider combine sites.
#' @param sata.gap integer maximum distance to consider combining sites.
#' @param details logical If TRUE, data columns will be appended to the metadata
#' of the GRanges object noting the original and adjusted position for the site.
#' FALSE by default.
#' @param quiet logical during processing, should messages be suppressed which
#' report findings? TRUE by default.
#'
#' @examples
#' gr <- gintools:::generate_test_granges()
#'
#' standardize_sites(gr)
#'
#' @author Christopher Nobles, Ph.D.
#' @export
#'

standardize_sites <- function(input.sites, counts.col = NULL, min.gap = 1L,
                              sata.gap = 5L, details = FALSE, quiet = TRUE){
  stopifnot(class(input.sites) == "GRanges")
  sites <- GenomicRanges::granges(input.sites)
  
  # Retain original order
  sites$ori.order <- seq_along(sites)
  
  # Identify counts or abundance info, assume 1 if not given, but check for
  # a column named counts first and use, otherwise error out if column is not
  # found.
  
  if(!is.null(counts.col)){
    if(counts.col %in% names(GenomicRanges::mcols(input.sites))){
      counts_pos <- grep(counts.col, names(GenomicRanges::mcols(input.sites)))
      sites$counts <- GenomicRanges::mcols(input.sites)[,counts_pos]
    }else{
      stop("Could not identify 'counts' column.")
    }
  }else{
    if(!quiet) message("Assuming abundance of 1 for each row of sites object.")
    sites$counts <- rep(1, length(input.sites))
  }
  
  # Start by reducing the sites object down to only the site's positions
  # and store the revmap for going back to the original sites object.
  if(!quiet){
    message(paste0("Standardizing ", length(sites), " positions."))
    message(
      "Generating initial graph by connecting positions with 1 nt difference.")
  }
  
  # Construct reduced sites object to initialize the processing
  red_sites <- GenomicRanges::reduce(
    GenomicRanges::flank(sites, -1, start = TRUE),
    min.gapwidth = 0L,
    with.revmap = TRUE)
  red_sites$site.id <- seq_along(red_sites)
  
  # Summarise count data for reduced site object
  red_counts <- GenomicRanges::mcols(sites) %>%
    as.data.frame(row.names = NULL)
  red_counts <- red_counts[unlist(red_sites$revmap),] %>%
    dplyr::mutate(grp = as.numeric(S4Vectors::Rle(
      values = seq_along(red_sites),
      lengths = S4Vectors::elementNROWS(red_sites$revmap)))) %>%
    dplyr::group_by(grp) %>%
    dplyr::summarise(abund = sum(counts)) %>%
    dplyr::ungroup()
  
  red_sites$abund <- red_counts$abund
  
  # Identify which sites are adjacent to each other
  red_hits <- GenomicRanges::as.data.frame(
    GenomicRanges::findOverlaps(
      red_sites, maxgap = min.gap - 1L, drop.self = TRUE
    )
  )
  
  # Organize data and filter for constructing directed graph
  red_hits <- red_hits %>%
    dplyr::mutate(
      hits.q = queryHits,
      hits.s = subjectHits,
      pos.q = GenomicRanges::start(red_sites[hits.q]),
      pos.s = GenomicRanges::start(red_sites[hits.s]),
      abund.q = red_sites[hits.q]$abund,
      abund.s = red_sites[hits.s]$abund,
      strand = as.vector(GenomicRanges::strand(red_sites[hits.q])),
      is.upstream = ifelse(strand == "+", pos.q < pos.s, pos.q > pos.s),
      keep = abund.q > abund.s,
      keep = ifelse(abund.q == abund.s, is.upstream, keep)) %>%
    dplyr::filter(keep)
  
  # Construct initial graph, the graph will be modified during other parts of
  # the function and can be used to identify clusters as well as 'sources'
  # within those clusters as it's a directed graph. Sources are identified and
  # used to identify the original position given a distribution.
  g <- igraph::make_empty_graph(n = length(red_sites), directed = TRUE) %>%
    igraph::add_edges(with(red_hits, vzip(hits.q, hits.s)))
  
  red_sites$clus.id <- igraph::clusters(g)$membership
  
  if(!quiet) message(paste0("Initial cluster count: ", igraph::clusters(g)$no))
  
  # When clusters are formed of positions with equivalent abundance, both
  # directional edges are created. This type of cluster does not have sources
  # or sinks. The following set resolves which directed edges to remove, but
  # imposes an upstream bias.
  
  #  g <- resolve_source_deficiency(red_sites, g, bias = "upstream")
  #  red_sites$clus.id <- clusters(g)$membership
  
  # Identify satalite positions that should be included in clusters up to 5nt
  # away. This portion of the function tries to reach out to up to 5nt from the
  # boundry of a cluster to see if there are any further "satalite" positions
  # that have been annotated. It does this by iterively increasing the size
  # from 2nt to 5nt by 1nt increments.
  
  if(!quiet){
    message(paste0(
      "Connecting satalite positions up to ", sata.gap, " nt apart."))
  }
  
  null <- lapply(2:sata.gap, function(gap){
    g <<- connect_satalite_vertices(red_sites, g, gap, "upstream")
    red_sites$clus.id <<- igraph::clusters(g)$membership
  })
  
  if(!quiet){
    message("Clusters after satalite connecting: ", igraph::clusters(g)$no)
  }
  
  # During these steps of expanding the graph, it's possible that clusters grew
  # larger than they needed to be. This could be identified by seeing multiple
  # sources within a single cluster. This may arrize by reaching out and
  # connecting the edge of one smaller cluster to a larger one. If this occures,
  # the path between the sources will be 'snipped' between the nodes that share
  # the largest distance appart from one another.
  
  g <- break_connecting_source_paths(red_sites, g, "upstream")
  red_sites$clus.id <- igraph::clusters(g)$membership
  
  if(!quiet){
    message(paste0("Clusters after clipping: ", igraph::clusters(g)$no))
  }
  
  # In the end, sources that are within the range of satalites (5L default),
  # should be grouped together. These sources should be connected by an edge,
  # pointing toward the larger source. The chosen source will have the highest
  # abundance of widths / sonic breaks / fragment lengths.
  
  if(!quiet){
    message("Connecting clusters with source nodes within ", sata.gap, " nt.")
  }
  
  g <- connect_adjacent_clusters(red_sites, g, gap = sata.gap, "upstream")
  red_sites$clus.id <- igraph::clusters(g)$membership
  
  if(!quiet) message(paste0("Final cluster count: ", igraph::clusters(g)$no))
  
  # If wide clusters are generated, these can confound the cluster source. For
  # this reason, the "true" source for the cluster will be resolved by
  # picking the source with the greatest number of fragment lengths, where ties
  # are decided by randon picking.
  
  g <- resolve_cluster_sources(red_sites, g, "upstream")
  red_sites$clus.id <- igraph::clusters(g)$membership
  
  # As cluster membership has been determined, the remaining section of the
  # function serves to consolidate the information and adjust the original
  # sites object with the standardized positions.
  
  src.nodes <- sources(g)
  
  clus.data <- data.frame(
    "clus.id" = seq_along(igraph::clusters(g)$csize),
    "chr" = GenomicRanges::seqnames(red_sites[src.nodes]),
    "strand" = GenomicRanges::strand(red_sites[src.nodes]),
    "position" = GenomicRanges::start(red_sites[src.nodes]),
    "width" = GenomicRanges::width(unlist(range(
      GenomicRanges::split(red_sites, red_sites$clus.id))))
  )
  
  # Index the sites object and add adjusted data to metadata
  sites <- sites[unlist(IRanges::as.list(red_sites$revmap))]
  sites$clus.id <- as.numeric(S4Vectors::Rle(
    igraph::clusters(g)$membership,
    S4Vectors::elementNROWS(red_sites$revmap)))
  sites$called.pos <- ifelse(
    GenomicRanges::strand(sites) == "+",
    GenomicRanges::start(sites),
    GenomicRanges::end(sites))
  sites$adj.pos <- clus.data[
    match(sites$clus.id, clus.data$clus.id), "position"]
  
  if(!quiet){
    message(paste0(
      "Cumulative displacement per range: ",
      sum(abs(sites$called.pos - sites$adj.pos))/length(sites)))
  }
  
  # Adjust range information within sites object
  GenomicRanges::ranges(sites) <- IRanges::IRanges(
    start = ifelse(
      GenomicRanges::strand(sites) == "+",
      sites$adj.pos,
      GenomicRanges::start(sites)),
    end = ifelse(
      GenomicRanges::strand(sites) == "+",
      GenomicRanges::end(sites),
      sites$adj.pos)
  )
  
  # Transfer range information from sites to output_sites object
  sites <- sites[order(sites$ori.order)]
  output_sites <- input.sites
  GenomicRanges::ranges(output_sites) <- GenomicRanges::ranges(sites)
  
  if(details){
    GenomicRanges::mcols(output_sites)$called.pos <- GenomicRanges::mcols(
      sites)$called.pos
    GenomicRanges::mcols(output_sites)$adj.pos <- GenomicRanges::mcols(
      sites)$adj.pos
  }
  
  output_sites
}


#' Zip vectors together into a single vector.
#'
#' \code{vzip} returns a single vector from input vectors, in order of input and
#' index in each vector.
#'
#' @description Similar to python zip, `vzip` takes input vectors and merges
#' them together by their input order and index. A simple example is two numeric
#' vectors, A = c(1,1,1) and B = c(2,2,2). The output of vzip(A,B) would simply
#' be a single vector of c(1,2,1,2,1,2). Any number of vectors can be input, but
#' each input vector must be of the same length. Output vector class depends on
#' input vector consensus.
#'
#' @usage
#' vzip(...)
#'
#' @param ... any number of vectors to zip together. Each vector must be of
#' equal length.
#'
#' @examples
#' A <- c(1,0,1)
#' B <- c(6,7,4)
#' vzip(A, B)
#'
#' @author Christopher Nobles, Ph.D.
#' @export
#'

vzip <- function(...){
  v <- list(...)
  if(length(unique(sapply(v, length))) > 1){
    stop("All input vectors are not of equal length.")
  }
  as.vector(t(matrix(unlist(v), ncol = length(v))))
}


#' Connect integration sites which lie a specific gap distance away from a
#' cluster of integration sites.
#'
#' \code{connect_satalite_vertices} returns a graph where nodes within 'gap'
#' distance from clusters are now connected to the edge of the clusters.
#'
#' @description Given a list of unique integration site positions (reduced
#' GRanges object) and a directed graph of connected components, this function
#' identifies each node within gap range of clusters and creates an edge to
#' connect the cluster to the 'satalite' node. Edges are drawn from the last
#' node in the cluster to the 'satalite' node, but directionality is determined
#' first by abundance and secondly by an upstream bias.
#'
#' @usage
#' connect_satalite_vertices(red.sites, graph, gap, bais)
#'
#' @param red.sites GRanges object which has been reduced to single nt positions
#' and contains the revmap from the original GRanges object. The object must
#' also contain a column for cluster membership (clus.id) and a column for
#' abundance (fragLengths).
#' @param graph a directed graph built from the red.sites object. Each node
#' corresponds to a row in the red.sites object.
#' @param gap integer nucleotide (nt) gap distance to consider joining to
#' clusters.
#' @param bias either "upsteam" or "downstream", designating which position to
#' choose if other decision metrics are tied.
#'
#' @examples
#' gr <- gintools:::generate_test_granges(stdev = 3)
#' red.sites <- reduce(
#'   flank(gr, -1, start = TRUE),
#'   min.gapwidth = 0L,
#'   with.revmap = TRUE)
#' red.sites$siteID <- seq_along(red.sites)
#' revmap <- as.list(red.sites$revmap)
#' red.sites$fragLengths <- lengths(revmap)
#' red.hits <- GenomicRanges::as.data.frame(
#'   findOverlaps(red.sites, maxgap = 0L, drop.self = TRUE))
#' red.hits <- red.hits %>%
#'   mutate(q_pos = start(red.sites[queryHits])) %>%
#'   mutate(s_pos = start(red.sites[subjectHits])) %>%
#'   mutate(q_fragLengths = red.sites[queryHits]$fragLengths) %>%
#'   mutate(s_fragLengths = red.sites[subjectHits]$fragLengths) %>%
#'   mutate(strand = unique(strand(
#'     c(red.sites[queryHits], red.sites[subjectHits])))) %>%
#'   mutate(is.upstream = ifelse(
#'     strand == "+",
#'     q_pos < s_pos,
#'     q_pos > s_pos)) %>%
#'   mutate(keep = q_fragLengths > s_fragLengths) %>%
#'   mutate(keep = ifelse(
#'     q_fragLengths == s_fragLengths,
#'     is.upstream,
#'     keep)) %>%
#'   filter(keep)
#' g <- make_empty_graph(n = length(red.sites), directed = TRUE) %>%
#'   add_edges(unlist(mapply(
#'     c, red.hits$queryHits, red.hits$subjectHits, SIMPLIFY = FALSE)))
#' red.sites$clus.id <- clusters(g)$membership
#'
#' connect_satalite_vertices(red.sites, g, gap = 2L, "upstream")
#'
#' @author Christopher Nobles, Ph.D.
#'
#' @importFrom magrittr %>%

connect_satalite_vertices <- function(red.sites, graph, gap, bias){
  clus_mem <- igraph::clusters(graph)$membership
  
  clus_ranges <- unlist(GenomicRanges::reduce(
    GenomicRanges::split(red.sites, clus_mem),
    min.gapwidth = (gap-1)
  ))
  
  sata_hits <- as.data.frame(
    GenomicRanges::findOverlaps(
      clus_ranges, maxgap = gap - 1L, drop.self = TRUE
    )
  )
  
  names(sata_hits) <- c("source.clus", "sata.clus")
  
  red_df <- GenomicRanges::as.data.frame(red.sites)
  
  if(nrow(sata_hits) > 0){
    clus_data <- red_df %>%
      dplyr::group_by(clus.id) %>%
      dplyr::summarize(
        clus.pos.mean = as.integer(mean(start)),
        min.abund = min(abund),
        sum.abund = sum(abund))
    
    if(bias == "upstream"){
      sata_hits <- sata_hits %>%
        dplyr::mutate(
          source.pos = clus_data[source.clus,]$clus.pos.mean,
          sata.pos = clus_data[sata.clus,]$clus.pos.mean,
          min.src.abund = clus_data[.$source.clus,]$min.abund,
          min.sat.abund = clus_data[.$sata.clus,]$min.abund,
          sum.src.abund = clus_data[.$source.clus,]$sum.abund,
          sum.sat.abund = clus_data[.$sata.clus,]$sum.abund,
          is.upstream = source.pos < sata.pos) %>%
        dplyr::filter(
          as.integer(min.src.abund) >= as.integer(min.sat.abund),
          sum.src.abund > sum.sat.abund,
          abs(source.clus - sata.clus) == 1)
    }else if(bias == "downstream"){
      sata_hits <- sata_hits %>%
        dplyr::mutate(
          source.pos = clus_data[source.clus,]$clus.pos.mean,
          sata.pos = clus_data[sata.clus,]$clus.pos.mean,
          min.src.abund = clus_data[.$source.clus,]$min.abund,
          min.sat.abund = clus_data[.$sata.clus,]$min.abund,
          sum.src.abund = clus_data[.$source.clus,]$sum.abund,
          sum.sat.abund = clus_data[.$sata.clus,]$sum.abund,
          is.downstream = source.pos > sata.pos) %>%
        dplyr::filter(
          as.integer(min.src.abund) >= as.integer(min.sat.abund),
          sum.src.abund > sum.sat.abund,
          abs(source.clus - sata.clus) == 1)
    }else{
      stop("No bias specified. Please choose either 'upstream' or 'downstream'.")
    }
    
    if(nrow(sata_hits) > 0){
      clus.map <- GenomicRanges::findOverlaps(clus_ranges, red.sites)
      clus.list <- split(
        S4Vectors::subjectHits(clus.map), S4Vectors::queryHits(clus.map))
      
      if(bias == "upstream"){
        sata_hits <- sata_hits %>%
          dplyr::mutate(
            source.node = ifelse(
              sata_hits$is.upstream,
              sapply(clus.list[sata_hits$source.clus], dplyr::last),
              sapply(clus.list[sata_hits$source.clus], dplyr::first)),
            sata.node = ifelse(
              is.upstream,
              sapply(clus.list[sata.clus], dplyr::first),
              sapply(clus.list[sata.clus], dplyr::last)))
      }else if(bias == "downstream"){
        sata_hits <- sata_hits %>%
          dplyr::mutate(
            source.node = ifelse(
              sata_hits$is.downstream,
              sapply(clus.list[sata_hits$source.clus], dplyr::first),
              sapply(clus.list[sata_hits$source.clus], dplyr::last)),
            sata.node = ifelse(
              is.downstream,
              sapply(clus.list[sata.clus], dplyr::last),
              sapply(clus.list[sata.clus], dplyr::first)))
      }
      
      sata.edges <- with(sata_hits, vzip(source.node, sata.node))
    }else{
      sata.edges <- c()
    }
  }else{
    sata.edges <- c()
  }
  igraph::add_edges(graph, sata.edges)
}


#' Break graph paths which connect sources.
#'
#' \code{break_connecting_source_paths} returns a graph where only one source
#' is present per cluster.
#'
#' @description Given a list of unique integration site positions (reduced
#' GRanges object) and a directed graph of connected components, this function
#' identifies clusters with multiple sources, the paths between those sources,
#' and removes edges along the path so that each cluster only has one source
#' node. Edge removal is first based on nucleotide distance (greater distance
#' prefered), then based on abundance (lowest abundance prefered), then on an
#' upstream bias (downstream connection will be removed when everything ties).
#'
#' @usage
#' break_connecting_source_paths(red.sites, graph, bias)
#'
#' @param red.sites GRanges object which has been reduced to single nt positions
#' and contains the revmap from the original GRanges object. The object must
#' also contain a column for cluster membership (clusID) and a column for
#' abundance (fragLengths).
#'
#' @param graph a directed graph built from the red.sites object. Each node
#' corresponds to a row in the red.sites object.
#'
#' @param bias either "upsteam" or "downstream", designating which position to
#' choose if other decision metrics are tied.
#'
#' @examples
#' gr <- generate_test_granges(stdev = 3)
#' red.sites <- reduce(
#'   flank(gr, -1, start = TRUE),
#'   min.gapwidth = 0L,
#'   with.revmap = TRUE)
#' red.sites$siteID <- seq(1:length(red.sites))
#' revmap <- as.list(red.sites$revmap)
#' red.sites$abundance <- sapply(revmap, length)
#' red.hits <- GenomicRanges::as.data.frame(
#'   findOverlaps(red.sites, maxgap = 0L, drop.self = TRUE))
#' red.hits <- red.hits %>%
#'   mutate(q_pos = start(red.sites[queryHits])) %>%
#'   mutate(s_pos = start(red.sites[subjectHits])) %>%
#'   mutate(q_abund = red.sites[queryHits]$abundance) %>%
#'   mutate(s_abund = red.sites[subjectHits]$abundance) %>%
#'   mutate(strand = unique(strand(
#'     c(red.sites[queryHits], red.sites[subjectHits])))) %>%
#'   mutate(is.upstream = ifelse(
#'     strand == "+",
#'     q_pos < s_pos,
#'     q_pos > s_pos)) %>%
#'   mutate(keep = q_abund > s_abund) %>%
#'   mutate(keep = ifelse(
#'     q_abund == s_abund,
#'     is.upstream,
#'     keep)) %>%
#'   filter(keep)
#' g <- make_empty_graph(n = length(red.sites), directed = TRUE) %>%
#'   add_edges(unlist(mapply(
#'     c, red.hits$queryHits, red.hits$subjectHits, SIMPLIFY = FALSE)))
#' red.sites$clusID <- clusters(g)$membership
#' g <- connect_satalite_vertices(red.sites, g, gap = 2L, bias = "upstream")
#' red.sites$clusID <- clusters(g)$membership
#' break_connecting_source_paths(red.sites, g, "upstream")
#'
#' @author Christopher Nobles, Ph.D.
#'
#' @importFrom magrittr %>%
#'

break_connecting_source_paths <- function(red.sites, graph, bias){
  src_nodes <- sources(graph)
  sources_p_clus <- IRanges::IntegerList(split(
    src_nodes, igraph::clusters(graph)$membership[src_nodes]))
  clus_w_multi_sources <- sources_p_clus[S4Vectors::elementNROWS(sources_p_clus) > 1]
  
  if(length(clus_w_multi_sources) > 0){
    adj_pairs <- do.call(c, lapply(clus_w_multi_sources, function(x){
      lapply(1:(length(x)-1), function(i) c(x[i], x[i+1]))
    }))
    
    snk_nodes <- sinks(graph)
    
    edges_to_edit <- data.frame(
      src.node.i = unlist(adj_pairs)[
        IRanges::start(IRanges::IntegerList(adj_pairs)@partitioning)],
      src.node.j = unlist(adj_pairs)[
        IRanges::end(IRanges::IntegerList(adj_pairs)@partitioning)]) %>%
      dplyr::mutate(
        src.node.i.abund = as.numeric(red.sites[src.node.i]$abund),
        src.node.j.abund = as.numeric(red.sites[src.node.j]$abund),
        sink.node = IRanges::start(
          IRanges::findOverlapPairs(
            IRanges::IRanges(src.node.i, src.node.j),
            IRanges::IRanges(snk_nodes, width = 1))@second))
    
    # Identify the nodes adjacent to sinks between connected sources
    # then filter adjacent pairs to identify which edge should be 'clipped'.
    # Filtering based first on adjacent node distance (edges with greater
    # distance get clipped), then abundance (lower abund gets clipped), then
    # biasing on upstream edges over downstream (downstream is clipped for
    # tie breaking).
    
    if(bias == "upstream"){
      target_edges <- dplyr::bind_rows(lapply(
        seq_len(nrow(edges_to_edit)), function(i){
          sink <- edges_to_edit[i, "sink.node"]
          path <- unlist(igraph::all_simple_paths(
            igraph::as.undirected(graph),
            edges_to_edit[i, "src.node.i"],
            edges_to_edit[i, "src.node.j"]))
          pos <- which(path == sink)
          data.frame(
            sink = rep(sink, 2),
            adj.node = c(path[pos-1], path[pos+1])) })) %>%
        dplyr::mutate(
          sink.pos = GenomicRanges::start(red.sites[sink]),
          adj.pos = GenomicRanges::start(red.sites[adj.node]),
          adj.abund = red.sites[adj.node]$abund,
          nt.dist = abs(sink.pos - adj.pos),
          strand = as.character(GenomicRanges::strand(red.sites[sink])),
          is.upstream = ifelse(
            strand == "+", sink.pos < adj.pos, sink.pos > adj.pos)) %>%
        dplyr::group_by(sink) %>%
        dplyr::filter(nt.dist == max(nt.dist)) %>%
        dplyr::filter(adj.abund == min(adj.abund)) %>%
        dplyr::mutate(
          group_size = n(),
          keep = ifelse(group_size == 1, TRUE, !is.upstream)) %>%
        dplyr::filter(keep) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    }else if(bias == "downstream"){
      target_edges <- dplyr::bind_rows(lapply(
        seq_len(nrow(edges_to_edit)), function(i){
          sink <- edges_to_edit[i, "sink.node"]
          path <- unlist(igraph::all_simple_paths(
            igraph::as.undirected(graph),
            edges_to_edit[i, "src.node.i"],
            edges_to_edit[i, "src.node.j"]))
          pos <- which(path == sink)
          data.frame(
            sink = rep(sink, 2),
            adj.node = c(path[pos-1], path[pos+1])) })) %>%
        dplyr::mutate(
          sink.pos = GenomicRanges::start(red.sites[sink]),
          adj.pos = GenomicRanges::start(red.sites[adj.node]),
          adj.abund = red.sites[adj.node]$abund,
          nt.dist = abs(sink.pos - adj.pos),
          strand = as.character(
            GenomicRanges::strand(red.sites[sink])),
          is.downstream = ifelse(
            strand == "+", sink.pos > adj.pos, sink.pos < adj.pos)) %>%
        dplyr::group_by(sink) %>%
        dplyr::filter(nt.dist == max(nt.dist)) %>%
        dplyr::filter(adj.abund == min(adj.abund)) %>%
        dplyr::mutate(
          group_size = n(),
          keep = ifelse(group_size == 1, TRUE, !is.downstream)) %>%
        dplyr::filter(keep) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    }else{
      stop("No bias specified. Please choose either 'upstream' or 'downstream'.")
    }
    
    break_edges <- with(target_edges, vzip(sink, adj.node))
    
    edge_ids_to_break <- igraph::get.edge.ids(
      graph, break_edges, directed = FALSE)
  }else{
    edge_ids_to_break <- c()
  }
  
  igraph::delete_edges(graph, edge_ids_to_break)
}




#' Identify sources in a directed graph
#'
#' \code{sources} returns a numerical vector of source nodes.
#'
#' @description From a directed graph, this function returns all source nodes, or
#' nodes which only act as head nodes and not as tail nodes.
#'
#' @usage
#' sources(graph)
#'
#' @param graph a directed graph (igraph).
#'
#' @examples
#' g <- make_graph(edges = c(1,2, 2,3, 1,3), directed = TRUE)
#' sources(g)
#'
#' @author Christopher Nobles, Ph.D.
#'

sources <- function(graph){
  if(!igraph::is_directed(graph)){
    message("Graph provided is not a directed graph.")
    srcs <- c()
  }else{
    srcs <- which(Matrix::colSums(
      igraph::get.adjacency(graph, sparse = TRUE)) == 0)
  }
  srcs
}




#' Connect adjacent clusters when their sources are within a specified distance.
#'
#' \code{connect_adjacent_clusters} returns a graph where adjacent clusters with
#' sources within the gap distance are joined.
#'
#' @description Given a list of unique integration site positions (reduced
#' GRanges object) and a directed graph of connected components, this function
#' identifies clusters where source nodes are within a gap distance away from
#' each other, and connects the source nodes with a directional edge, based
#' first on abundance, and secondly on an upstream bias for tie breaking.
#'
#' @usage
#' connect_adjacent_clusters(red.sites, graph, gap, bais)
#'
#' @param red.sites GRanges object which has been reduced to single nt positions
#' and contains the revmap from the original GRanges object. The object must
#' also contain a column for cluster membership (clusID) and a column for
#' abundance (fragLengths).
#'
#' @param graph a directed graph built from the red.sites object. Each node
#' corresponds to a row in the red.sites object.
#'
#' @param bias either "upsteam" or "downstream", designating which position to
#' choose if other decision metrics are tied.
#'
#' @examples
#' gr <- gintools:::generate_test_granges(stdev = 3)
#' red.sites <- reduce(
#'   flank(gr, -1, start = TRUE),
#'   min.gapwidth = 0L,
#'   with.revmap = TRUE)
#' red.sites$siteID <- seq_along(red.sites)
#' revmap <- as.list(red.sites$revmap)
#' red.sites$fragLengths <- lengths(revmap)
#' red.hits <- GenomicRanges::as.data.frame(
#'   findOverlaps(red.sites, maxgap = 0L, drop.self = TRUE))
#' red.hits <- red.hits %>%
#'   mutate(q_pos = start(red.sites[queryHits])) %>%
#'   mutate(s_pos = start(red.sites[subjectHits])) %>%
#'   mutate(q_abund = red.sites[queryHits]$abundance) %>%
#'   mutate(s_abund = red.sites[subjectHits]$abundance) %>%
#'   mutate(strand = unique(strand(
#'     c(red.sites[queryHits], red.sites[subjectHits])))) %>%
#'   mutate(is.upstream = ifelse(
#'     strand == "+",
#'     q_pos < s_pos,
#'     q_pos > s_pos)) %>%
#'   mutate(keep = q_abund > s_abund) %>%
#'   mutate(keep = ifelse(
#'     q_abund == s_abund,
#'     is.upstream,
#'     keep)) %>%
#'   filter(keep)
#' g <- make_empty_graph(n = length(red.sites), directed = TRUE) %>%
#'   add_edges(unlist(mapply(
#'     c, red.hits$queryHits, red.hits$subjectHits, SIMPLIFY = FALSE)))
#' red.sites$clusID <- clusters(g)$membership
#' g <- connect_satalite_vertices(red.sites, g, gap = 2L, "upstream")
#' red.sites$clusID <- clusters(g)$membership
#' g <- break_connecting_source_paths(red.sites, g, "upstream")
#' red.sites$clusID <- clusters(g)$membership
#'
#' connect_adjacent_clusters(red.sites, g, gap = 5L, "upstream")
#'
#' @author Christopher Nobles, Ph.D.
#'
#' @importFrom magrittr %>%
#'

connect_adjacent_clusters <- function(red.sites, graph, gap, bias){
  src_nodes <- sources(graph)
  near_sources <- GenomicRanges::findOverlaps(
    red.sites[src_nodes],
    maxgap = gap - 1L,
    drop.self = TRUE
  )
  
  if(length(near_sources) > 0){
    # Identify sources of clusters within the largets satalite gap distance
    # and identify the directionality of the edge to create based first on
    # abundance (source will likely have greater abundance) and then by
    # upstream bias (more likely the origin site is upstream of drifting mapped
    # reads).
    if(bias == "upstream"){
      near_src_df <- data.frame(
        node.i = src_nodes[S4Vectors::queryHits(near_sources)],
        node.j = src_nodes[S4Vectors::subjectHits(near_sources)]) %>%
        dplyr::mutate(
          abund.i = red.sites[node.i]$abund,
          abund.j = red.sites[node.j]$abund,
          pos.i = GenomicRanges::start(red.sites[node.i]),
          pos.j = GenomicRanges::start(red.sites[node.j]),
          strand = as.character(GenomicRanges::strand(red.sites[node.i])),
          is.upstream = ifelse(strand == "+", pos.i < pos.j, pos.i > pos.j))
    }else if(bias == "downstream"){
      near_src_df <- data.frame(
        node.i = src_nodes[S4Vectors::queryHits(near_sources)],
        node.j = src_nodes[S4Vectors::subjectHits(near_sources)]) %>%
        dplyr::mutate(
          abund.i = red.sites[node.i]$abund,
          abund.j = red.sites[node.j]$abund,
          pos.i = GenomicRanges::start(red.sites[node.i]),
          pos.j = GenomicRanges::start(red.sites[node.j]),
          strand = as.character(GenomicRanges::strand(red.sites[node.i])),
          is.downstream = ifelse(strand == "+", pos.i > pos.j, pos.i < pos.j))
    }else{
      stop("No bias specified. Please choose either 'upstream' or 'downstream'.")
    }
    
    redundant_graph <- igraph::make_graph(
      edges = with(near_src_df, vzip(
        1:nrow(near_src_df),
        match(paste(node.i, node.j), paste(node.j, node.i)))))
    redundant_groups <- igraph::clusters(redundant_graph)$membership
    
    if(bias == "upstream"){
      near_src_df <- dplyr::mutate(
        near_src_df, redundant.grp = redundant_groups) %>%
        dplyr::filter(abund.i >= abund.j) %>%
        dplyr::group_by(redundant.grp) %>%
        dplyr::mutate(
          group_size = n(),
          keep = ifelse(group_size == 1, TRUE, is.upstream)) %>%
        dplyr::filter(keep)
    }else if(bias == "downstream"){
      near_src_df <- dplyr::mutate(
        near_src_df, redundant.grp = redundant_groups) %>%
        dplyr::filter(abund.i >= abund.j) %>%
        dplyr::group_by(redundant.grp) %>%
        dplyr::mutate(
          group_size = n(),
          keep = ifelse(group_size == 1, TRUE, is.downstream)) %>%
        dplyr::filter(keep)
    }
    
    edges_to_connect_near_srcs <- with(near_src_df, vzip(node.i, node.j))
  }else{
    edges_to_connect_near_srcs <- c()
  }
  
  igraph::add.edges(graph, edges_to_connect_near_srcs)
}




#' Resolve primary sources from clusters with multiple souce nodes.
#'
#' \code{resolve_cluster_sources} returns a graph where each cluster only
#' has a single primary source node.
#'
#' @description Given a list of unique integration site positions (reduced
#' GRanges object) and a directed graph of connected components, this function
#' identifies clusters where multiple source nodes exist and then identifies
#' which source should be considered the primary source node, first based on
#' abundance and then
#'
#' @usage
#' resolve_cluster_sources(red.sites, graph)
#'
#' @param red.sites GRanges object which has been reduced to single nt positions
#' and contains the revmap from the original GRanges object. The object must
#' also contain a column for cluster membership (clusID) and a column for
#' abundance (fragLengths).
#'
#' @param graph a directed graph built from the red.sites object. Each node
#' corresponds to a row in the red.sites object.
#'
#' @param bias either "upsteam" or "downstream", designating which position to
#' choose if other decision metrics are tied.
#'
#' @examples
#' gr <- gintools:::generate_test_granges(stdev = 3)
#' red.sites <- reduce(
#'   flank(gr, -1, start = TRUE),
#'   min.gapwidth = 0L,
#'   with.revmap = TRUE)
#' red.sites$siteID <- seq_along(red.sites)
#' revmap <- as.list(red.sites$revmap)
#' red.sites$abundance <- lengths(revmap)
#' red.hits <- GenomicRanges::as.data.frame(
#'   findOverlaps(red.sites, maxgap = 0L, drop.self = TRUE))
#' red.hits <- red.hits %>%
#'   mutate(q_pos = start(red.sites[queryHits])) %>%
#'   mutate(s_pos = start(red.sites[subjectHits])) %>%
#'   mutate(q_abund = red.sites[queryHits]$abundance) %>%
#'   mutate(s_abund = red.sites[subjectHits]$abundance) %>%
#'   mutate(strand = unique(strand(
#'     c(red.sites[queryHits], red.sites[subjectHits])))) %>%
#'   mutate(is.upstream = ifelse(
#'     strand == "+",
#'     q_pos < s_pos,
#'     q_pos > s_pos)) %>%
#'   mutate(keep = q_abund > s_abund) %>%
#'   mutate(keep = ifelse(
#'     q_abund == s_abund,
#'     is.upstream,
#'     keep)) %>%
#'   filter(keep)
#' g <- make_empty_graph(n = length(red.sites), directed = TRUE) %>%
#'   add_edges(unlist(mapply(
#'     c, red.hits$queryHits, red.hits$subjectHits, SIMPLIFY = FALSE)))
#' red.sites$clusID <- clusters(g)$membership
#' g <- connect_satalite_vertices(red.sites, g, gap = 2L, "upstream")
#' red.sites$clusID <- clusters(g)$membership
#' g <- break_connecting_source_paths(red.sites, g, "upstream")
#' red.sites$clusID <- clusters(g)$membership
#' g <- connect_adjacent_clusters(red.sites, g, gap = 5L, "upstream")
#' red.sites$clusID <- clusters(g)$membership
#'
#' resolve_cluster_sources(red.sites, g, "upstream")
#'
#' @author Christopher Nobles, Ph.D.
#'
#' @importFrom magrittr %>%

resolve_cluster_sources <- function(red.sites, graph, bias){
  src_nodes <- sources(graph)
  sources_p_clus <- IRanges::IntegerList(split(
    src_nodes, igraph::clusters(graph)$membership[src_nodes]))
  clus_w_multi_sources <- sources_p_clus[S4Vectors::elementNROWS(sources_p_clus) > 1]
  
  if(length(clus_w_multi_sources) > 0){
    if(bias == "upstream"){
      resolve_df <- data.frame(
        node = unlist(clus_w_multi_sources),
        clus = as.numeric(S4Vectors::Rle(
          values = seq_along(clus_w_multi_sources),
          lengths = S4Vectors::elementNROWS(clus_w_multi_sources)))) %>%
        dplyr::mutate(abund = red.sites[node]$abund) %>%
        dplyr::group_by(clus) %>%
        dplyr::mutate(
          top.abund = abund == max(abund),
          strand = as.character(GenomicRanges::strand(red.sites[node])),
          pos = GenomicRanges::start(red.sites[node])) %>%
        dplyr::group_by(clus, top.abund) %>%
        dplyr::mutate(
          grp.size = n(),
          is.upstream = ifelse(strand == "+", pos == min(pos), pos == max(pos)),
          src = ifelse(
            top.abund == TRUE & grp.size > 1, is.upstream, top.abund)) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    }else if(bias == "downstream"){
      resolve_df <- data.frame(
        node = unlist(clus_w_multi_sources),
        clus = as.numeric(S4Vectors::Rle(
          values = seq_along(clus_w_multi_sources),
          lengths = S4Vectors::elementNROWS(clus_w_multi_sources)))) %>%
        dplyr::mutate(
          abund = red.sites[node]$abund) %>%
        dplyr::group_by(clus) %>%
        dplyr::mutate(
          top.abund = abund == max(abund),
          strand = as.character(GenomicRanges::strand(red.sites[node])),
          pos = GenomicRanges::start(red.sites[node])) %>%
        dplyr::group_by(clus, top.abund) %>%
        dplyr::mutate(
          grp.size = n(),
          is.downstream = ifelse(
            strand == "+", pos == max(pos), pos == min(pos)),
          src = ifelse(
            top.abund == TRUE & grp.size > 1, is.downstream, top.abund)) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    }else{
      stop("No bias specified. Please choose either 'upstream' or 'downstream'.")
    }
    
    src_nodes <- resolve_df[resolve_df$src == TRUE,]$node
    snk_nodes <- lapply(seq_along(clus_w_multi_sources), function(i){
      resolve_df[resolve_df$clus == i & resolve_df$src == FALSE,]$node
    })
    
    # Accomidates multiple sinks for singular sources in a cluster.
    resolve_edges <- do.call(c, lapply(seq_along(src_nodes), function(i){
      src <- src_nodes[i]
      snk <- snk_nodes[[i]]
      vzip(rep(src, length(snk)), snk) }))
  }else{
    resolve_edges <- c()
  }
  
  igraph::add.edges(graph, resolve_edges)
}


#' Identify sinks in a directed graph
#'
#' \code{sinks} returns a numerical vector of sink nodes.
#'
#' @description From a directed graph, this function returns all sink nodes, or
#' nodes which only act as tail nodes and not as head nodes.
#'
#' @usage
#' sinks(graph)
#'
#' @param graph a directed graph (igraph).
#'
#' @examples
#' g <- make_graph(edges = c(1,2, 2,3, 1,3), directed = TRUE)
#' sinks(g)
#'
#' @author Christopher Nobles, Ph.D.
#'

sinks <- function(graph){
  if(!igraph::is_directed(graph)){
    message("Graph provided is not a directed graph.")
    snks <- c()
  }else{
    snks <- which(Matrix::rowSums(
      igraph::get.adjacency(graph, sparse = TRUE)) == 0)
  }
  snks
}





