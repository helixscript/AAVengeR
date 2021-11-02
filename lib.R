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
    if(file.exists(f)) break
    Sys.sleep(seconds)
  }
  return(TRUE)
}


qualTrimReads <- function(f, chunkSize, label, ouputDir){
  
  # Create a pointer like object to the read file.
  strm <- FastqStreamer(f, n = as.integer(chunkSize))
  
  # Extract chunks from file and write out chunks with numeric suffixes, eg. anchorReads.5
  n <- 1
  repeat {
    fq <- yield(strm)
    if(length(fq) == 0) break
    
    fq <- trimTailw(fq, 2, opt$demultiplex_qualtrim_code, 5)
    fq <- fq[width(fq) >= opt$demultiplex_qualtrim_minLength]
    
    if(length(fq) > 0) writeFastq(fq, file = file.path(ouputDir, paste0(label, '.', n)), compress = FALSE)
    
    n <- n + 1
  }
}

syncReads <-function(...){
  arguments <- list(...)

  # Create a list of read IDs common to all read arguments.
  n <- Reduce(base::intersect, lapply(arguments, names))

  lapply(arguments, function(x){
    x <- x[names(x) %in% n];
    x[order(names(x))]
  })
}


unpackUniqueSampleID <- function(d){
  d$subject   <- unlist(lapply(strsplit(d$uniqueSample, '~'), '[[', 1))
  d$sample    <- unlist(lapply(strsplit(d$uniqueSample, '~'), '[[', 2))
  d$replicate <- unlist(lapply(strsplit(d$uniqueSample, '~'), '[[', 3))
  d
}


loadSamples <- function(){
  samples <- readr::read_tsv(opt$sampleConfigFile, col_types = readr::cols())
  if(nrow(samples) == 0) stop('Error - no lines of information was read from the sample configuration file.')
  
  if(! 'subject' %in% names(samples))   samples$subject <- 'subject'
  if(! 'replicate' %in% names(samples)) samples$replicate <- 1
  if(any(grepl('~|\\|', paste(samples$subject, samples$sample, samples$replicate)))) stop('Error -- tildas (~) and pipes (|) are reserved characters and can not be used in the subject, sample, or replicate sample configuration columns.')
  
  samples$uniqueSample <- paste0(samples$subject, '~', samples$sample, '~', samples$replicate)
  
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
  if(length(s) == 1 | n_distinct(s) == 1) return(list(0, s[1]))

  if(length(s) > opt$buildFragments_representativeSeqCalc_maxReads){
    set.seed(1)
    s <- sample(s, opt$buildFragments_representativeSeqCalc_maxReads)
  }

  f <- tmpFile()

  # Align sequences in order to handle potential indels.
  inputFile <- file.path(opt$outputDir, 'tmp', paste0(f, '.fasta'))
  s <- Biostrings::DNAStringSet(s)
  names(s) <- paste0('s', 1:length(s))
  Biostrings::writeXStringSet(s, file = inputFile)
  outputFile <- file.path(opt$outputDir, 'tmp', paste0(f, '.representativeSeq.muscle'))

  system(paste(opt$command_muscle, '-quiet -maxiters 1 -diags -in ', inputFile, ' -out ', outputFile))

  if(! file.exists(outputFile)) waitForFile(outputFile)

  s <- as.character(ShortRead::readFasta(outputFile)@sread)
  
  file.remove(outputFile)

  # Create an all vs. all edit distnace matrix.
  m <- as.matrix(stringdist::stringdistmatrix(s, nthread = opt$buildFragments_representativeSeqCalc_CPUs))

  # Create a data frame of string indcies and the sum of their edit distances to all other strings
  # and then order this data frame such that the strings with the small edit distnaces to other
  # strings are at the top of the data frame.
  maxDiffPerNT <- data.frame(n = 1:length(s), diffs = apply(m, 1, sum)) %>% dplyr::arrange(diffs)

  # The function returns both the representative sequence (sequence with the lowest edit distnaces to others)
  # and a metric max edit distance of any string / num characters in the selected representaive.
  #
  #  123456789012345
  #  AGTCAGCTAGCTAGC  max edit distance to other LTRs: 3,  metric 3/15 = 0.2

  # Remove 5% of the most dissimilar reads becasue we do not want an odd-ball read or two to skew the 
  # returned metric. With percentReads = 95, dissimilar reads will not be removed unless more than 20 
  # reads are provided.
  rows <- maxDiffPerNT[1:ceiling(length(s) * (percentReads/100)),]$n
  m <- m[rows, rows]
  s <- s[rows]

  d <- apply(m, 1, sum)
  list(max(apply(m, 1, max) / nchar(s)), gsub('-', '', s[which(d == min(d))[1]]))
}


standardizationSplitVector <- function(d, v){
  if(v == 'replicate'){
    return(d$uniqueSample)
  } else if (v == 'sample'){
    return(d$sample)
  } else {
    return(d$subject)
  }
}


standardizedFragments <- function(frags, opt, cluster){

  source(file.path(opt$softwareDir, 'lib.standardizePositions.R'))
  
  g <- GenomicRanges::makeGRangesFromDataFrame(frags, keep.extra.columns = TRUE)
  g$s <- standardizationSplitVector(g, opt$standardizeFragments_standardizeSitesBy)


  g <- unlist(GenomicRanges::GRangesList(parallel::parLapply(cluster, split(g, g$s), function(x){
  #g <- unlist(GRangesList(lapply(split(g, g$s), function(x){
         source(file.path(opt$softwareDir, 'lib.R'))
         x$intSiteRefined <- FALSE
         out <- tryCatch({
                           o <- standardize_sites(x, counts.col = 'reads')
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

  g$s <- standardizationSplitVector(g, opt$standardizeFragments_standardizeBreaksBy)
  g <- unlist(GenomicRanges::GRangesList(parallel::parLapply(cluster, split(g, g$s), function(x){

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

  parallel::stopCluster(cluster)
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
  
  if(! all(names(x) == a$id)) stop('There was an ordering error during the Golay correction step.')
  
  r <- DNAStringSet(a$seq)
  names(r) <- a$id
  r
}


blast2rearangements <- function(x, minAlignmentLength = 10, minPercentID = 95, CPUs = 20){
  library(data.table)
  
  if(nrow(x) == 0) return(data.frame())
  x <- subset(x, alignmentLength >= minAlignmentLength & gapopen <= 1 & pident >= minPercentID)
  if(nrow(x) == 0) return(data.frame())
  
  cluster <- parallel::makeCluster(CPUs)
  
  # Here we create a spliting variable across the blast data frame
  # being careful not to split read ids into different chunks.
  a <- floor(n_distinct(x$qname) / CPUs)
  b <- 1
  
  z <- dplyr::select(x,  qname, qstart, qend, sstart, send, evalue) %>% group_split(qname)
  z <- as.data.table(bind_rows(mapply(function(x, n){ x$n <- n; x }, z, dplyr::ntile(1:length(z), CPUs), SIMPLIFY = FALSE)))
  
  r <- bind_rows(parallel::parLapply(cluster, split(z, z$n), function(b){
    #r <- rbindlist(lapply(split(z, z$n), function(b){
    library(dplyr)
    library(IRanges)
    library(data.table)
    
    rbindlist(lapply(split(b, b$qname), function(b2){
      # Bin start and end positions to help mitigate sequencing and aligner error.
      # Binned positions are shifted to the lower end of 3 NT intervals.
      breaks <- seq(1, max(b2$qstart), by = 3)
      b2$qstart_binned <- as.integer(as.character(cut(b2$qstart, breaks = c(breaks, Inf), include.lowest = TRUE, labels = breaks)))
      
      # Binned end positions are shifted to the upper end of 3 NT intervals.
      breaks <- seq(min(b2$qend), max(b2$qend), by = 3)
      b2$qend_binned <- as.integer(as.character(cut(b2$qend, breaks = c(breaks, Inf), include.lowest = TRUE, labels = breaks, right = FALSE)))
      
      # Sort BLAST results by query start position and evalue (low to high).
      b2 <- arrange(b2, qstart_binned, evalue)
      b2$strand <- ifelse(b2$send < b2$sstart, '-', '+')
      
      # Alignment to the negative strand will result in the subject end to come before the start.
      # Switch it back so that they are sequential.
      b2$sstart2 <- ifelse(b2$sstart > b2$send, b2$send, b2$sstart)
      b2$send2   <- ifelse(b2$sstart > b2$send, b2$sstart, b2$send)
      
      # Create IRanges
      ir <- IRanges(start = b2$qstart_binned, end = b2$qend_binned)
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
  }))
  
  parallel::stopCluster(cluster)
  r
}


captureLTRseqsLentiHMM <- function(reads, hmm){
  
  # The passed HMM is expected to cover at least 100 NT of the end of the LTR
  # being sequences out of and the HMM is expected to end in CA.
  
  outputFile <- file.path(opt$outputDir, 'tmp', tmpFile())
  writeXStringSet(reads, outputFile)
  comm <- paste0(opt$command.hmmsearch, ' --tblout ', outputFile, '.tbl --domtblout ', outputFile, '.domTbl ', 
                 hmm, ' ', outputFile, ' > ', outputFile, '.hmmsearch')
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
  
  h <- readLines(hmm)
  hmmLength <- as.integer(unlist(strsplit(h[grepl('^LENG', h)], '\\s+'))[2])
  hmmName <- unlist(strsplit(h[grepl('^NAME', h)], '\\s+'))[2]
  
  # Subset HMM results such that alignments start at the start of reads, the end of the HMM
  # which contains the CA is includes and the alignment has a significant alignment scores.
  o <- subset(o, targetStart <= opt$prepReads_HMMmaxStartPos & 
                hmmEnd == hmmLength & 
                fullEval <= as.numeric(opt$prepReads_HMMminEval))
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

