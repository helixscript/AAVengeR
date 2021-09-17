tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }

mixAndChunkSeqs <- function(reads, n){
  reads <- reads[order(width(reads))]
  a <- round(length(reads) / 2)
  b <- unique(c(rbind(c(1:a),  rev(c((a+1):length(reads))))))
  reads <- reads[b]
  split(reads, ceiling(seq_along(reads)/n))
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
    
    fq <- trimTailw(fq, 2, opt$qualtrim_code, 5)
    fq <- fq[width(fq) >= opt$qualtrim_minLength]
    
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




alignReads.BLAT <- function(x, db){
  f <- tmpFile()
  writeFasta(x, file = file.path(opt$outputDir, 'alignedReads', 'blat', paste0(f, '.fasta')))
  
  comm <- paste0(opt$command_blat, ' ', db, ' ',
                 file.path(opt$outputDir, 'alignedReads', 'blat', paste0(f, '.fasta')), ' ',
                 file.path(opt$outputDir, 'alignedReads', 'blat', paste0(f, '.psl')),
                 ' -tileSize=11 -stepSize=9 -minIdentity=90 -out=psl -noHead')
  system(comm)
  file.remove(file.path(opt$outputDir, 'alignedReads', 'blat', paste0(f, '.fasta')))
  
  if(file.exists(file.path(opt$outputDir, 'alignedReads', 'blat', paste0(f, '.psl')))){
    b <- parseBLAToutput(file.path(opt$outputDir, 'alignedReads', 'blat', paste0(f, '.psl')))
    file.remove(file.path(opt$outputDir,  'alignedReads', 'blat', paste0(f, '.psl')))
    b
  } else {
    return(tibble())
  }
}


parseBLAToutput <- function(f){
  if(! file.exists(f) | file.info(f)$size == 0) return(tibble::tibble())
  b <- readr::read_delim(f, delim = '\t', col_names = FALSE, col_types = readr::cols())
  names(b) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand',
                'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')
  
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



standardizedFragments <- function(frags, opt){

  cluster <- parallel::makeCluster(opt$standardizeFragments_CPUs)
  parallel::clusterExport(cluster, 'opt')

  g <- GenomicRanges::makeGRangesFromDataFrame(frags, keep.extra.columns = TRUE)
  g$s <- standardizationSplitVector(g, opt$standardizeFragments_standardizeSitesBy)


  g <- unlist(GenomicRanges::GRangesList(parallel::parLapply(cluster, split(g, g$s), function(x){
  #g <- unlist(GRangesList(lapply(split(g, g$s), function(x){
         source(file.path(opt$softwareDir, 'lib.R'))
         x$intSiteRefined <- FALSE
         out <- tryCatch({
                           o <- gintools::standardize_sites(x)
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
                            o <- gintools::refine_breakpoints(x, counts.col = 'reads')
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
  
  f <- file.path(opt$outputDir, 'demultiplex', 'tmp', paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') )
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




#
# file_exists <- function(f){
#   if(is.na(f)) return(FALSE)
#   file.exists(f)
# }
#
# createReadLengthPlots <- function(anchorReads, adriftReads, file){
#
#   ggplot2::ggsave(ggplot2::ggplot(dplyr::bind_rows(list(tibble(readLength = width(anchorReads), type = 'anchor'),
#                                                         tibble(readLength = width(adriftReads), type = 'adrift'))), aes(readLength)) +
#                     ggplot2::geom_histogram(bins=100) +
#                     ggplot2::facet_grid(.~type), file = file.path(config$outputDir, 'logs', file))
# }
#
#
# findPosIDposIDclusters <- function(x, clusters){
#
#   for(i in 1:length(clusters)){
#       if(any(x %in% clusters[[i]])){
#         return(i)
#       }
#     }
#
#   return(NA)
# }
#
#

#
#
# removeReadFragsWithSameRandomID <- function(frags, config){
#   logMsg(config, 'Removing read pairs which share random linker IDs between samples.', config$logFile)
#
#   randomIDs <- bind_rows(lapply(list.files(file.path(config$outputDir, 'tmp'),
#                                            pattern = 'randomIDs',
#                                            full.names = TRUE), function(x){ load(x); randomIDs }))
#
#   o <- left_join(unnest(frags, readIDs), randomIDs, by = c('readIDs' = 'readID')) %>%
#        #mutate(fragID = paste(subject, sample, strand, start, end)) %>%
#        group_by(randomSeqID) %>%
#        mutate(pass = ifelse(n_distinct(paste(subject, sample)) > 1, FALSE, TRUE)) %>%
#        ungroup() %>%
#        dplyr::select(-fragID)
#
#
#    return(frags)
# }
#
#
#

#
#
#
#
# logMsg <- function(config, msg, logFile = '~/log', append = TRUE){
#   library(lubridate)
#   # Create log time elapsed stamp using the start time in the config object, eg.  [15.5 minutes] log message.
#   te <- paste0('[',
#                sprintf("%.1f",
#                        as.numeric(difftime(ymd_hms(format(Sys.time(), "%y-%m-%d %H:%M:%S")),
#                                            config$startTime,
#                                            units = 'min'))),
#                ' minutes]')
#   write(paste0(te, '\t', msg), file = logFile, append = append)
# }
#
#
# checkConfigFilePaths <- function(files){
#   files <- unique(files)
#   files <- files[! is.na(files)]
#
#   invisible(sapply(unique(files), function(file){
#     if(! file.exists(file)){
#       logMsg(config, paste0('Error. "', file, '" does not exist.'), config$logFile)
#       stop(paste0('Stopping AVVengeR -- "', file, '" not found error.'))
#     }
#   }))
# }
#
#

#
#
#

#
#
#
#

#
#
#
#
# createReadChunks <- function(f, chunkSize, label, ouputDir){
#   # Create a pointer like object to the read file.
#   strm <- FastqStreamer(f, n = as.integer(chunkSize))
#
#   # Extract chunks from file and write out chunks with numeric suffixes, eg. anchorReads.5
#   n <- 1
#   repeat {
#     fq <- yield(strm)
#     if(length(fq) == 0) break
#     writeFastq(fq, file = file.path(ouputDir, paste0(label, '.', n)), compress = FALSE)
#     n <- n + 1
#   }
# }
#


#
#
#

#
# captureLTRseqs <- function(reads, seqs){
#
#   # Group reads by sequence and store read ids in id.list so they can be expanded after analysis.
#   d <- group_by(tibble(ids = names(reads),
#                        read = as.character(reads)), read) %>%
#     summarise(id.list = paste0(ids, collapse=',')) %>%
#     ungroup() %>%
#     mutate(id = paste0('s', 1:n()),
#            qlength = nchar(read))
#
#   f <- tmpFile()
#   write(paste0('>', d$id, '\n', d$read), file = file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))
#
#   seqsNames <- unlist(lapply(unlist(strsplit(seqs, ',')), function(x){
#     x <- unlist(strsplit(x, ':'))
#     x[1]
#   }))
#
#   # Split the expeted ITR/LTR sequences and write them to a single tmp FASTA file.
#   invisible(lapply(unlist(strsplit(seqs, ',')), function(x){
#     x <- unlist(strsplit(x, ':'))
#     write(paste0('>', x[1], '\n', x[2]), file =  file.path(config$outputDir, 'tmp', paste0(f, '.db')), append = TRUE)
#   }))
#
#   # Create a blastn database from the expected ITR/LTR sequences.
#   system(paste0(config$command.makeblastdb, ' -in ', file.path(config$outputDir, 'tmp', paste0(f, '.db')), ' -dbtype nucl -out ', file.path(config$outputDir, 'tmp', f)), ignore.stderr = TRUE)
#   waitForFile(file.path(file.path(config$outputDir, 'tmp', paste0(f, '.nin'))))
#
#
#   # Blast the collapsed anchor reads against the expected ITR/LTR sequences.
#   system(paste0(config$command.blastn, ' -word_size 7 -evalue 50 -outfmt 6 -query ',
#                 file.path(config$outputDir, 'tmp', paste0(f, '.fasta')), ' -db ',
#                 file.path(config$outputDir, 'tmp', f),
#                 ' -num_threads 5 -out ', file.path(config$outputDir, 'tmp', paste0(f, '.blast'))),
#          ignore.stdout = TRUE, ignore.stderr = TRUE)
#
#   waitForFile(file.path(config$outputDir, 'tmp', paste0(f, '.blast')), seconds = 1)
#
#   # Handle case that the returned blast result file is empty.
#   if(file.info(file.path(config$outputDir, 'tmp', paste0(f, '.blast')))$size == 0){
#     invisible(file.remove(list.files(file.path(config$outputDir, 'tmp'), pattern = f, full.names = TRUE)))
#     return(list(reads = DNAStringSet(), LTRs = data.frame(id = NA, LTRname = NA, LTRseq = NA)))
#   }
#
#   b <- read.table(file.path(config$outputDir, 'tmp', paste0(f, '.blast')), sep = '\t', header = FALSE)
#   names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
#   b$alignmentLength <- b$qend - b$qstart + 1
#   b <- left_join(b, d, by = c('qname' = 'id'))
#
#   b <- subset(b, pident  >= config$anchorReads.identification.minPercentSeqID &
#                  qstart  <= config$anchorReads.identification.maxAlignmentStart &
#                  gapopen <= config$anchorReads.identification.maxGapOpen)
#
#   # Handle case where all results were filtered out.
#   if(nrow(b) == 0){
#     invisible(file.remove(list.files(file.path(config$outputDir, 'tmp'), pattern = f, full.names = TRUE)))
#     return(list(reads = DNAStringSet(), LTRs = data.frame(id = NA, LTRname = NA, LTRseq = NA)))
#   }
#
#   # For each read alignment, find the longest alignment and use the ordering of the LTR sequences as a tie breaker.
#   b <- bind_rows(lapply(split(b, b$read), function(x){
#     x <- subset(x, alignmentLength ==  max(x$alignmentLength))
#     x$n <- match(x$sseqid, seqsNames)
#     arrange(x, n) %>% slice(1)
#   }))
#
#  # Expand reads back to single read rows.
#  b$id.list <- strsplit(b$id.list, ',')
#  b <- tidyr::unnest(b, id.list)
#
#  # Extract ITR/LTR sequences from reads.
#  b$LTRseq <- substr(b$read, b$qstart, b$qend)
#
#  # Subset reads to those with significant hits, reorder to match incoming.
#  reads <- reads[names(reads) %in% b$id.list]
#  reads <- reads[match(b$id.list , names(reads))]
#
#  b <- select(b, id.list, sseqid, LTRseq)
#  names(b) <- c('id', 'LTRname', 'LTRseq')
#
#  invisible(file.remove(list.files(file.path(config$outputDir, 'tmp'), pattern = f, full.names = TRUE)))
#
#  return(list(reads = reads, LTRs = b))
# }
#
# captureLTRseqsLentiHMM <- function(reads, hmm){
#
#   # The passed HMM is expected to cover at least 100 NT of the end of the LTR
#   # being sequences out of and the HMM is expected to end in CA.
#
#   outputFile <- file.path(config$outputDir, 'tmp', tmpFile())
#   writeXStringSet(reads, outputFile)
#   comm <- paste0(config$command.hmmsearch, ' --tblout ', outputFile, '.tbl --domtblout ', outputFile, '.domTbl ',
#                 hmm, ' ', outputFile, ' > ', outputFile, '.hmmsearch')
#   system(comm)
#
#   r <- readLines(paste0(outputFile, '.domTbl'))
#   r <- r[!grepl('^\\s*#', r)]
#   r <- strsplit(r, '\\s+')
#
#   o <- bind_rows(lapply(r, function(x) data.frame(t(x))))
#
#   if(nrow(o) == 0) return(list())
#
#   names(o) <- c('targetName', 'targetAcc', 'tlen', 'queryName', 'queryAcc', 'queryLength', 'fullEval',
#                 'fullScore', 'fullBias', 'domNum', 'totalDoms', 'dom_c-Eval', 'dom_i-Eval', 'domScore',
#                 'domBias', 'hmmStart', 'hmmEnd', 'targetStart', 'targetEnd', 'envStart', 'envEnd',
#                 'meanPostProb',  'desc')
#   write.table(o, sep = '\t', file = paste0(outputFile, '.domTbl2'), col.names = TRUE, row.names = FALSE, quote = FALSE)
#
#   o <- readr::read_delim(paste0(outputFile, '.domTbl2'), '\t', col_types = readr::cols())
#
#   h <- readLines(hmm)
#   hmmLength <- as.integer(unlist(strsplit(h[grepl('^LENG', h)], '\\s+'))[2])
#   hmmName <- unlist(strsplit(h[grepl('^NAME', h)], '\\s+'))[2]
#
#   # Subset HMM results such that alignments start at the start of reads, the end of the HMM
#   # which contains the CA is includes and the alignment has a significant alignment scores.
#   o <- subset(o, targetStart <= config$anchorReads.captureLTRseqs.HMMmaxStartPos &
#                  hmmEnd == hmmLength &
#                  fullEval <= as.numeric(config$anchorReads.captureLTRseqs.HMMminEval))
#   if(nrow(o) == 0) return(list())
#
#   reads2 <- reads[names(reads) %in% o$targetName]
#   rm(reads)
#   gc()
#   reads2 <- reads2[match(o$targetName, names(reads2))]
#
#   # Make sure all HMMs alignments result in an CA in the target sequences.
#   reads2 <- reads2[as.character(subseq(reads2, o$targetEnd-1, o$targetEnd)) == 'CA']
#   if(length(reads2) == 0) return(list())
#
#   r <- list()
#   r[['reads']] <- reads2
#   r[['LTRs']]  <- tibble(id = names(reads2),
#                          LTRname = hmmName,
#                          LTRseq = as.character(subseq(reads2, 1, o$targetEnd)))
#   r
# }
#
#
# captureStaticLTRseq <- function(){
#
#   return()
# }
#
#
#
# golayCorrection <- function(x){
#   library(ShortRead)
#   library(dplyr)
#
#   f <- file.path(config$outputDir, 'tmp', paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') )
#   writeFasta(x, file = f)
#
#   system(paste(config$command.python2, file.path(config$softwareDir, 'bin', 'golayCorrection', 'processGolay.py'), f))
#   file.remove(f)
#
#   corrected <- readFasta(paste0(f, '.corrected'))
#   file.remove(paste0(f, '.corrected'))
#
#   a <- left_join(tibble(id = names(x), seq = as.character(x)),
#                  tibble(id = as.character(corrected@id), seq2 = as.character(sread(corrected))), by = 'id')
#   a <- rowwise(a) %>% mutate(editDist = as.integer(adist(seq, seq2))) %>% ungroup()
#
#   i <- which(a$editDist <= 2)
#   a[i,]$seq <- a[i,]$seq2
#
#   if(! all(names(x) == a$id)) stop('There was an ordering error during the Golay correction step.')
#
#   r <- DNAStringSet(a$seq)
#   names(r) <- a$id
#   r
# }
#
#
# collateSampleReads <- function(label){
#   v <- tibble(file = list.files(file.path(config$outputDir, 'tmp'), pattern = paste0('\\.', label, '\\.')),
#               sample = unlist(lapply(strsplit(file, paste0('\\.', label, '\\.')), '[[', 1)))
#   invisible(lapply(split(v, v$sample), function(x){
#     system(paste('cat ', paste0(file.path(config$outputDir, 'tmp', x$file), collapse = ' '), ' > ',
#                  paste0(file.path(config$outputDir, 'sampleReads', paste0(x$sample[1], '.', label, '.fasta')))))
#   }))
# }
#
#
#

#

#
#
#
# trimLeadingSeq <- function(x, seq){
#   f <- tmpFile()
#
#   writeFasta(x,  file = file.path(config$outputDir, 'tmp', f))
#
#   system(paste0(config$command.cutadapt3, ' -f fasta  -e 0.15 -g ', seq, ' --overlap 2 ',
#                 file.path(config$outputDir, 'tmp', f), ' > ',
#                 file.path(config$outputDir, 'tmp', paste0(f, '.out'))),
#          ignore.stderr = TRUE)
#
#   s <- Biostrings::readDNAStringSet(file.path(config$outputDir, 'tmp', paste0(f, '.out')))
#   file.remove(file.path(config$outputDir, 'tmp', f))
#   file.remove(file.path(config$outputDir, 'tmp', paste0(f, '.out')))
#   s
# }
#
#
# # Trim static over read sequences.
# trimOverReadSeq <- function(x, seq, logFile = '~/log'){
#   f <- tmpFile()
#
#   logFile <- paste0(logFile, '.static.', seq)
#   writeFasta(x,  file = file.path(config$outputDir, 'logs', 'cutadapt', f))
#
#   system(paste0(config$command.cutadapt3, ' -f fasta  -e 0.15 -a ', seq, ' --overlap 2 ',
#         ###        '--info-file=', file.path(config$outputDir, 'logs', 'cutadapt', logFile), ' ',
#                 file.path(config$outputDir, 'logs', 'cutadapt', f), ' > ',
#                 file.path(config$outputDir, 'logs', 'cutadapt', paste0(f, '.out'))),
#          ignore.stderr = TRUE)
#
#  ### parseCutadaptLog(file.path(config$outputDir, 'logs', 'cutadapt', logFile))
#
#   x <- Biostrings::readDNAStringSet(file.path(config$outputDir, 'logs', 'cutadapt', paste0(f, '.out')))
#   file.remove(file.path(config$outputDir, 'logs', 'cutadapt', f))
#   file.remove(file.path(config$outputDir, 'logs', 'cutadapt', paste0(f, '.out')))
#   x
# }
#
#
# # Trim dynamic over read sequences.
# trimOverReadSeq2 <- function(reads, o, logFile = '~/log'){
#
#   # Create adapter sequences (RC of last 12 NT of LTRs) for cutadapt to remove.
#   o$LTRs$LTRseqAdapter <- Biostrings::DNAStringSet(substr(o$LTRs$LTRseq, nchar(o$LTRs$LTRseq) - config$anchorReads.identification.minLength, nchar(o$LTRs$LTRseq)))
#   o$LTRs$LTRseqAdapter <- as.character(Biostrings::reverseComplement(o$LTRs$LTRseqAdapter))
#
#   o$LTRs <- o$LTRs[match(names(reads), o$LTRs$id),]
#   if(! all(names(reads) == o$LTRs$id)) stop('ITR/LTR object sort error')
#
#   mcols(reads) <- data.frame(adapter = o$LTRs$LTRseqAdapter)
#
#   Reduce('append', lapply(split(reads, o$LTRs$LTRseqAdapter), function(x){
#          adapter <- as.character(mcols(x)$adapter)[1]
#          f <- tmpFile()
#          Biostrings::writeXStringSet(x, file = file.path(config$outputDir, 'logs', 'cutadapt', f))
#          cutadptLogFile <- file.path(config$outputDir, 'logs', 'cutadapt', paste0(logFile, '.dynamic.', adapter))
#
#          system(paste0(config$command.cutadapt3, ' -f fasta -e 0.30 -a ', adapter, ' --overlap 2 ',
#           ###            '--info-file=', cutadptLogFile, ' ',
#                        file.path(config$outputDir, 'logs', 'cutadapt', f), ' > ',
#                        file.path(config$outputDir, 'logs', 'cutadapt', paste0(f, '.out'))),
#                        ignore.stderr = TRUE)
#
#           ### parseCutadaptLog(cutadptLogFile)
#
#          s <- Biostrings::readDNAStringSet(file.path(config$outputDir, 'logs', 'cutadapt', paste0(f, '.out')))
#          invisible(file.remove(file.path(config$outputDir, 'logs', 'cutadapt', f)))
#          invisible(file.remove(file.path(config$outputDir, 'logs', 'cutadapt', paste0(f, '.out'))))
#          s
#       }))
# }
#

#
#
# getVectorReadIDs <- function(reads, config, vectorFile){
#   d <- group_by(tibble(ids = names(reads),
#                        read = as.character(reads)), read) %>%
#     summarise(id.list = paste0(ids, collapse=',')) %>%
#     ungroup() %>%
#     mutate(id = paste0('s', 1:n()),
#            qlength = nchar(read))
#
#     f <- tmpFile()
#     write(paste0('>', d$id, '\n', d$read), file = file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))
#
#     system(paste0(config$command.makeblastdb, ' -in ', vectorFile, ' -dbtype nucl -out ', file.path(config$outputDir, 'tmp', f)), ignore.stderr = TRUE)
#
#     waitForFile(file.path(file.path(config$outputDir, 'tmp', paste0(f, '.nin'))))
#
#     system(paste0(config$command.blastn, ' -word_size 10 -evalue 10 -outfmt 6 -query ',
#                   file.path(config$outputDir, 'tmp', paste0(f, '.fasta')), ' -db ',
#                   file.path(config$outputDir, 'tmp', f),
#                   ' -num_threads 5 -out ', file.path(config$outputDir, 'tmp', paste0(f, '.blast'))),
#            ignore.stdout = TRUE, ignore.stderr = TRUE)
#
#     waitForFile(file.path(config$outputDir, 'tmp', paste0(f, '.blast')), seconds = 1)
#
#     if(file.info(file.path(config$outputDir, 'tmp', paste0(f, '.blast')))$size == 0) return(character(length = 0))
#
#     b <- read.table(file.path(config$outputDir, 'tmp', paste0(f, '.blast')), sep = '\t', header = FALSE)
#     names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
#     b <- left_join(b, d, by = c('qname' = 'id'))
#     b$pcoverage <- ((b$qend - b$qstart) / b$qlength)*100
#
#     b <- subset(b, pident >= config$alignment.seqFilter.minPercentID & pcoverage >= config$alignment.seqFilter.minPercentQueryCoverage & gapopen <= 1)
#
#     invisible(file.remove(list.files(file.path(config$outputDir, 'tmp'), pattern = f, full.names = TRUE)))
#     unique(unlist(strsplit(b$id.list, ',')))
# }
#
#
# createFragReadAlignments <- function(config, samples, frags){
#   invisible(lapply(split(samples, samples$subject), function(x){
#
#     f <- list.files(file.path(config$outputDir, 'sampleReads'), pattern = paste0(unique(x$sample), collapse = '|'), full.names = TRUE)
#
#     anchorReads <- Reduce('append', lapply(f[grep('anchor', f)], readFasta))
#     adriftReads <- Reduce('append', lapply(f[grep('adrift', f)], readFasta))
#
#     adriftReads <- adriftReads[as.character(adriftReads@id) %in% frags$ids.anchorReads]
#     anchorReads <- anchorReads[as.character(anchorReads@id) %in% frags$ids.anchorReads]
#
#     adriftReads@id <- BStringSet(sub('\\|.+$', '', as.character(adriftReads@id)))
#     anchorReads@id <- BStringSet(sub('\\|.+$', '', as.character(anchorReads@id)))
#
#     if(! all(as.character(adriftReads@id) == as.character(anchorReads@id))) stop('Read id mismatch error.')
#
#     writeFasta(adriftReads, file = file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.adriftReads.fasta')))
#     writeFasta(anchorReads, file = file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.anchorReads.fasta')))
#
#     system(paste0(config$command.bwa, ' mem -M ', x$refGenomeBWAdb[1], ' ',
#                   file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.adriftReads.fasta')), ' ',
#                   file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.anchorReads.fasta')), ' > ',
#                   file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.sam'))))
#
#     system(paste0(config$command.samtools, ' view -S -b ',
#                   file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.sam')), ' > ',
#                   file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.bam'))))
#
#     invisible(file.remove(file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.sam'))))
#
#     system(paste0(config$command.samtools, ' sort -o ',
#                   file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.sorted.bam')), ' ',
#                   file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.bam'))))
#
#     system(paste0(config$command.samtools, ' index ',
#                   file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.sorted.bam'))))
#   }))
# }
#
#
#
# createIntUCSCTrack <- function(d, abundCuts = c(5,10,50),
#                                posColors = c("#8C9DFF", "#6768E3", "#4234C7", "#1D00AB"),
#                                negColors = c("#FF8C8C", "#E35D5D", "#C72E2E", "#AB0000"),
#                                title = 'intSites', outputFile = 'track.ucsc', visibility = 1,
#                                padSite = 0,
#                                siteLabel = NA){
#
#   # Check function inputs.
#   if(length(posColors) != length(negColors))
#     stop('The pos and neg color vectors are not the same length.')
#
#   if(length(abundCuts) != length(posColors) - 1)
#     stop('The number of aundance cut offs must be one less than the number of provided colors.')
#
#   if(! all(c('position', 'end', 'strand', 'seqnames', 'estAbund') %in% names(d)))
#     stop("The expected column names 'start', 'end', 'strand', 'seqnames', 'estAbund' were not found.")
#
#   if(is.na(siteLabel) | ! siteLabel %in% names(d))
#     stop('The siteLabel parameter is not defined or can not be found in your data.')
#
#
#   # Cut the abundance data. Abundance bins will be used to look up color codes.
#   # We flank the provided cut break points with 0 and Inf in order to bin all values outside of breaks.
#   cuts <- cut(d$estAbund, breaks = c(0, abundCuts, Inf), labels = FALSE)
#
#   posColors <- apply(grDevices::col2rgb(posColors), 2, paste0, collapse = ',')
#   negColors <- apply(grDevices::col2rgb(negColors), 2, paste0, collapse = ',')
#
#
#   # Create data fields needed for track table.
#   d$score <- 0
#   d$color <- ifelse(d$strand == '+', posColors[cuts], negColors[cuts])
#
#   # Pad the site n NTs to increase visibility.
#   if(padSite > 0){
#     d$start <- floor(d$start - padSite/2)
#     d$end   <- ceiling(d$end + padSite/2)
#   }
#
#   # Define track header.
#   trackHead <- sprintf("track name='%s' description='%s' itemRgb='On' visibility=%s",
#                        title, title, visibility)
#
#   # Write out track table.
#   write(trackHead, file = outputFile, append = FALSE)
#   write.table(d[, c('seqnames', 'position', 'position', siteLabel, 'score', 'strand', 'start', 'end', 'color')],
#               sep = '\t', col.names = FALSE, row.names = FALSE, file = outputFile, append = TRUE, quote = FALSE)
# }
#
#
# createFragUCSCTrack <- function(d, readCuts = c(5,10,50),
#                                 posColors = c("#8C9DFF", "#6768E3", "#4234C7", "#1D00AB"),
#                                 negColors = c("#FF8C8C", "#E35D5D", "#C72E2E", "#AB0000"),
#                                 title = 'fragments', outputFile = 'track.ucsc', visibility = 1,
#                                 padSite = 0,
#                                 label = NA){
#
#   # Check function inputs.
#   if(length(posColors) != length(negColors))
#     stop('The pos and neg color vectors are not the same length.')
#
#   if(length(readCuts) != length(posColors) - 1)
#     stop('The number of aundance cut offs must be one less than the number of provided colors.')
#
#   if(! all(c('start', 'end', 'strand', 'seqnames', 'reads') %in% names(d)))
#     stop("The expected column names 'start', 'end', 'strand', 'seqnames', 'estAbund' were not found.")
#
#   if(is.na(label) | ! label %in% names(d))
#     stop('The label parameter is not defined or can not be found in your data.')
#
#
#   # Cut the abundance data. Abundance bins will be used to look up color codes.
#   # We flank the provided cut break points with 0 and Inf in order to bin all values outside of breaks.
#   cuts <- cut(d$reads, breaks = c(0, readCuts, Inf), labels = FALSE)
#
#   posColors <- apply(grDevices::col2rgb(posColors), 2, paste0, collapse = ',')
#   negColors <- apply(grDevices::col2rgb(negColors), 2, paste0, collapse = ',')
#
#
#   # Create data fields needed for track table.
#   d$score <- 0
#   d$color <- ifelse(d$strand == '+', posColors[cuts], negColors[cuts])
#
#   # Pad the site n NTs to increase visibility.
#   if(padSite > 0){
#     d$start <- floor(d$start - padSite/2)
#     d$end   <- ceiling(d$end + padSite/2)
#   }
#
#   # Define track header.
#   trackHead <- sprintf("track name='%s' description='%s' itemRgb='On' visibility=%s",
#                        title, title, visibility)
#
#   # Write out track table.
#   write(trackHead, file = outputFile, append = FALSE)
#   write.table(d[, c('seqnames', 'start', 'end', label, 'score', 'strand', 'start', 'end', 'color')],
#               sep = '\t', col.names = FALSE, row.names = FALSE, file = outputFile, append = TRUE, quote = FALSE)
# }
#
#
#
