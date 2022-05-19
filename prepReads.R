library(ShortRead)
library(dplyr)
library(parallel)
library(lubridate)
options(stringsAsFactors = FALSE)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)

invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

source(file.path(opt$softwareDir, 'lib.R'))

write(c(paste(lubridate::now(), 'Creating required directories.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs'))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'trimmed'))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'unique'))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'dupTables'))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'anchorReadAlignments'))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'final'))

write(c(paste(now(), 'Reading sample data.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
samples <- loadSamples()

if('vectorFastaFile' %in% names(samples)){
  samples$vectorFastaFile <- file.path(opt$softwareDir, 'data', 'vectors', samples$vectorFastaFile)

  if(! all(sapply(unique(samples$vectorFastaFile), file.exists))){
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


cluster <- makeCluster(opt$prepReads_CPUs)
clusterExport(cluster, c('opt', 'samples', 'tmpFile', 'waitForFile'))

# Create a tibble of anchor read files from the demultiplexing output directory.
d <- tibble(file = list.files(file.path(opt$outputDir, opt$prepReads_inputDir), full.names = TRUE, pattern = 'anchor'))

# Extract unique sample ids from file paths.
d$uniqueSample <- unlist(lapply(d$file, function(x){ 
  o <- unlist(strsplit(x, '/'))
  unlist(strsplit(o[length(o)], '\\.'))[1]
})) 

d <- left_join(d, select(samples, uniqueSample, adriftRead.linker.seq, vectorFastaFile), by = 'uniqueSample')

# Extract the common linker -- make its length a setting.
d$adapter <- substr(d$adriftRead.linker.seq, nchar(d$adriftRead.linker.seq) - opt$prepReads_commonLinkerAdapterLength, nchar(d$adriftRead.linker.seq))
d$adapter <- as.character(reverseComplement(DNAStringSet(d$adapter)))

write(c(paste(now(), 'Trimming adapter sequences.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

invisible(parLapply(cluster, split(d, d$file), function(x){
#invisible(lapply(split(d, d$file), function(x){  
  source(file.path(opt$softwareDir, 'lib.R'))
  library(Biostrings)
  
  system(paste0(opt$command_cutadapt, ' -e 0.15 -a ', x$adapter, ' --overlap 2 ',
                x$file, ' > ', file.path(opt$outputDir, opt$prepReads_outputDir, lpe(x$file))), ignore.stderr = TRUE)
  
  # Read in trimmed R1 sequences.
  r1 <- readDNAStringSet(file.path(opt$outputDir, opt$prepReads_outputDir, lpe(x$file)))
  file.remove(file.path(opt$outputDir, opt$prepReads_outputDir, lpe(x$file)))
  r1 <- r1[width(r1) >= opt$prepReads_minAnchorReadLength]
  
  if(length(r1) == 0) return(1)
  
  # Read in corresponding adrift read file and trim linker sequence.
  r2 <- readDNAStringSet(sub('anchorReads', 'adriftReads', x$file))
  r2 <- subseq(r2, nchar(x$adriftRead.linker.seq) + 1)
  r2 <- r2[width(r2) >= opt$prepReads_minAdriftReadLength]
  
  if(length(r2) == 0) return(1)
  
  i <- base::intersect(names(r1), names(r2))
  if(length(i) == 0) return(1)
  
  r1 <- r1[names(r1) %in% i]
  r2 <- r2[names(r2) %in% i]
  
  f <- file.path(opt$outputDir, opt$prepReads_outputDir, 'trimmed', lpe(x$file))
  writeXStringSet(r1, f)
  writeXStringSet(r2, sub('anchorReads', 'adriftReads', f))
  return(0)
}))


# Create unique FASTA files and record duplicate read pairs.

files <- list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'trimmed'), pattern = 'anchorReads', full.names = TRUE)

write(c(paste(now(), 'Creating unique FASTA sequence pairs.')), file = file.path(opt$outputDir, 'log'), append = TRUE)

invisible(parLapply(cluster, files, function(x){
#invisible(lapply(files, function(x){  
  library(Biostrings)
  library(dplyr)
  source(file.path(opt$softwareDir, 'lib.R'))
  message(x)
  
  a <- readDNAStringSet(x)
  b <- readDNAStringSet(sub('anchorReads', 'adriftReads', x))
  
  if(length(a) != length(b)){
    message('Error -- ', x)
    
    return()
  }
  
  d <- tibble(id = names(a), seq = paste0(as.character(a), as.character(b)))
  
  if(any(duplicated(d$seq))){
    o <- d[duplicated(d$seq),]
  
    o1 <- subset(d, seq %in% o$seq)
    o2 <- subset(d, ! seq %in% o$seq)
  
    # Here we create a data frame that stores the ids of duplicate reads that we will be using (id)
    # and the equivalent reads that we will not be using (id2)
    r <- bind_rows(lapply(split(o1, o1$seq), function(x){
           x <- x[order(x$id),]
           tibble(id = x$id[1], n = n_distinct(x$id) - 1, id2 = x$id[2:nrow(x)])
         }))
    
    saveRDS(r, file.path(opt$outputDir, opt$prepReads_outputDir, 'dupTables', paste0(lpe(x), '.rds')))
    
    a <- a[! names(a) %in% r$id2]
    b <- b[! names(b) %in% r$id2]
  }
  
  if(! all(names(a) == names(b))){
    write(c(paste(now(), 'Error -Read names did not match after creating unique FASTA files.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE) 
  }
  
  writeXStringSet(a, file.path(opt$outputDir, opt$prepReads_outputDir, 'unique', lpe(x)))
  writeXStringSet(b, file.path(opt$outputDir, opt$prepReads_outputDir, 'unique', sub('anchorReads', 'adriftReads', lpe(x))))
}))

write(c(paste(now(), 'Capturing duplicated reads')), file = file.path(opt$outputDir, 'log'), append = TRUE)
o <- bind_rows(lapply(list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'dupTables'), full.names = TRUE), readRDS))
saveRDS(o, file.path(opt$outputDir, opt$prepReads_outputDir, 'duplicateReads.rds'))
    
invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

mappings <- tibble()
vectorHits <- tibble()

clusterExport(cluster, c('tmpFile', 'waitForFile', 'opt', 'lpe'))
  
# Align the ends of anchor reads to the vector to identify vector reads which should be removed.

write(c(paste(now(), 'Aligning the ends of anchor reads to vector file.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
invisible(lapply(split(d, d$vectorFastaFile), function(x){
         invisible(file.remove(list.files(file.path(opt$outputDir, opt$vectorFilter_outputDir, 'dbs'), full.names = TRUE)))
  
          system(paste0(opt$command_makeblastdb, ' -in ', x$vectorFastaFile[1], ' -dbtype nucl -out ', file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
          waitForFile(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd.nin'))
  
          files <- sapply(x$file, lpe)
          anchorReadFiles <- file.path(opt$outputDir, opt$prepReads_outputDir, 'unique', files[grepl('anchorRead', files)])
          anchorReadFiles <- anchorReadFiles[file.exists(anchorReadFiles)]
        
          file_n <- 1
          invisible(lapply(anchorReadFiles, function(file){  
            message(file_n, ' / ', length(anchorReadFiles))
            file_n <<- file_n + 1
          
            reads <- readDNAStringSet(file)
    
            # Here we truncate reads to a handful of NTs at their ends to test if they originate from the vector.
            reads <- subseq(reads, (width(reads) - opt$prepReads_vectorAlignmentTestLength) + 1 , width(reads))
            
            b <- bind_rows(parLapply(cluster, mixAndChunkSeqs(reads, opt$prepReads_vectorAlignmentChunkSize), blastReads))
            
            if(nrow(b) > 0){
              readLengths <- tibble(file = lpe(file), qname = names(reads), qlength = width(reads))
              b <- dplyr::left_join(b, readLengths, by = 'qname') 
            
              b$alignmentLength <- b$qend - b$qstart + 1
              b <- dplyr::filter(b, pident >= opt$prepReads_minAlignmentPercentID, alignmentLength >= floor(opt$prepReads_vectorAlignmentTestLength * 0.90), gapopen <= 1)
            }
            
            saveRDS(b, sub('fasta$', 'ends.rds', file.path(opt$outputDir, opt$prepReads_outputDir, 'anchorReadAlignments', lpe(file))))
        }))
}))

# Retrieve the anchor read ends alignments and create a tibble with most significant hit for each read.
write(c(paste(now(), 'Defining vector hits.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
vectorHits <- bind_rows(parLapply(cluster, list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'anchorReadAlignments'), 
                                            pattern = 'ends.rds', full.names = TRUE), function(x){
                  library(dplyr)
                  o <- readRDS(x)
                  
                  if(nrow(o) == 0) return(tibble())
                  o$start <- ifelse(o$sstart > o$send, o$send, o$sstart)
                  o$end <- ifelse(o$sstart > o$send,  o$sstart, o$send)
                  o$strand <- ifelse(o$sstart > o$send, '-', '+')
                  group_by(o, qname) %>% top_n(1, wt = pident) %>% arrange(desc(strand)) %>% dplyr::slice(1) %>% ungroup()
                }))
  
  
saveRDS(vectorHits, file.path(opt$outputDir, opt$prepReads_outputDir, 'vectorHits.rds'))
invisible(file.remove(list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'anchorReadAlignments'), full.names = TRUE)))
  
if(! 'leaderSeqHMM' %in% names(samples)){
    # Now align the full anchor reads to the vector excluding those in vectorHits.
    write(c(paste(now(), 'Aligning full length anchor reads to the vector sequence.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
    invisible(lapply(split(d, d$vectorFastaFile), function(x){
      invisible(file.remove(list.files(file.path(opt$outputDir, opt$vectorFilter_outputDir, 'dbs'), full.names = TRUE)))
    
      system(paste0(opt$command_makeblastdb, ' -in ', x$vectorFastaFile[1], ' -dbtype nucl -out ', file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
      waitForFile(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd.nin'))
    
      files <- sapply(x$file, lpe)
      anchorReadFiles <- file.path(opt$outputDir, opt$prepReads_outputDir, 'unique', files[grepl('anchorRead', files)])
      anchorReadFiles <- anchorReadFiles[file.exists(anchorReadFiles)]
    
      file_n <- 1
      invisible(lapply(anchorReadFiles, function(file){  
        message(file_n, ' / ', length(anchorReadFiles))
        file_n <<- file_n + 1
      
        reads <- readDNAStringSet(file)
      
        # Exclude reads that have ends which align to vector sequences.
        reads <- reads[! names(reads) %in% vectorHits$qname]
        if(length(reads) == 0) return()
      
        b <- bind_rows(parLapply(cluster, mixAndChunkSeqs(reads, opt$prepReads_vectorAlignmentChunkSize), blastReads))
      
        readLengths <- tibble(file = lpe(file), qname = names(reads), qlength = width(reads))
        b <- dplyr::left_join(b, readLengths, by = 'qname') 
        b$alignmentLength <- b$qend - b$qstart + 1
        b <- dplyr::filter(b, pident >= opt$prepReads_minAlignmentPercentID, alignmentLength >= opt$prepReads_minAlignmentLength, gapopen <= 1)
      
        saveRDS(b, sub('fasta$', 'rds', file.path(opt$outputDir, opt$prepReads_outputDir, 'anchorReadAlignments', lpe(file))))
      }))
    }))

    invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))
    stopCluster(cluster)

    f <- list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'anchorReadAlignments'), full.names = TRUE)
    
    write(c(paste(now(), 'Creating read maps from local alignments to the vector file.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
    
    mappings <- bind_rows(lapply(f, function(x){
       b <- readRDS(x)
       
       if(nrow(b) == 0) return(tibble())

       b <- select(b, qname, evalue, alignmentLength, sseqid, qstart, qend, sstart, send)
       
       # Split read alignments into CPU chunks while making sure read ids do not span multiple chunks. 
       z <- split(b, b$qname)
       z <- bind_rows(mapply(function(a, n){ a$n <- n; a}, z, ntile(1:length(z), opt$prepReads_CPUs), SIMPLIFY = FALSE))
       cluster <- makeCluster(opt$prepReads_CPUs)
       clusterExport(cluster, 'opt')
       
       mappings <- bind_rows(parLapply(cluster, split(z, z$n), function(b){
       #mappings <- bind_rows(lapply(split(z, z$n), function(b){   
                        library(GenomicRanges)
                        library(dplyr)

                        # It would be easier to define the alignment region as the min. qStart and max. qEnd
                        # from the blast result. This would allow us to capture NTs not found in the vector file
                        # between NTs which do align but this has the potential to remove too many NTs from 
                        # the start of reads. We can still capture these rearangements since they will not align
                        # to the genome though the extra NTs have a greater potential to shift genomic alignments.
                        # Better to potentially shift alignments than remove genomic NTs.
         
                        bind_rows(lapply(split(b, b$qname), function(a){
                          g <- makeGRangesFromDataFrame(a, ignore.strand = TRUE, seqnames.field = 'sseqid',  start.field = 'qstart', end.field = 'qend')
                          if(length(g) == 0) return(tibble())
                          g <- GenomicRanges::reduce(g, min.gapwidth = 4, ignore.strand = TRUE) # Allow merging if ranges separated by <= 3 NTs.
                          g <- g[start(g) <= opt$prepReads_buildReadMaps_minMapStartPostion]
                          if(length(g) == 0) return(tibble())
                          g <- g[width(g) == max(width(g))][1]
                          return(tibble(id = a$qname[1], leaderMapping.qStart = 1, leaderMapping.qEnd = end(g),
                                        leaderMapping.sStart = NA, leaderMapping.sEnd = NA))
                          
                       }))
                    }))
       
       stopCluster(cluster)
       mappings
    }))
    
    if(nrow(mappings) > 0 & nrow(vectorHits) > 0){
      mappings <- subset(mappings, ! id %in% vectorHits$qname)
    }
    
    saveRDS(mappings, file.path(opt$outputDir, opt$prepReads_outputDir, 'alignments.rds'))
} else {
    
  # If a leaderSeqHMM column is present then we run the anchor reads through the HMM, select reads with 
  # significant HMM hits, and rewrite the mapping object created earlier so that we can use the same code
  # in the alignment module.

  write(c(paste(now(), 'Using leader sequence HMM to define mappings.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
  hmmResults <- bind_rows(lapply(split(d, 1:nrow(d)), function(x){ 
                 source(file.path(opt$softwareDir, 'lib.R'))
                 library(dplyr)
                 library(Biostrings)
                 message(x$file)
                 captureLTRseqsLentiHMM(readDNAStringSet(x$file), subset(samples, uniqueSample == x$uniqueSample)$leaderSeqHMM)
              }))
 
 if(nrow(hmmResults) > 0){
   mappings <- tibble(id = hmmResults$id, leaderMapping.qStart = 1, leaderMapping.qEnd = nchar(hmmResults$LTRseq), 
                      leaderMapping.sStart = NA, leaderMapping.sEnd = NA)
 } else {
   write(c(paste(now(), "Error - no reads matched the HMMs.")), file = file.path(opt$outputDir, 'log'), append = TRUE)
   q(save = 'no', status = 1, runLast = FALSE) 
 }
 
 if(nrow(mappings) > 0 & nrow(vectorHits) > 0){
   mappings <- subset(mappings, ! id %in% vectorHits$qname)
 }
 
 if(nrow(mappings) == 0){
   write(c(paste(now(), "Error - no read mappings remain after vector filter.")), file = file.path(opt$outputDir, 'log'), append = TRUE)
   q(save = 'no', status = 1, runLast = FALSE) 
 }
 
 stopCluster(cluster)  
  
 saveRDS(mappings, file.path(opt$outputDir, opt$prepReads_outputDir, 'alignments.rds'))
}


write(c(paste(now(), "Creating final FASTA files (unique with vector hits removed).")), file = file.path(opt$outputDir, 'log'), append = TRUE)

# Make a list of FASTA files containing unique reads.
files <- list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'unique'), full.names = TRUE)

invisible(lapply(files[grepl('anchorReads', files)], function(file){
  a <- readDNAStringSet(file)
  b <- readDNAStringSet(sub('anchorReads', 'adriftReads', file))
  
  a <- a[! names(a) %in% vectorHits$qname]  
  b <- b[! names(b) %in% vectorHits$qname]
  
  if(length(a) == 0 | length(b) == 0) return()
  
  # Limit reads to those with leader seq HMM hits since we know what
  # the leader sequence should look like because the user provided an HMM.
  if('leaderSeqHMM' %in% names(samples)){
    if(nrow(hmmResults) > 0){
      if(any(names(a) %in% hmmResults$id)) a <- a[names(a) %in% hmmResults$id]
      if(any(names(b) %in% hmmResults$id)) b <- b[names(b) %in% hmmResults$id]
    }
  }
  
  if(length(a) == 0 | length(b) == 0) return()
  
  if(! all(names(a) == names(b))){
    write(c(paste(now(), 'Errror -- read ids do not match before final write out.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
    q(save = 'no', status = 1, runLast = FALSE) 
  }
  
  writeFasta(a, file.path(opt$outputDir, opt$prepReads_outputDir, 'final', lpe(file)))
  writeFasta(b, file.path(opt$outputDir, opt$prepReads_outputDir, 'final', sub('anchorReads', 'adriftReads', lpe(file))))
}))

invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

if(file.exists(file.path(opt$outputDir, opt$demultiplex_outputDir, 'readAttritionTbl.tsv'))){
  
  write(c(paste(now(), "Building read attrition table.")), file = file.path(opt$outputDir, 'log'), append = TRUE)
  
  o <- readr::read_tsv(file.path(opt$outputDir, opt$demultiplex_outputDir, 'readAttritionTbl.tsv'))
  o$preppedReads <- 0
    
  for(x in o$sample){
    file <- file.path(opt$outputDir, opt$prepReads_outputDir, 'final', paste0(x, '.anchorReads.fasta'))
    if(file.exists(file)) o[o$sample == x,]$preppedReads <- length(Biostrings::readDNAStringSet(file))
  }
  
  # Too long
  t <- length(ShortRead::readFastq(opt$demultiplex_anchorReadsFile))
  o <- tibble::add_column(o, .after = 'sample', 'totalReads' = t)
  o$preppedReadsPercentTotal <- (o$preppedReads / o$totalReads)*100
  readr::write_tsv(o, file.path(opt$outputDir, opt$prepReads_outputDir, 'readAttritionTbl.tsv'))
}

if(! opt$prepReads_keepIntermediateFiles){
  unlink(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs'), recursive = TRUE)
  unlink(file.path(opt$outputDir, opt$prepReads_outputDir, 'trimmed'), recursive = TRUE)
  unlink(file.path(opt$outputDir, opt$prepReads_outputDir, 'unique'), recursive = TRUE)
  unlink(file.path(opt$outputDir, opt$prepReads_outputDir, 'dupTables'), recursive = TRUE)
  unlink(file.path(opt$outputDir, opt$prepReads_outputDir, 'anchorReadAlignments'), recursive = TRUE)
}

q(save = 'no', status = 0, runLast = FALSE) 
