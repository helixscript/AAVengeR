library(ShortRead)
library(dplyr)
library(parallel)
options(stringsAsFactors = FALSE)

opt <- yaml::read_yaml('config.yml')

invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$prepReads_outputDir))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs'))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'trimmed'))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'unique'))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'dupTables'))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'anchorReadAlignments'))
dir.create(file.path(opt$outputDir, opt$prepReads_outputDir, 'final'))

samples <- loadSamples()

cluster <- makeCluster(opt$prepReads_CPUs)
clusterExport(cluster, c('opt', 'samples'))

d <- tibble(file = list.files(file.path(opt$outputDir, opt$prepReads_inputDir), full.names = TRUE, pattern = 'anchor'))
d$uniqueSample <- unlist(lapply(d$file, function(x){ 
  o <- unlist(strsplit(x, '/'))
  unlist(strsplit(o[length(o)], '\\.'))[1]
})) 

d <- left_join(d, select(samples, uniqueSample, adriftRead.linker.seq, vectorFastaFile), by = 'uniqueSample')

# Trim anchor read overeading and adriftRead linker sequences.

d$adapter <- substr(d$adriftRead.linker.seq, nchar(d$adriftRead.linker.seq) - 9, nchar(d$adriftRead.linker.seq))
d$adapter <- as.character(reverseComplement(DNAStringSet(d$adapter)))

invisible(parLapply(cluster, split(d, d$file), function(x){
#invisible(lapply(split(d, d$file), function(x){  
  source(file.path(opt$softwareDir, 'lib.R'))
  library(Biostrings)
  system(paste0(opt$command_cutadapt, ' -f fasta  -e 0.15 -a ', x$adapter, ' --overlap 2 ',
                x$file, ' > ', file.path(opt$outputDir, opt$prepReads_outputDir, lpe(x$file))), ignore.stderr = TRUE)
  
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
invisible(parLapply(cluster, files, function(x){
  library(Biostrings)
  library(dplyr)
  source(file.path(opt$softwareDir, 'lib.R'))
  
  
  a <- readDNAStringSet(x)
  b <- readDNAStringSet(sub('anchorReads', 'adriftReads', x))
  
  d <- tibble(id = names(a), seq = paste0(as.character(a), as.character(b)))
  
  if(any(duplicated(d$seq))){
    o <- d[duplicated(d$seq),]
  
    o1 <- subset(d, seq %in% o$seq)
    o2 <- subset(d, ! seq %in% o$seq)
  
    # Here we create a data frame that stores the ids of duplicate reads that we will be using (id)
    # and the equivelnt reads that we will not be using (id2)
    r <- bind_rows(lapply(split(o1, o1$seq), function(x){
           x <- x[order(x$id),]
           tibble(id = x$id[1], n = n_distinct(x$id) - 1, id2 = x$id[2:nrow(x)])
         }))
    
    saveRDS(r, file.path(opt$outputDir, opt$prepReads_outputDir, 'dupTables', paste0(lpe(x), '.rds')))
    
    a <- a[! names(a) %in% r$id2]
    b <- b[! names(b) %in% r$id2]
  }
  
  if(! names(a) == names(b)) stop('Read names did not match after creating unique FASTA files.')
  
  writeXStringSet(a, file.path(opt$outputDir, opt$prepReads_outputDir, 'unique', lpe(x)))
  writeXStringSet(b, file.path(opt$outputDir, opt$prepReads_outputDir, 'unique', sub('anchorReads', 'adriftReads', lpe(x))))
}))


o <- bind_rows(lapply(list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'dupTables'), full.names = TRUE), readRDS))
saveRDS(o, file.path(opt$outputDir, opt$prepReads_outputDir, 'duplicateReads.rds'))
    

# Align anchor reads to vector sequences.
# These alignments will be used to determine which parts of the start of anchor reads align to the vector 
# and should be removed during alignments to the genome. This data will be overriden by HMM results if HHM libraries 
# are provided in the sample data. These alignments will also be used to test if the end of anchor reads align to 
# vectors and signal their removal. 

invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

invisible(lapply(split(d, d$vectorFastaFile), function(x){
       invisible(file.remove(list.files(file.path(opt$outputDir, opt$vectorFilter_outputDir, 'dbs'), full.names = TRUE)))
  
        system(paste0(opt$command_makeblastdb, ' -in ', x$vectorFastaFile[1], ' -dbtype nucl -out ', file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
        waitForFile(file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd.nin'))
  
        files <- sapply(x$file, lpe)
        anchorReadFiles <- file.path(opt$outputDir, opt$prepReads_outputDir, 'unique', files[grepl('anchorRead', files)])
        anchorReadFiles <- anchorReadFiles[file.exists(anchorReadFiles)]
        
        file_n <- 1
        bind_rows(lapply(anchorReadFiles, function(file){  
          library(Biostrings)
          library(dplyr)
          source(file.path(opt$softwareDir, 'lib.R'))
    
          message(file_n, ' / ', length(anchorReadFiles))
          file_n <<- file_n + 1
          
          reads <- readDNAStringSet(file)
    
          b <- bind_rows(parLapply(cluster, mixAndChunkSeqs(reads, opt$prepReads_vectorAlignmentChunkSize), function(a){
                 source(file.path(opt$softwareDir, 'lib.R'))
                 f <- tmpFile()
                 writeXStringSet(a,  file.path(opt$outputDir, 'tmp', paste0(f, '.fasta')))
    
                 system(paste0(opt$command_blastn, ' -word_size 6 -evalue 10 -outfmt 6 -query ',
                              file.path(opt$outputDir, 'tmp', paste0(f, '.fasta')), ' -db ',
                              file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd'),
                              ' -num_threads 2 -out ', file.path(opt$outputDir, 'tmp', paste0(f, '.blast'))),
                        ignore.stdout = TRUE, ignore.stderr = TRUE)
    
                 waitForFile(file.path(opt$outputDir, 'tmp', paste0(f, '.blast')))
                 if(file.info(file.path(opt$outputDir, 'tmp', paste0(f, '.blast')))$size == 0) return(tibble())
    
                 b <- read.table(file.path(opt$outputDir, 'tmp', paste0(f, '.blast')), sep = '\t', header = FALSE)
                 names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
                 b
               }))
    
          if(nrow(b) == 0) return(tibble())

          readLengths <- tibble(file = lpe(file), qname = names(reads), qlength = width(reads))
          b <- dplyr::left_join(b, readLengths, by = 'qname') 
          
          b$alignmentLength <- b$qend - b$qstart + 1
          b <- dplyr::filter(b, pident >= opt$prepReads_minAlignmentPercentID, alignmentLength >= opt$prepReads_minAlignmentLength, gapopen <= 1)
          
          saveRDS(b, file.path(opt$outputDir, opt$prepReads_outputDir, 'anchorReadAlignments', lpe(file)))
  }))
}))
  
  
invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))


# Create list of anchor reads with ends that align to the vector and tables of alignments need for trimming reads before genomic alignments.
o <- lapply(list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'anchorReadAlignments'), full.names = TRUE), function(x){
  
       message(date(), ' - processing ', x)
       b <- readRDS(x)
       
       # Alignments already %id filtered. Consider that there may be an accumilation of mismatches near the end.
       # Split reads on query lengths then find alignments which include the end of query sequences. 
       vectorIDs <- unique(unname(unlist(lapply(split(b, b$qlength), function(a){
                        a <- subset(a, qstart <= (a$qlength[1] - opt$prepReads_vectorAlignmentTestLength) & 
                                       qend >= a$qlength[1] - opt$prepReads_vectorAlignmentTestLength_drift)
                        if(nrow(a) > 0){
                          return(a$qname)
                        } else {
                          return(NULL)
                        }
                     }))))
  
       b0 <- b
       b <- subset(b, ! qname %in% vectorIDs)
       message('Removed ', n_distinct(b0$qname) - n_distinct(b$qname), ' reads aligning to the vector')
       
       b <- select(b, qname, evalue, alignmentLength, sseqid, qstart, qend, sstart, send)
 
       # Capturing concatenated lists of hit sequence ids greatly slowed the code below and was omitted.
       mappings <- group_by(b, qname) %>%
         filter(evalue == min(evalue)) %>%
         filter(alignmentLength == max(alignmentLength)) %>%
         summarise(id = qname[1],
                   leaderMapping.qStart = min(qstart), leaderMapping.qEnd = max(qend),
                   leaderMapping.sStart = min(sstart), leaderMapping.sEnd = max(send)) %>%
         ungroup() %>%
         select(-qname)
      
       list(vectorIDs = vectorIDs, mappings = mappings)
})

mappings <- bind_rows(lapply(o, '[[', 'mappings'))
saveRDS(mappings, file.path(opt$outputDir, opt$prepReads_outputDir, 'alignments.rds'))

vectorIDs <- unlist(lapply(o, '[[', 'vectorIDs'))
saveRDS(mappings, file.path(opt$outputDir, opt$prepReads_outputDir, 'vectorIDs.rds'))


# If a leaderSeqHMM column is present then we run the anchor reads through the HMM, select reads with 
# significant HMM hits, and rewrite the mapping object created earlier so that we can use the same code
# in the alignment module.

if('leaderSeqHMM' %in% names(samples)){
 hmmResults <- bind_rows(parLapply(cluster, split(d, 1:nrow(d)), function(x){ 
                 source(file.path(opt$softwareDir, 'lib.R'))
                 library(dplyr)
                 library(Biostrings)
                 captureLTRseqsLentiHMM(readDNAStringSet(x$file), subset(samples, uniqueSample == x$uniqueSample)$leaderSeqHMM)
              }))
 
 
 saveRDS(mappings, file.path(opt$outputDir, opt$prepReads_outputDir, 'blast_alignments.rds'))
 
 mappings <- tibble(id = hmmResults$id, leaderMapping.qStart = 1, leaderMapping.qEnd = nchar(hmmResults$LTRseq), 
             leaderMapping.sStart = NA, leaderMapping.sEnd = NA)
 
 saveRDS(mappings, file.path(opt$outputDir, opt$prepReads_outputDir, 'alignments.rds'))
}


# Make a list of FASTA files containing unique reads.
files <- list.files(file.path(opt$outputDir, opt$prepReads_outputDir, 'unique'), full.names = TRUE)

invisible(lapply(files[grepl('anchorReads', files)], function(file){
  a <- readDNAStringSet(file)
  b <- readDNAStringSet(sub('anchorReads', 'adriftReads', file))
  
  a <- a[! names(a) %in% vectorIDs]
  b <- b[! names(b) %in% vectorIDs]
  
  if(length(a) == 0 | length(b) == 0) return()
  
  # Limit reads to those with leader seq HMM hits since we know what
  # the leader sequence should look like because the user provided an HMM.
  if('leaderSeqHMM' %in% names(samples) & nrow(hmmResults) > 0){
    a <- a[names(a) %in% hmmResults$id]
    b <- b[names(b) %in% hmmResults$id]
  }
  
  if(length(a) == 0 | length(b) == 0) return()
  
  if(! all(names(a) == names(b))) stop('Error -- read ids do not match before final write out.')
  
  writeFasta(a, file.path(opt$outputDir, opt$prepReads_outputDir, 'final', lpe(file)))
  writeFasta(b, file.path(opt$outputDir, opt$prepReads_outputDir, 'final', sub('anchorReads', 'adriftReads', lpe(file))))
}))

invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), full.names = TRUE)))

stopCluster(cluster)



