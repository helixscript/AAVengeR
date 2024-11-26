#!/usr/bin/Rscript

# AAVengeR/predictPCRartifacts.R
# John K. Everett, Ph.D.
#
# This script predicts potential PCR artifacts for AAV analyses.
# Two types of predictions are provided.
# 1. alignment between the ends of ITR remnants and gDNA preceding integration sites.
# 2. alignment between the gDNA following integration sites and vector sequences.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(parallel))

set.seed(1)

source(file.path(yaml::read_yaml(commandArgs(trailingOnly=TRUE)[1])$softwareDir, 'lib.R'))
opt <- loadConfig()
optionsSanityCheck()

createOuputDir()
dir.create(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir))
dir.create(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs'))

# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)
write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)

quitOnErorr <- function(msg){
  updateLog(msg)
  message(msg)
  message(paste0('See log for more details: ', opt$defaultLogFile))
  q(save = 'no', status = 1, runLast = FALSE) 
}

sites <- readRDS(file.path(opt$outputDir, opt$predictPCRartifacts_inputFile))

if(! opt$predictPCRartifacts_addAfter %in% names(sites)){
  write(c(paste(now(), paste0('   Error - ', opt$predictPCRartifacts_addAfter, ' is not a column in your input data frame.'))), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}

sites$uniqueSite <- paste0(sites$trial, '~', sites$subject, '~', sites$sample, '~', sites$posid)

cluster <- makeCluster(opt$predictPCRartifacts_CPUs)
clusterExport(cluster, c('opt', 'tmpFile'))

sites$refGenome <- file.path(opt$softwareDir, 'data', 'referenceGenomes', 'blat', paste0(sites$refGenome, '.2bit'))
sites$vectorFastaFile <- file.path(opt$softwareDir, 'data', 'vectors', sites$vector)
sites$vectorFastaFile <- gsub('^\\s*|\\s*$', '', sites$vectorFastaFile)

sites <- bind_rows(lapply(split(sites, sites$refGenome), function(x){
  
       x$n <- ntile(1:nrow(x), opt$predictPCRartifacts_CPUs)
       
       r <- bind_rows(parLapply(cluster, split(x, x$n), function(a){
              library(dplyr)
              library(Biostrings)
              library(rtracklayer)
         
              g <- import(TwoBitFile(a$refGenome[1]))
         
              bind_rows(lapply(split(a, a$uniqueSite), function(x2){
                us <- x2$uniqueSite
                
                o <- unlist(strsplit(us, '~'))
                x2 <- o[length(o)]
             
                o <- unlist(strsplit(x2, '[\\+\\-]'))
                strand <- stringr::str_extract(x2, '[\\-\\+]')
             
                o[2] <- sub('\\.\\S+$', '', o[2])
                pos <- as.integer(o[2])
             
                if(! o[1] %in% names(g)) return(tibble(uniqueSite = us, PCRartifact1 = NA))
                s <- g[names(g) == o[1]]

                leaderSeq <- subset(x, posid == x2)$repLeaderSeq
                i <- ifelse(nchar(leaderSeq) < opt$predictPCRartifacts_adjacentSeqLength, nchar(leaderSeq),opt$predictPCRartifacts_adjacentSeqLength)
                t <- substr(leaderSeq, nchar(leaderSeq)-i+1, nchar(leaderSeq))
                
                if(strand == '-'){
                  endPos <- pos + opt$predictPCRartifacts_adjacentSeqLength
                  startPos <- pos + 1
                  
                  if(startPos > width(s)) startPos <- width(s)
                  if(endPos > width(s))   endPos <- width(s)
                  if(startPos == endPos) startPos <- startPos - 1
                  
                  r <- tibble(uniqueSite = us, leaderSeqSeg = t, seq = as.character(reverseComplement( subseq(s, start = startPos, end = endPos))))   
                } else {
                  startPos <- pos - opt$predictPCRartifacts_adjacentSeqLength
                  endPos <- pos-1
                  
                  if(startPos < 1) startPos <- 1
                  if(endPos < 1)   endPos <- 1
                  if(startPos == endPos) endPos <- endPos + 1
                  
                  r <- tibble(uniqueSite = us, leaderSeqSeg = t, seq = as.character(subseq(s, start = startPos, end = endPos)))   
                }
    
                f <- tmpFile()
                aln <- Biostrings::pairwiseAlignment(r$leaderSeqSeg, r$seq, gapOpening = opt$predictPCRartifacts_gapOpeningPenalty, gapExtension = opt$predictPCRartifacts_gapExtensionPenalty)
                writePairwiseAlignments(aln, f)
                aln <- readLines(f)
                invisible(file.remove(f))
          
                a <- sub('^\\S+\\s+\\d\\s+', '', aln[23])
                a <- sub('\\s+\\d+$', '', a)
             
                b <- sub('^\\S+\\s+\\d\\s+', '', aln[25])
                b <- sub('\\s+\\d+$', '', b)
             
                o <- strsplit(c(a, b), '')
             
                m <- paste0(unlist(lapply(1:nchar(a), function(n){
                       b <- '.'
                       n1 <- o[[1]][n]
                       n2 <- o[[2]][n]
                    
                       if( n1 == n2) b <- '*'
                       if(n1 == '-' | n2 == '-') b <- '-'
                   
                       b
                     })), collapse = '')
             
                # Look at the right half of the alignment and report it if there are >= opt$predictPCRartifacts_minReportHalfMatches matches.
                t <- substr(m, ceiling(opt$predictPCRartifacts_adjacentSeqLength / 2)+1, nchar(m))
                
                if((stringr::str_count(t, '\\*') - stringr::str_count(t, '\\-')) < opt$predictPCRartifacts_minReportHalfMatches){
                  o <- NA
                } else if (stringr::str_count(m, '\\*') < opt$predictPCRartifacts_minReportMatches) {
                  o <- NA
                } else {
                  o <- paste0(stringr::str_count(m, '\\*'), ' | ', m)
                }
                
                tibble(uniqueSite = us, PCRartifact1 = o)
              }))
            }))
       
       left_join(x, r, by = 'uniqueSite') %>% relocate(PCRartifact1, .after = opt$predictPCRartifacts_addAfter)
}))

o <- bind_rows(lapply(split(sites, sites$refGenome), function(x){
  
  g <- import(rtracklayer::TwoBitFile(x$refGenome[1]))
  
  r <- bind_rows(lapply(x$uniqueSite, function(x2){
    us <- x2
    o <- unlist(strsplit(x2, '~'))
    x2 <- o[length(o)]
    
    o <- unlist(strsplit(x2, '[\\+\\-]'))
    strand <- stringr::str_extract(x2, '[\\-\\+]')
    
    o[2] <- sub('\\.\\S+$', '', o[2])
    pos <- as.integer(o[2])
    
    if(! o[1] %in% names(g)) return(tibble(uniqueSite = us, strand = strand, seq = paste0(rep('N', opt$predictPCRartifacts_adjacentSeqLength), collapse = '')))
    
    s <- g[names(g) == o[1]]
    
    if(strand == '-'){
      startPos <- pos+1-opt$predictPCRartifacts_adjacentSeqLength
      endPos <- pos
      
      if(startPos < 1) startPos <- 1
      if(endPos < 1) endPos <- 1
      if(startPos == endPos) endPos <- endPos + 1
      
      r <- tibble(uniqueSite = us, strand = strand, seq = as.character(reverseComplement( subseq(s, start = startPos, end = endPos) )))   
    } else {
      startPos <- pos
      endPos <- pos + opt$predictPCRartifacts_adjacentSeqLength - 1
      
      if(startPos > width(s)) startPos <- width(s)
      if(endPos > width(s))   endPos <- width(s)
      if(startPos == endPos) startPos <- startPos - 1
      
      r <- tibble(uniqueSite = us, strand = strand, seq = as.character(subseq(s, start = startPos, end = endPos)))   
    }
    
    r
  }))
  
  left_join(x, r, by = 'uniqueSite')
}))

o2 <- bind_rows(lapply(split(o, o$vectorFastaFile), function(a){
  vectorFastaFile <- a$vectorFastaFile[1]
    
  #browser()
  
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs'), full.names = TRUE)))
  system(paste0('makeblastdb -in ', vectorFastaFile, ' -dbtype nucl -out ', file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
  
  cluster <- makeCluster(opt$predictPCRartifacts_CPUs)
  clusterExport(cluster, c('opt'))
  
  r <- bind_rows(parLapply(cluster, split(a, 1:nrow(a)), function(x){
  #r <- bind_rows(lapply(split(a, 1:nrow(a)), function(x){
    library(dplyr)
    library(Biostrings)
    library(stringr)
    source(file.path(opt$softwareDir, 'lib.R'))
        
    f <- tmpFile()
    
    write(c('>seq', x$seq), file = file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', f), append = FALSE)
    
    system(paste0('blastn -dust no -soft_masking false -evalue 50 -outfmt 6 -word_size ', opt$predictPCRartifacts_wordSize, 
                  ' -gapopen ', opt$predictPCRartifacts_gapOpeningPenalty, ' -gapextend ', opt$predictPCRartifacts_gapExtensionPenalty, ' -query ',
                  file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', f), ' -db ',
                  file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', 'd'),
                  ' -out ', file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', paste0(f, '.blast'))), ignore.stdout = TRUE, ignore.stderr = TRUE)
        
    if(file.info(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', paste0(f, '.blast')))$size > 0){
      
      b <- read.table(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', paste0(f, '.blast')), sep = '\t', header = FALSE)
      
      names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
      b$matches = b$length - b$mismatch
      
      b <- dplyr::filter(b, matches >= opt$predictPCRartifacts_minReportHalfMatches) %>% dplyr::slice_max(bitscore, with_ties = TRUE)
            
      if(nrow(b) > 0){
                
        b <- dplyr::slice_max(b, matches, with_ties = FALSE)
        
        f2 <- tmpFile()
        v <- readDNAStringSet(x$vectorFastaFile)
        
        b$sseqid <- gsub('^\\s*|\\s*$', '', b$sseqid)
        names(v) <- gsub('^\\s*|\\s*$', '', names(v))
        
        v <- v[names(v) == b$sseqid]
        
        ### if(grepl('mm9', a$refGenome[1])) browser()
        
        writeXStringSet(v, file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', paste0(f2, '.fasta')))
        
        system(paste0('makeblastdb -in ', file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', paste0(f2, '.fasta')), ' -dbtype nucl -out ', file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', f2)), ignore.stderr = TRUE)
                
        system(paste0('blastn -dust no -soft_masking false -evalue 50 -outfmt 1 -word_size ', opt$predictPCRartifacts_wordSize, 
                      ' -gapopen ', opt$predictPCRartifacts_gapOpeningPenalty, ' -gapextend ', opt$predictPCRartifacts_gapExtensionPenalty, ' -query ',
                      file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', f), ' -db ',
                      file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', f2),
                      ' -out ', file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', paste0(f, '.blast2'))), ignore.stdout = TRUE, ignore.stderr = TRUE)
                
        if(file.info(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', paste0(f, '.blast2')))$size > 0){
                    
          p <- readLines(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', paste0(f, '.blast2')))
          invisible(file.remove(list.files(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs'), pattern = f2, full.names = TRUE)))
        
          if(length(p) == 0){
            invisible(file.remove(list.files(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs'), pattern = f, full.names = TRUE)))
            return(tibble(uniqueSite = x$uniqueSite, PCRartifact2 = NA))
          }
                    
          q <- p[grepl('Query_1', p)]
          
          if(length(q) == 0){
            invisible(file.remove(list.files(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs'), pattern = f, full.names = TRUE)))
            return(tibble(uniqueSite = x$uniqueSite, PCRartifact2 = NA))
          }
                  
          dashes <- stringr::str_count(str_extract(q, '[ATCGN\\-]+'), '-')
                    
          i <- stringr::str_locate(q, '[ATCGN\\-]+')
          q <- as.integer(unlist(stringr::str_extract_all(sub('Query_1\\s+', '', q), '\\d+')))
        
          p <- p[grepl(paste0('\\.\\s+', b$send, '$'), p)][1]
          p <- substr(p, i[1], i[2])
          p <- gsub('\\.', '*', p)
          p <- gsub('\\s', '.', p)
          p <- unlist(strsplit(p, ''))
          names(p) <- q[1]:(q[2]+dashes)
        
          p <- paste0(unlist(lapply(1:(opt$predictPCRartifacts_adjacentSeqLength+dashes), function(n){
            n <- as.character(n)
            if(n %in% names(p)) return(p[names(p) == n])
            return('.')
          })), collapse = '')
        
          # Look at the right half of the alignment and report it if there are >= opt$predictPCRartifacts_minReportHalfMatches matches.
          t <- substr(p, ceiling(opt$predictPCRartifacts_adjacentSeqLength / 2)+1, nchar(p))
          
          if((stringr::str_count(t, '\\*') - stringr::str_count(t, '\\-')) < opt$predictPCRartifacts_minReportHalfMatches){
            p <- NA
          } else if (stringr::str_count(p, '\\*') < opt$predictPCRartifacts_minReportMatches) {
            p <- NA
          } else {
            p <- paste0(stringr::str_count(p, '\\*'), ' | ', p)
          }
        
          invisible(file.remove(list.files(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs'), pattern = f, full.names = TRUE)))
          return(tibble(uniqueSite = x$uniqueSite, PCRartifact2 = p))
        } else {
          invisible(file.remove(list.files(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs'), pattern = f, full.names = TRUE)))
          return(tibble(uniqueSite = x$uniqueSite, PCRartifact2 = NA))
        }
      } else {
        invisible(file.remove(list.files(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs'), pattern = f, full.names = TRUE)))
        return(tibble(uniqueSite = x$uniqueSite, PCRartifact2 = NA))
      }
    } else {
      invisible(file.remove(list.files(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs'), pattern = f, full.names = TRUE)))
      return(tibble(uniqueSite = x$uniqueSite, PCRartifact2 = NA))
    }
  }))
  
  stopCluster(cluster)
  r
}))

sites <- left_join(sites, o2, by = 'uniqueSite') %>% relocate(PCRartifact2, .after = PCRartifact1) %>% select(-uniqueSite)
sites <- arrange(sites, desc(sonicLengths))

sites$vectorFastaFile <- sapply(sites$vectorFastaFile, lpe)

sites$refGenome <- sapply(sites$refGenome, lpe)
sites$refGenome <- sub('\\.2bit$', '', sites$refGenome)

if(tolower(opt$databaseConfigGroup) != 'none'){
  suppressPackageStartupMessages(library(RMariaDB))
  uploadSitesToDB(sites)
}

saveRDS(sites, file = file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'sites.rds'))
openxlsx::write.xlsx(sites, file = file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'sites.xlsx'))
readr::write_tsv(sites, file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'sites.tsv.gz'))

q(save = 'no', status = 0, runLast = FALSE) 
