library(dplyr)
library(lubridate)
library(rtracklayer)
library(Biostrings)
library(parallel)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir))

sites <- readRDS(file.path(opt$outputDir, opt$predictPCRartifacts_inputFile))


if(! opt$predictPCRartifacts_addAfter %in% names(sites)){
  write(c(paste(now(), paste0('   Error - ', opt$predictPCRartifacts_addAfter, ' is not a column in your input data frame.'))), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}


sites$uniqueSite <- paste0(sites$trial, '~', sites$subject, '~', sites$sample, '~', sites$posid)

o <- bind_rows(lapply(split(sites, sites$refGenome), function(x){

       g <- import(rtracklayer::TwoBitFile(file.path(opt$softwareDir, 'data', 'blatDBs', paste0(x$refGenome[1], '.2bit'))))
       
       r <- bind_rows(lapply(x$uniqueSite, function(x2){
             message(x2)

             us <- x2
             o <- unlist(strsplit(x2, '~'))
             x2 <- o[length(o)]
             
             o <- unlist(strsplit(x2, '[\\+\\-]'))
             strand <- stringr::str_extract(x2, '[\\-\\+]')
             
             o[2] <- sub('\\.\\d+$', '', o[2])
             pos <- as.integer(o[2])
             
             if(! o[1] %in% names(g)) return(tibble(uniqueSite = us, strand = strand, 
                                                    downGenome = paste0(rep('N', opt$predictPCRartifacts_adjacentSeqLength), collapse = ''), 
                                                    upGenome = paste0(rep('N', opt$predictPCRartifacts_adjacentSeqLength), collapse = '')))
             
             if(strand == '-'){
               i1 <- pos - opt$predictPCRartifacts_adjacentSeqLength + 1
               i2 <- pos + opt$predictPCRartifacts_adjacentSeqLength
               if(i1 < 1) i1 <- 1 
               s <- g[names(g) == o[1]]
               if(i2 > width(s)) i2 <- width(s)

               d <- subseq(s, i1, pos)
               u <- subseq(s, pos+1, i2)
               r <-  tibble(uniqueSite = us, strand = strand, downGenome = as.character(d), upGenome = as.character(u))   
             } else {
               i1 <- pos - opt$predictPCRartifacts_adjacentSeqLength
               i2 <- pos + opt$predictPCRartifacts_adjacentSeqLength - 1
               if(i1 < 1) i1 <- 1 
               s <- g[names(g) == o[1]]
               if(i2 > width(s)) i2 <- width(s)
               
               d <- subseq(s, i1, pos-1)
               u <- subseq(s, pos, i2)
               r <-  tibble(uniqueSite = us, strand = strand, downGenome = as.character(d), upGenome = as.character(u))   
             }
             
             r
       }))
       
       left_join(x, r, by = 'uniqueSite')
})) %>% relocate(strand, downGenome, upGenome, .after = posid) %>% arrange(posid) 

dir.create(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs'))

o2 <- bind_rows(lapply(split(o, o$vectorFastaFile), function(a){
  vectorFastaFile <- file.path(opt$softwareDir, 'data', 'vectors', a$vectorFastaFile[1])
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs'), full.names = TRUE)))
  system(paste0(file.path(opt$softwareDir, 'bin', 'makeblastdb'), ' -in ', vectorFastaFile, ' -dbtype nucl -out ', file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)

  #browser()
  
  cluster <- makeCluster(opt$predictPCRartifacts_CPUs)
  clusterExport(cluster, c('opt'))
  
  r <- bind_rows(parLapply(cluster, split(a, 1:nrow(a)), function(x){
  #r <- bind_rows(lapply(split(a, 1:nrow(a)), function(x){
         library(dplyr)
         library(Biostrings)
         source(file.path(opt$softwareDir, 'lib.R'))
    
         a <- unlist(strsplit(x$posid, '[\\+\\-]'))
         a[2] <- sub('\\.\\d+$', '', a[2])
         pos <- as.integer(a[2])
         strand <- stringr::str_extract(x$posid, '[\\+\\-]')
    
         if(strand == '+'){
           seq <- x$upGenome
         } else {
           seq <- as.character(reverseComplement(DNAString(x$downGenome)))
         }
    
         f <- tmpFile()
         
         write(c('>seq', seq), file = file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', f), append = FALSE)
    
         #if(x$posid == 'chr11-112293416') browser()
         
         system(paste0(file.path(opt$softwareDir, 'bin', 'blastn'), ' -dust no -soft_masking false -word_size 6 -evalue 10 -outfmt 6 -query ',
                       file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', f), ' -db ',
                       file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', 'd'),
                       ' -out ', file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', paste0(f, '.blast'))), ignore.stdout = TRUE, ignore.stderr = TRUE)
    
         if(file.info(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', paste0(f, '.blast')))$size > 0){
           b <- read.table(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', paste0(f, '.blast')), sep = '\t', header = FALSE)
           
           
           
           if(file.exists(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', f))) file.remove(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', f))
           if(file.exists(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', paste0(f, '.blast')))) file.remove(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', paste0(f, '.blast')))
     
           names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
           b$len <- b$qend - b$qstart + 1
           b$plen <- (b$qend - b$qstart + 1) / opt$predictPCRartifacts_adjacentSeqLength
           b <- dplyr::filter(b, len >= opt$predictPCRartifacts_minAlignmentLength, pident >= opt$predictPCRartifacts_minPercentSeqID) %>% dplyr::slice_max(bitscore, with_ties = FALSE)

           if(nrow(b) == 1){
            return(tibble(uniqueSite = x$uniqueSite, neighboringSeq = seq, pcrArtifactAlnLength = b$len, pcrArtifactPercentSeqID = b$pident))
           } else {
            return(tibble(uniqueSite = x$uniqueSite, neighboringSeq = seq, pcrArtifactAlnLength = NA, pcrArtifactPercentSeqID = NA))
         }
     } else {
       if(file.exists(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', f))) file.remove(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', f))
       if(file.exists(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', paste0(f, '.blast')))) file.remove(file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, 'dbs', paste0(f, '.blast')))
       return(tibble(uniqueSite = x$uniqueSite, neighboringSeq = seq, pcrArtifactAlnLength = NA, pcrArtifactPercentSeqID = NA))
    }
  }))
  
  stopCluster(cluster)
  r
}))
  
sites <- left_join(sites, o2, by = 'uniqueSite')  %>% relocate(neighboringSeq, pcrArtifactAlnLength, pcrArtifactPercentSeqID, .after = opt$predictPCRartifacts_addAfter)

saveRDS(sites, file = file.path(opt$outputDir, opt$predictPCRartifacts_outputDir, opt$predictPCRartifacts_outputFile))
openxlsx::write.xlsx(sites, file = paste0(sub('\\.\\S+$', '', opt$predictPCRartifacts_outputFile), '.xlsx'))

q(save = 'no', status = 0, runLast = FALSE) 