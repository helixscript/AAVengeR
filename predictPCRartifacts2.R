library(dplyr)
library(lubridate)
library(rtracklayer)
library(Biostrings)
library(parallel)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir))

sites <- readRDS(file.path(opt$outputDir, opt$predictPCRartifacts2_inputFile))


if(! opt$predictPCRartifacts2_addAfter %in% names(sites)){
  write(c(paste(now(), paste0('   Error - ', opt$predictPCRartifacts2_addAfter, ' is not a column in your input data frame.'))), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}

sites$uniqueSite <- paste0(sites$trial, '~', sites$subject, '~', sites$sample, '~', sites$posid)


o <- bind_rows(lapply(split(sites, sites$refGenome), function(x){

       if(! grepl('\\.2bit$', x$refGenome[1])) x$refGenome <- paste0(x$refGenome, '.2bit')
       
       g <- import(rtracklayer::TwoBitFile(file.path(opt$softwareDir, 'data', 'blatDBs', x$refGenome[1])))
       
       r <- bind_rows(lapply(x$uniqueSite, function(x2){
             us <- x2
             o <- unlist(strsplit(x2, '~'))
             x2 <- o[length(o)]
             
             o <- unlist(strsplit(x2, '[\\+\\-]'))
             strand <- stringr::str_extract(x2, '[\\-\\+]')
             
             o[2] <- sub('\\.\\d+$', '', o[2])
             pos <- as.integer(o[2])
             
            
             if(! o[1] %in% names(g)) return(tibble(uniqueSite = us, strand = strand, seq = paste0(rep('N', opt$predictPCRartifacts2_adjacentSeqLength), collapse = '')))
               
             s <- g[names(g) == o[1]]
             
             if(strand == '-'){
               r <- tibble(uniqueSite = us, strand = strand, seq = as.character(reverseComplement( subseq(s, pos+1-opt$predictPCRartifacts2_adjacentSeqLength, pos) )))   
             } else {
               r <- tibble(uniqueSite = us, strand = strand, seq = as.character(subseq(s, pos, pos+opt$predictPCRartifacts2_adjacentSeqLength-1)))   
             }
             
             r
       }))
       
       left_join(x, r, by = 'uniqueSite')
}))




dir.create(file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs'))

o2 <- bind_rows(lapply(split(o, o$vectorFastaFile), function(a){
  vectorFastaFile <- file.path(opt$softwareDir, 'data', 'vectors', a$vectorFastaFile[1])
  invisible(file.remove(list.files(file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs'), full.names = TRUE)))
  system(paste0('makeblastdb -in ', vectorFastaFile, ' -dbtype nucl -out ', file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)
  
  cluster <- makeCluster(opt$predictPCRartifacts2_CPUs)
  clusterExport(cluster, c('opt'))
  
  r <- bind_rows(parLapply(cluster, split(a, 1:nrow(a)), function(x){
  #r <- bind_rows(lapply(split(a, 1:nrow(a)), function(x){
         library(dplyr)
         library(Biostrings)
         library(stringr)
         source(file.path(opt$softwareDir, 'lib.R'))
    
         f <- tmpFile()
         
         write(c('>seq', x$seq), file = file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs', f), append = FALSE)
         
         system(paste0('blastn -dust no -soft_masking false -evalue 50 -outfmt 6 -word_size ', opt$predictPCRartifacts2_wordSize, 
                       ' -gapopen ', opt$predictPCRartifacts2_gapOpeningPenalty, ' -gapextend ', opt$predictPCRartifacts2_gapExtensionPenalty, ' -query ',
                       file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs', f), ' -db ',
                       file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs', 'd'),
                       ' -out ', file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs', paste0(f, '.blast'))), ignore.stdout = TRUE, ignore.stderr = TRUE)
         
        
         if(file.info(file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs', paste0(f, '.blast')))$size > 0){

           b <- read.table(file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs', paste0(f, '.blast')), sep = '\t', header = FALSE)
           
           names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
           b$matches = b$length - b$mismatch
           
           b <- dplyr::filter(b, matches >= opt$predictPCRartifacts2_minReportMatches, gapopen <= opt$predictPCRartifacts2_maxAlnGaps) %>% dplyr::slice_max(bitscore, with_ties = TRUE)

           if(nrow(b) > 0){
          
             b <- dplyr::slice_max(b, matches, with_ties = FALSE)
             
             f2 <- tmpFile()
             v <- readDNAStringSet(file.path(opt$softwareDir, 'data', 'vectors', x$vectorFastaFile))
             v <- v[names(v) == b$sseqid]
             writeXStringSet(v, file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs', paste0(f2, '.fasta')))
             system(paste0('makeblastdb -in ', file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs', paste0(f2, '.fasta')), ' -dbtype nucl -out ', file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs', f2)), ignore.stderr = TRUE)

             system(paste0('blastn -dust no -soft_masking false -evalue 50 -outfmt 1 -word_size ', opt$predictPCRartifacts2_wordSize, 
                           ' -gapopen ', opt$predictPCRartifacts2_gapOpeningPenalty, ' -gapextend ', opt$predictPCRartifacts2_gapExtensionPenalty, ' -query ',
                           file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs', f), ' -db ',
                           file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs', f2),
                           ' -out ', file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs', paste0(f, '.blast2'))), ignore.stdout = TRUE, ignore.stderr = TRUE)
 
             p <- readLines(file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs', paste0(f, '.blast2')))
             invisible(file.remove(list.files(file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs'), pattern = f2, full.names = TRUE)))
             
             q <- p[grepl('Query_1', p)]
             dashes <- stringr::str_count(str_extract(q, '[ATCGN\\-]+'), '-')
             i <- stringr::str_locate(q, '[ATCGN\\-]+')
             q <- as.integer(unlist(stringr::str_extract_all(sub('Query_1\\s+', '', q), '\\d+')))
             
             p <- p[grepl(paste0('\\.\\s+', b$send, '$'), p)][1]
             p <- substr(p, i[1], i[2])
             p <- gsub('\\.', '*', p)
             p <- gsub('\\s', '.', p)
             p <- unlist(strsplit(p, ''))
             names(p) <- q[1]:(q[2]+dashes)
             
             p <- paste0(unlist(lapply(1:(opt$predictPCRartifacts2_adjacentSeqLength+dashes), function(n){
               n <- as.character(n)
               if(n %in% names(p)) return(p[names(p) == n])
               return('.')
             })), collapse = '')
             
             p <- paste0(b$matches, ' | ', p)
             
             invisible(file.remove(list.files(file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs'), pattern = f, full.names = TRUE)))
          
             return(tibble(uniqueSite = x$uniqueSite, PCRartifact2 = p))
           } else {
             invisible(file.remove(list.files(file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs'), pattern = f, full.names = TRUE)))
             return(tibble(uniqueSite = x$uniqueSite, PCRartifact2 = NA))
         }
     } else {
       invisible(file.remove(list.files(file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'dbs'), pattern = f, full.names = TRUE)))
       return(tibble(uniqueSite = x$uniqueSite, PCRartifact2 = NA))
    }
  }))
  
  stopCluster(cluster)
  r
}))
  
sites <- left_join(sites, o2, by = 'uniqueSite')  %>% relocate(PCRartifact2, .after = opt$predictPCRartifacts2_addAfter)

saveRDS(sites, file = file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'sites.rds'))
openxlsx::write.xlsx(sites, file =file.path(opt$outputDir, opt$predictPCRartifacts2_outputDir, 'sites.xlsx'))

q(save = 'no', status = 0, runLast = FALSE) 