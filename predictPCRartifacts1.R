library(dplyr)
library(lubridate)
library(rtracklayer)
library(Biostrings)
library(parallel)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$predictPCRartifacts1_outputDir))

sites <- readRDS(file.path(opt$outputDir, opt$predictPCRartifacts1_inputFile))

if(! opt$predictPCRartifacts1_addAfter %in% names(sites)){
  write(c(paste(now(), paste0('   Error - ', opt$predictPCRartifacts1_addAfter, ' is not a column in your input data frame.'))), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}

sites$uniqueSite <- paste0(sites$trial, '~', sites$subject, '~', sites$sample, '~', sites$posid)

cluster <- makeCluster(opt$predictPCRartifacts1_CPUs)
clusterExport(cluster, c('opt', 'tmpFile'))

o <- bind_rows(lapply(split(sites, sites$refGenome), function(x){
  
       x$n <- ntile(1:nrow(x), opt$predictPCRartifacts1_CPUs)
       
       r <- bind_rows(parLapply(cluster, split(x, x$n), function(a){
       #r <- bind_rows(lapply(split(x, x$n), function(a){
              library(dplyr)
              library(Biostrings)
              library(rtracklayer)
         
              if(! grepl('\\.2bit$', a$refGenome[1])) a$refGenome <- paste0(a$refGenome, '.2bit')
              #message(file.path(opt$softwareDir, 'data', 'blatDBs', a$refGenome[1]))
              g <- import(TwoBitFile(file.path(opt$softwareDir, 'data', 'blatDBs', a$refGenome[1])))
         
              bind_rows(lapply(split(a, a$uniqueSite), function(x2){
                us <- x2$uniqueSite
                o <- unlist(strsplit(us, '~'))
                x2 <- o[length(o)]
             
                o <- unlist(strsplit(x2, '[\\+\\-]'))
                strand <- stringr::str_extract(x2, '[\\-\\+]')
             
                o[2] <- sub('\\.\\d+$', '', o[2])
                pos <- as.integer(o[2])
             
                if(! o[1] %in% names(g)) return(tibble(uniqueSite = us, PCRartifact1 = NA))
                s <- g[names(g) == o[1]]

                leaderSeq <- subset(x, posid == x2)$repLeaderSeq
                i <- ifelse(nchar(leaderSeq) < opt$predictPCRartifacts1_adjacentSeqLength, nchar(leaderSeq),opt$predictPCRartifacts1_adjacentSeqLength)
                t <- substr(leaderSeq, nchar(leaderSeq)-i+1, nchar(leaderSeq))
             
                if(strand == '-'){
                  r <- tibble(uniqueSite = us, leaderSeqSeg = t, seq = as.character(reverseComplement( subseq(s, pos+1, pos+opt$predictPCRartifacts1_adjacentSeqLength))))   
                } else {
                  r <- tibble(uniqueSite = us, leaderSeqSeg = t, seq = as.character(subseq(s, pos-opt$predictPCRartifacts1_adjacentSeqLength, pos-1)))   
                }
    
                f <- tmpFile()
                aln <- Biostrings::pairwiseAlignment(r$leaderSeqSeg, r$seq, gapOpening = opt$predictPCRartifacts1_gapOpeningPenalty, gapExtension = opt$predictPCRartifacts1_gapExtensionPenalty)
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
             
                score <- stringr::str_count(m, '\\*')
                o <- paste0(score, ' | ', m)
                if(score < opt$predictPCRartifacts1_minReportMatches) o <- NA
                
                tibble(uniqueSite = us, PCRartifact1 = o)
              }))
            }))
       
       left_join(x, r, by = 'uniqueSite') %>%  relocate(PCRartifact1, .after = opt$predictPCRartifacts1_addAfter) %>% select(-uniqueSite)
}))

saveRDS(o, file = file.path(opt$outputDir, opt$predictPCRartifacts1_outputDir, 'sites.rds'))
openxlsx::write.xlsx(o, file = file.path(opt$outputDir, opt$predictPCRartifacts1_outputDir, 'sites.xlsx'))

q(save = 'no', status = 0, runLast = FALSE) 