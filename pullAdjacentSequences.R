library(dplyr)
library(lubridate)
library(rtracklayer)
library(Biostrings)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$pullAdjacentSequences_outputDir))

sites <- readRDS(file.path(opt$outputDir, opt$pullAdjacentSequences_inputFile))

sites$uniqueSite <- paste0(sites$trial, '~', sites$subject, '~', sites$sample, '~', sites$posid)

o <- bind_rows(lapply(split(sites, sites$refGenome), function(x){

       g <- import(rtracklayer::TwoBitFile(file.path(opt$softwareDir, 'data', 'blatDBs', paste0(x$refGenome[1], '.2bit'))))
       
       r <- bind_rows(lapply(x$uniqueSite, function(x2){
             us <- x2
             o <- unlist(strsplit(x2, '~'))
             x2 <- o[length(o)]
             
             o <- unlist(strsplit(x2, '[\\+\\-]'))
             o[2] <- sub('\\.\\d+$', '', o[2])
             pos <- as.integer(o[2])
             
             i1 <- pos - opt$pullAdjacentSequences_seqLength
             i2 <- pos + opt$pullAdjacentSequences_seqLength
             if(i1 < 1) i1 <- 1 
             s <- g[names(g) == o[1]]
             if(i2 > width(s)) i2 <- width(s)
            
             d <- subseq(s, i1, i1 + (opt$pullAdjacentSequences_seqLength-1))
             u <- subseq(s, pos, pos + opt$pullAdjacentSequences_seqLength-1)
             tibble(uniqueSite = us, strand = stringr::str_extract(x2, '[\\+\\-]'), downGenome = as.character(d), upGenome = as.character(u))   
       }))
       
       left_join(x, r, by = 'uniqueSite')
})) %>% relocate(strand, downGenome, upGenome, .after = posid) %>% arrange(posid) %>% select(-uniqueSite)

saveRDS(o, file = file.path(opt$outputDir, opt$pullAdjacentSequences_outputDir, opt$pullAdjacentSequences_outputFile))
openxlsx::write.xlsx(o, file = file.path(opt$outputDir, opt$pullAdjacentSequences_outputDir, 'sites.xlsx'))

q(save = 'no', status = 0, runLast = FALSE) 