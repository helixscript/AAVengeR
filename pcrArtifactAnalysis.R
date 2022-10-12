library(dplyr)
library(lubridate)
library(rtracklayer)
library(Biostrings)
options(stringsAsFactors = FALSE)

vectorFastaFile <- '/home/ubuntu/software/AAVengeR/data/vectors/Encoded_ITR_to_ITR.fasta'
#vectorFastaFile <- '/home/ubuntu/software/AAVengeR/data/vectors/Sabatino_vectors.fasta'

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)

dir.create(file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir))
dir.create(file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs'))
sites <- readRDS(file.path(opt$outputDir, opt$pcrArtifactAnalysis_inputFile))

sites$uniqueSite <- paste0(sites$trial, '~', sites$subject, '~', sites$sample, '~', sites$posid)

system(paste0(opt$command_makeblastdb, ' -in ', vectorFastaFile, ' -dbtype nucl -out ', file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs', 'd')), ignore.stderr = TRUE)

o <- bind_rows(lapply(split(sites, 1:nrow(sites)), function(x){
         a <- unlist(strsplit(x$posid, '[\\+\\-]'))
         a[2] <- sub('\\.\\d+$', '', a[2])
         pos <- as.integer(a[2])
         strand <- stringr::str_extract(x$posid, '[\\+\\-]')
        
         if(strand == '+'){
           seq <- x$upGenome
         } else {
           seq <- as.character(reverseComplement(DNAString(x$downGenome)))
         }

         write(c('>seq', seq), file = file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs', 'seq'), append = FALSE)
         
         system(paste0(opt$command_blastn, ' -word_size 6 -evalue 10 -outfmt 6 -query ',
                       file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs', 'seq'), ' -db ',
                       file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs', 'd'),
                       ' -out ', file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs', 'seq.blast')), ignore.stdout = TRUE, ignore.stderr = TRUE)
         
         if(file.info(file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs', 'seq.blast'))$size > 0){
           b <- read.table(file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs', 'seq.blast'), sep = '\t', header = FALSE)
           if(file.exists(file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs', 'seq'))) file.remove(file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs', 'seq'))
           if(file.exists(file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs', 'seq.blast'))) file.remove(file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs', 'seq.blast'))
           names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
           b$plen <- ((b$qend - b$qstart + 1) / nchar(x$upGenome))*100
           b <- filter(b, pident >= opt$pcrArtifactAnalysis_minPercentSeqID & plen >= opt$pcrArtifactAnalysis_minPercentLength) %>% dplyr::slice_max(bitscore, with_ties = FALSE)
           if(nrow(b) == 1){
             return(tibble(uniqueSite = x$uniqueSite, possiblePCRartifact = TRUE))
           } else {
             return(tibble(uniqueSite = x$uniqueSite, possiblePCRartifact = FALSE))
           }
         } else {
           if(file.exists(file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs', 'seq'))) file.remove(file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs', 'seq'))
           if(file.exists(file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs', 'seq.blast'))) file.remove(file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'dbs', 'seq.blast'))
           return(tibble(uniqueSite = x$uniqueSite, possiblePCRartifact = FALSE))
         }
}))

sites <- select(left_join(sites, o, by = 'uniqueSite'), -uniqueSite) %>% relocate(possiblePCRartifact, .after = posid)

saveRDS(sites, file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, opt$pcrArtifactAnalysis_outputFile))
openxlsx::write.xlsx(sites, file = file.path(opt$outputDir, opt$pcrArtifactAnalysis_outputDir, 'sites.xlsx'))
