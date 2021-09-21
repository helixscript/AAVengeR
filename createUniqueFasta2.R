library(ShortRead)
library(dplyr)
library(parallel)
options(stringsAsFactors = FALSE)

opt <- yaml::read_yaml('config.yml')

opt$inputDir <- file.path(opt$outputDir, opt$createUniqueFasta_inputDir)

dir.create(file.path(opt$outputDir, opt$createUniqueFasta_outputDir))

d <- data.frame(file = list.files(opt$inputDir, pattern = '*.fasta$'))
d$sample <- unlist(lapply(strsplit(d$file, '\\.'), '[', 1))

cluster <- makeCluster(opt$createUniqueFasta_CPUs)
clusterExport(cluster, 'opt')

date()
invisible(lapply(split(d, d$sample), function(x){
  library(ShortRead)
  library(dplyr)
  
  message(x$sample)
  
  if(nrow(x) != 2) stop('Error - 2 files were not found for sample ', x$sample[1])
  x <- x[order(x$file),]
  
  adriftReads <- readFasta(file.path(opt$inputDir, x[1,]$file))
  anchorReads <- readFasta(file.path(opt$inputDir, x[2,]$file))
  
  message(length(anchorReads), ' / ', length(adriftReads), ' reads')
  
  adriftReads.ids <- as.character(adriftReads@id)
  anchorReads.ids <- as.character(anchorReads@id)
  
  adriftReads <- adriftReads@sread
  anchorReads <- anchorReads@sread
  
  names(adriftReads) <- adriftReads.ids
  names(anchorReads) <- anchorReads.ids
  
  # Consider using concatenated DNA strings and sequence hashes.
  seqs <- paste(as.character(anchorReads), as.character(adriftReads))
  
  if(any(duplicated(seqs))){
    a <- data.frame(id = names(anchorReads), seqs = seqs)
    
    # Consider using the parLapply here.
    b <- bind_rows(parLapply(cluster, split(a, a$seqs), function(s){
           if(nrow(s) > 1){
              return(data.frame(sample = x$sample[1], id = s$id[1], n = (nrow(s) - 1), id_list = I(list(as.character(s$id[2:nrow(s)])))))
           } else {
             return(data.frame())
           }
         }))
  
    u <- unlist(b$id_list)
    anchorReads <- anchorReads[! names(anchorReads) %in% u]
    adriftReads <- adriftReads[! names(adriftReads) %in% u]
  
    b <- bind_rows(lapply(1:nrow(b), function(x){
           x <- b[x,]
           x$id_list <- paste0(unlist(x$id_list), collapse = ',')
           x
         }))
    
    write.table(b, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = file.path(opt$outputDir, opt$createUniqueFasta_outputDir, paste0(x$sample[1], '.dupReadPairMap.tsv')))
  }
  
  writeFasta(adriftReads, file.path(opt$outputDir, opt$createUniqueFasta_outputDir, x[1,]$file))
  writeFasta(anchorReads, file.path(opt$outputDir, opt$createUniqueFasta_outputDir, x[2,]$file))
}))
date()


b <- bind_rows(lapply(list.files(file.path(opt$outputDir, opt$createUniqueFasta_outputDir), pattern = 'dupReadPairMap', full.names = TRUE), function(x){
       read.table(x, sep = '\t', header = TRUE)
     }))

invisible(file.remove(list.files(file.path(opt$outputDir, opt$createUniqueFasta_outputDir), pattern = 'dupReadPairMap', full.names = TRUE)))
write.table(b, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = file.path(opt$outputDir, opt$createUniqueFasta_outputDir, 'dupReadPairMap.tsv'))


                                                                                                                               