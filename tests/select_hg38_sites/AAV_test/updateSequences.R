library(dplyr)
library(ShortRead)

set.seed(1)
letters <- c('A', 'T', 'C', 'G')

o <- readr::read_tsv('simulatedITRs.tsv')
o$ITRseq2 <- lapply(o$ITRseq, function(x){

  if(! grepl('N', x)){
    return(x)
  } else{
    o <- unlist(lapply(unlist(strsplit(x, '')), function(x2){
      if(x2 == 'N'){
        return(letters[sample(1:4, 1)])
      } else {
        return(x2)
      }
    }))
    
    return(paste0(o, collapse = ''))
  }
})

o$posid <- paste0(o$chromosome, o$strand, o$position)

R2 <- readFastq('../U5_test/U5_syn_R2.fastq.gz')
d <- tibble(id = as.character(R2@id), seq = as.character(R2@sread))
d$posid <- stringr::str_extract(d$id, '^[^_]+')
d$seq <- sub('GAAAATCTCTAGCA', '', d$seq)
d <- left_join(d, select(o, posid, ITRseq2), by = 'posid')
d$seq <- paste0(d$ITRseq2, d$seq)

R2 <- ShortRead::ShortReadQ(sread = DNAStringSet(d$seq), 
                            id = BStringSet(d$id), 
                            quality = BStringSet(sapply(nchar(d$seq), function(x) paste0(rep('?', x), collapse = ''))))

writeFastq(R2, 'AAV_syn_R2.fastq.gz', compress = TRUE)
file.copy('../U5_test/U5_syn_I1.fastq.gz', 'AAV_syn_I1.fastq.gz')
file.copy('../U5_test/U5_syn_R1.fastq.gz', 'AAV_syn_R1.fastq.gz')

