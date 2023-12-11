library(dplyr)
args = commandArgs(trailingOnly=TRUE)

x <- readr::read_tsv(args[1], col_names = FALSE, col_types = 'ccciicccc')
x$gene_id <- stringr::str_extract(x$X9, 'gene_id\\s"([^"]+)')
x$gene_id <- substr(x$gene_id, stringr::str_locate(x$gene_id, '"') + 1, nchar(x$gene_id))

x$gene_ref <- stringr::str_extract(x$X9, 'transcript_id\\s"([^"]+)')
x$gene_ref <- sub('\\.\\d+$', '', substr(x$gene_ref, stringr::str_locate(x$gene_ref, '"') + 1, nchar(x$gene_ref)))

t <- n_distinct(x$gene_id)
n <- 1
r <- bind_rows(lapply(split(x, x$gene_id), function(k){
       if(n %% 100 == 0) message(n, '/', t) 
       n <<- n + 1
       
       tibble(bin = 0, name = k$gene_ref[1], chrom = k$X1[1], strand = k$X7[1], 
              txStart = min(subset(k, X3 == 'transcript')$X4), 
              txEnd = max(subset(k, X3 == 'transcript')$X5),
              cdsStart = min(subset(k, X3 == 'CDS')$X4),  
              cdsEnd = max(subset(k, X3 == 'CDS')$X5), 
              exonCount = nrow(subset(k, X3 == 'exon')), 
              exonsStarts = paste0(subset(k, X3 == 'exon')$X4, collapse = ','),
              exonsEnds = paste0(subset(k, X3 == 'exon')$X5, collapse = ','),
              score = 0, name2 = k$gene_id[1], cdsStartStat = 'unk', cdsEndStat = 'unk', exonFrames = 'unk')
     }))

r$cdsStart <- ifelse(is.infinite(r$cdsStart), NA, r$cdsStart)
r$cdsEnd <- ifelse(is.infinite(r$cdsEnd), NA, r$cdsEnd)

readr::write_tsv(r, args[2], col_names = FALSE)
