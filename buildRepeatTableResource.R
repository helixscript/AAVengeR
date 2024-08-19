library(Biostrings)
library(parallel)

# This script calls repeatMasker to create the expected repeat table data file
# needed for the annotateRepeats module. Requires an environment with repeatMasker
# installed and configured.

# Setup.
c <- 20
g <- '/home/everett/rheMac10.std.fasta'
d <- '/home/everett/repeatMasker/rheMac10'
s <- 'macaque'
o <- '/home/everett/rheMac10.repeatTable.gz'

cluster <- makeCluster(c)

comm <- paste0('RepeatMasker -s -species "', s, '"')

dir.create(d)
o <- readDNAStringSet(g)
s <- split(o, names(o))
clusterExport(cluster, c('d', 'comm'))

invisible(parLapply(cluster, s, function(x){
  library(Biostrings)
  f <- paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp')
  writeXStringSet(x, file = f)
  dir.create(paste0(d, '/', names(x)))
  system(paste0(comm, ' -dir ', d, '/', names(x), ' ', f))
  file.remove(f)
}))

stopCluster(cluster)

r <- bind_rows(lapply(system(paste0('find ', d, '/ -name "*.tmp.out"'), intern = TRUE), function(x){
  o <- readLines(x)
  f <- tempfile(tmpdir = '.')
  write(o[4:length(o)], file = f)
  r <- readr::read_table(f, col_names = FALSE, comment = '#')
  r$X1 <- as.integer(r$X1)
  r <- r[! is.na(r$X1),]
  file.remove(f)
  r
}))

names(r) <-  c('SW_score', 'percent_div', 'percent_del', 'percent_ins', 'query_seq', 'query_start', 
               'query_end', 'query_after', 'strand', 'repeat_name', 'repeat_class', 'repeat_start', 
               'repeat_end', 'repeat_after',  'ID',  'alt')

readr::write_tsv(r, o)
