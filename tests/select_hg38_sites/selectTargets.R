library(dplyr)
library(Biostrings)
library(parallel)

set.seed(1)

dbPath    <- '../../data/referenceGenomes/blat/GCA_009914755.4.2bit'
outputDir <- './'

nSites <- 1000

g <- rtracklayer::import.2bit(dbPath)
g <- g[! names(g) %in% c('chrM', 'chrY')]

targets <- bind_rows(lapply(sample(names(g), nSites, replace = TRUE), function(x){
  
  # Retrieve 10 1000 NT wide fragments.
  pos <- sample(10000:(width(g[names(g) == x])-10000), 10)
  frags <- Reduce('append', lapply(pos, function(p) subseq(g[names(g) == x], start = p, width = 10000)))
  names(frags) <- paste('pos', pos)  
  
  # Remove fragments with any Ns and select the first one.
  frags <- frags[! sapply(as.character(frags), function(x) stringr::str_detect(x, 'N'))]
  frag  <- frags[1]
  
  # Define the position as the center of the 10,000 NT frag without Ns.
  tibble(chromosome = x, 
         strand = sample(c('+', '-'), 1), 
         position = as.integer(stringr::str_extract(names(frag), '\\d+')) + 5000,
         desc = 'none')
}))

readr::write_tsv(targets, file.path(outputDir, 'U5_target_candidates.tsv'))   

targets$strand <- ifelse(targets$strand == '+', '-', '+')
targets$position <- targets$position - 4

readr::write_tsv(targets, file.path(outputDir, 'U3_target_candidates.tsv'))
