library(yaml)
library(dplyr)
library(parallel)
library(GenomicRanges)

opt <- read_yaml('config.yml')
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, 'standardizeFragments'))

frags <-readRDS(file.path(opt$outputDir, 'buildFragments', 'fragments.rds'))
frags <- unpackUniqueSampleID(frags)

frags <-standardizedFragments(frags, opt)
frags <- subset(frags, intSiteRefined == TRUE & breakPointRefined == TRUE)

cluster <- makeCluster(opt$standardizeFragments_CPUs)
clusterExport(cluster, c('opt'))

frags <- bind_rows(parLapply(cluster, split(frags, paste0(frags$uniqueSample, frags$seqnames, frags$strand, frags$start, frags$end)), function(x){
           source(file.path(opt$softwareDir, 'lib.R'))
           library(dplyr)
  
            if(nrow(x) > 1){
              r <- representativeSeq(x$repLeaderSeq)
              i <- stringdist::stringdist(r[[2]], x$repLeaderSeq) / nchar(r[[2]]) <= opt$buildFragments_maxLeaderSeqDiffScore
              if(any(! i)) return(data.frame())
    
              x <- x[i,]
              x$reads <- sum(x$reads)
              x <- x[1,]
              x$repLeaderSeq <- r[[2]]
            }
  
            x
        }))

stopCluster(cluster)

