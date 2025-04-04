#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(stringdist))
options(warn=-1)

# Parse command line arguments.
parser <- ArgumentParser()
parser$add_argument("-s", "--sitesFile",    type="character",  default='sites.rds',            help="Path to AAVengeR sites output file (rds).", metavar="")
parser$add_argument("-m", "--multiHitFile", type="character",  default='multiHitClusters.rds', help="Path to AAVengeR multiJitCluster output file (rds).", metavar="")
parser$add_argument("-t", "--truthFile",    type="character",  default='truth.tsv',            help="Path to synthetic data truth file (tsv).", metavar="")
parser$add_argument("-o", "--outputDir",    type="character",  default='out',                  help="Path to output directory.", metavar="")
parser$add_argument("-w", "--siteWidth",    type="integer",    default=5,                      help="Number of NTs to expand truth positions during evaluation.", metavar="")
args <- parser$parse_args()

# Test for input errors.
if(! file.exists(args$sitesFile)){
  message('Error - the sites file could not be found.')
  q(save = "no", status = 1, runLast = FALSE)
}

if(! file.exists(args$multiHitFile)){
  message('Error - the multiHits file could not be found.')
  q(save = "no", status = 1, runLast = FALSE)
}

if(! file.exists(args$truthFile)){
  message('Error - the truth file could not be found.')
  q(save = "no", status = 1, runLast = FALSE)
}

if(! dir.exists(args$outputDir)){
  dir.create(args$outputDir)
}

if(! dir.exists(args$outputDir)){
  message('Error - the output directory could not be created')
  q(save = "no", status = 1, runLast = FALSE)
}

# Read in truth file and AAVengeR output file.
t <- readr::read_tsv(args$truthFile, show_col_types = FALSE)

s <- readRDS(args$sitesFile)
s$posid <- sub('\\.\\d+$', '', s$posid)

m <- readRDS(args$multiHitFile)



posidToGrange <- function(x, expand = 0){
  oo <- tibble(x = x)
  oo <- tidyr::separate(oo, x, c('seqname', 'start'), sep = '[\\+\\-]', remove = FALSE)
  oo$strand <- stringr::str_extract(oo$x, '[\\+\\-]')
  oo$end <- oo$start
  makeGRangesFromDataFrame(oo) + expand
}


tg <- posidToGrange(t$posid, expand = args$siteWidth)
rg <- posidToGrange(s$posid, expand = args$siteWidth)

i <- findOverlaps(tg, rg)

if(any(duplicated(queryHits(i)))){
  
  ti <- unique(queryHits(i)[duplicated(queryHits(i))])
  si <- unique(subset(data.frame(i), queryHits %in% ti)$subjectHits)
  
  message('Identified instances where there was chatter around the expected site.')
  message('Removing ', length(ti), ' truth sites, ', sprintf("%.2f%%", (length(ti) / nrow(t))*100), ' of ', nrow(t), ' sites removed.')
  
  t <- t[-ti,]
  s <- s[-si,]
  
  tg <- posidToGrange(t$posid, expand = args$siteWidth)
  rg <- posidToGrange(s$posid, expand = args$siteWidth)
  
  i <- findOverlaps(tg, rg)
}

recoveredLeaderSeqs <- s[subjectHits(i),]$repLeaderSeq
truthLeaderSeqs     <- t[queryHits(i),]$leaderSeq

d <- mapply(function(a, b){ stringdist(a, b) }, truthLeaderSeqs, recoveredLeaderSeqs, SIMPLIFY = FALSE)

ii <- ! 1:nrow(t) %in% queryHits(i)

matchingMultiHitSites <- 0

if(any(ii)){
  missingSites <- t[ii,]
  multiHitSitesPositions <- unlist(m$posids)
  multiHitSitesPositions <- sub('\\.\\d+$', '' , multiHitSitesPositions)
  
  xg <- posidToGrange(missingSites$posid, expand = args$siteWidth)
  mg <- posidToGrange(multiHitSitesPositions, expand = args$siteWidth)
  
  matchingMultiHitSites <- n_distinct(queryHits(findOverlaps(xg, mg)))
}    

tab1 <- tibble(nSitesExpected = n_distinct(t$posid), 
               nSitesRecovered = n_distinct(s$posid), 
               percentUniqueRecovery = sprintf("%.1f%%", (n_distinct(s$posid) / n_distinct(t$posid))*100),
               percentTotalRecovery = sprintf("%.1f%%", ((n_distinct(s$posid) + matchingMultiHitSites) / n_distinct(t$posid))*100),
               leaderSeqDistMean = sprintf("%.2f", mean(unlist(d))),
               leaderSeqDistSD = sprintf("%.2f", mean(unlist(d))))


tab2 <- bind_rows(lapply(split(t, paste(t$trial, t$subject, t$sample)), function(x){
          r <- subset(s, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1])
          rg <- posidToGrange(r$posid, expand = args$siteWidth)
  
          bind_rows(lapply(split(x, 1:nrow(x)), function(xx){
            tg <- posidToGrange(xx$posid, expand = args$siteWidth)
    
            i <- findOverlaps(tg, rg)
    
            if(length(i) == 1){                                    
              rPos <- as.integer(unlist(strsplit(r[subjectHits(i),]$posid, '[\\+\\-]'))[2])
              tPos <- as.integer(unlist(strsplit(xx$posid, '[\\+\\-]'))[2])
      
              xx$posDiff  <- tPos - rPos
              xx$readDiff <- xx$nReads - r[subjectHits(i),]$reads
              xx$fragDiff <- xx$nFrags - r[subjectHits(i),]$sonicLengths
              xx$leaderSeqDist <- stringdist(xx$leaderSeq, r[subjectHits(i),]$repLeaderSeq)
            } else {
              xx$posDiff  <- NA
              xx$readDiff <- NA
              xx$fragDiff <- NA
              xx$leaderSeqDist <- NA
            }
            
            xx
          }))
        }))

readr::write_tsv(tab1, file.path(args$outputDir, 'table1.tsv'))
readr::write_tsv(tab2, file.path(args$outputDir, 'table2.tsv'))

q(save = 'no', status = 0, runLast = FALSE)
