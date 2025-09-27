#!/usr/bin/env Rscript

for (p in c('argparse', 'dplyr', 'ShortRead', 'GenomicRanges', 'stringdist', 'stringr')) suppressPackageStartupMessages(library(p, character.only = TRUE))
options(warn=-1)

# Parse command line arguments.
parser <- ArgumentParser()
parser$add_argument("-s", "--sitesFile",    type="character",  default='sites.rds',            help="Path to AAVengeR sites output file (rds).", metavar="")
parser$add_argument("-m", "--multiHitFile", type="character",  default='multiHitClusters.rds', help="Path to AAVengeR multiJitCluster output file (rds).", metavar="")
parser$add_argument("-t", "--truthFile",    type="character",  default='truth.tsv',            help="Path to synthetic data truth file (tsv).", metavar="")
parser$add_argument("-o", "--outputDir",    type="character",  default='out',                  help="Path to output directory.", metavar="")
parser$add_argument("-w", "--siteWidth",    type="integer",    default=5,                      help="Number of NTs to expand truth positions during evaluation.", metavar="")
args <- parser$parse_args()

posidToGrange <- function(x, expand = 0){
  oo <- tibble(x = x)
  oo <- tidyr::separate(oo, x, c('seqname', 'start'), sep = '[\\+\\-]', remove = FALSE)
  oo$strand <- stringr::str_extract(oo$x, '[\\+\\-]')
  oo$end <- oo$start
  makeGRangesFromDataFrame(oo) + expand
}

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

tab <- bind_rows(lapply(split(t, paste(t$trial, t$subject, t$sample)), function(k){
         sampleSites <- subset(s, trial == k$trial[1] & subject == k$subject[1] & sample == k$sample[1])
         sg <- posidToGrange(sampleSites$posid, expand = args$siteWidth)
         sg$posid <- sampleSites$posid
         sg$reads <- sampleSites$reads
         sg$sonicLengths <- sampleSites$sonicLengths
         sg$repLeaderSeq <- sampleSites$repLeaderSeq
        
         sampleMultiHits <- subset(m, trial == k$trial[1] & subject == k$subject[1] & sample == k$sample[1])
         mg <- posidToGrange(sub('\\.\\d+$', '', unlist(sampleMultiHits$posids)), expand = args$siteWidth)
  
         bind_rows(lapply(split(k, 1:nrow(k)), function(x){
           g <- posidToGrange(x$posid, expand = args$siteWidth)
  
           i1 <- findOverlaps(g, sg)
           i2 <- findOverlaps(g, mg)
  
           x$found <- length(i1) > 0 | length(i2) > 0
           x$splitSite <- FALSE
  
          if(length(i1) >= 1){
            
            
            tPos <- as.integer(unlist(strsplit(x$posid, '[\\+\\-]'))[2])
            
            if(length(i1) > 1){
              x$splitSite <- TRUE
              positions <- as.integer(unlist(lapply(str_split(sg[subjectHits(i1)]$posid, '[\\+\\-]'), '[[', 2)))
              absDiff <- abs(positions - tPos)
              i <- which(absDiff == min(absDiff))
              
              if(length(i) == 1){
                message('selecting min dist of overlapping sites.')
                i1 <- i1[i]
              } else {
                reads <- sg[subjectHits(i1)]$reads
                i <- which(reads == max(reads))
                
                if(length(i) == 1){
                  message('selecting max reads of overlapping sites.')
                  i1 <- i1[i]
                } else {
                  message('selecting first site of overlapping sites.')
                  i1 <- i1[1]
                }
              }
            } 
            
            sPos <- as.integer(unlist(str_split(sg[subjectHits(i1)]$posid, '[\\+\\-]'))[2])
      
            x$posDiff  <- tPos - sPos
            x$readDiff <- x$nReads - sg[subjectHits(i1)]$reads
            x$fragDiff <- x$nFrags - sg[subjectHits(i1)]$sonicLengths
            x$leaderSeqDist <- stringdist(x$leaderSeq, sg[subjectHits(i1)]$repLeaderSeq)
          } else {
            x$posDiff  <- NA
            x$readDiff <- NA
            x$fragDiff <- NA
            x$leaderSeqDist <- NA
          }
    
          x
        }))
     }))

sumTab <- tibble(nSitesExpected = n_distinct(t$posid), 
                 nSitesRecovered = n_distinct(s$posid), 
                 percentUniqueRecovery = sprintf("%.2f%%", (n_distinct(tab[! is.na(tab$posDiff),]$posid) / n_distinct(t$posid))*100),
                 percentTotalRecovery = sprintf("%.2f%%",  (n_distinct(tab[tab$found == TRUE,]$posid) / n_distinct(t$posid))*100),
                 leaderSeqDistMean = sprintf("%.2f", mean(tab$leaderSeqDist, na.rm = TRUE)),
                 leaderSeqDistSD = sprintf("%.2f", sd(tab$leaderSeqDist, na.rm = TRUE)))

readr::write_tsv(sumTab, file.path(args$outputDir, 'table1.tsv'))
readr::write_tsv(tab, file.path(args$outputDir, 'table2.tsv'))

q(save = 'no', status = 0, runLast = FALSE)
