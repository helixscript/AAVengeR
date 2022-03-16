library(dplyr)
library(ggplot2)
library(Biostrings)
library(GenomicRanges)
library(parallel)
source('/home/ubuntu/software/AAVengeR/lib.R')

opt <- yaml::read_yaml('/data/project/Encoded/220222_M03249_0243_000000000-JHTY6/config2.yml')
db <- '/data/project/Encoded/vector_ITR_to_ITR.2bit'
blat <- '/home/ubuntu/software/blat'
d <- tibble(file = list.files('/data/AAVengeR_runs/Encoded2/prepReads/unique', pattern = 'anchor', full.names = TRUE))
d$sample <- unlist(lapply(strsplit(d$file, '~'), '[[', 2))

r <- bind_rows(lapply(split(d, d$sample), function(x){
  
   message(x$sample[1])
   tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }
   f <- tmpFile()
  
   r <- readDNAStringSet(x$file)
   sampleReads <- length(r)
   writeXStringSet(unique(r), f)
   
   system(paste0(blat, ' ', db, ' ', f, ' ', paste0(f, '.psl'), ' -tileSize=10 -stepSize=8 -minIdentity=98 -out=psl -noHead'))
  
   b <- parseBLAToutput(paste0(f, '.psl'))
  
   invisible(file.remove(list.files(pattern = f)))
  
   if(nrow(b) == 0) return(tibble())
  
   b$sample <- x$sample[1]
   b$sampleReads <- sampleReads
   
   r <- dplyr::filter(b, alignmentPercentID >= 98, tNumInsert  <= 1, qNumInsert <= 1, 
                      tBaseInsert <= 2, qBaseInsert <= 2, qStart >= 25) %>%
        dplyr::select(sample, sampleReads, qName, strand, qSize, qStart, qEnd, tName, tSize, 
                      tStart, tEnd, queryPercentID, tAlignmentWidth, queryWidth, alignmentPercentID, percentQueryCoverage) %>%
        group_by(qName) %>%
        arrange(desc(tAlignmentWidth)) %>%
        dplyr::slice(1) %>%
        ungroup()
   r
}))

d <- distinct(tibble(seqnames = 'chr1', subject = r$sample, replicate = 1, reads = 1, 
                     qStart = r$qStart, strand = r$strand, fragStart = r$tStart, fragEnd = r$tEnd))

d <- bind_rows(lapply(split(d, paste(d$subject, d$fragStart, d$fragEnd, d$strand)), function(x){
       x$reads <- nrow(x)
       x[1,]
      }))

cluster <- makeCluster(20)

ds <- standardizedFragments(d, opt, cluster)
ds <- filter(ds, intSiteRefined  == TRUE, breakPointRefined == TRUE) %>% select(subject, strand, start, end) %>% distinct()
ds$width <- ds$end - ds$start + 1
ds <- arrange(ds, subject, start, width)
ds <- bind_rows(lapply(split(ds, ds$subject), function(x){
        x$y <- 1:nrow(x)
        x
       }))

p <- ggplot()+
geom_rect(data=ds[grepl('umor', ds$subject),], mapping=aes(xmin=start, xmax=end, ymin=y, ymax=y+0.5)) +
  theme_bw()+
  labs(x = 'Genome - ITR to ITR', y = 'Unique fragments') +
  facet_wrap(.~subject, scales = 'free_y')+
  theme(plot.margin=grid::unit(c(10,10,10,10), "mm"))

ggsave(p, file = 'p.pdf', width = 10, units = 'in')


