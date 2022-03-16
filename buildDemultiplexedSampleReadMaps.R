library(dplyr)
library(ggplot2)
library(scales)
library(gridExtra)
library(png)

n <- 50000
w <- 200

configFile <- '/data/project/Encoded/220222_M03249_0243_000000000-JHTY6/config2.yml'

opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))


d <- tibble(file = list.files(file.path(opt$outputDir, 'prepReads', 'final'), pattern = 'anchor', full.names = TRUE))
d$sample <- unlist(lapply(unlist(lapply(d$file, lpe)), function(x){ o <- unlist(strsplit(x, '~'))[2] }))


createPlot <- function(ds, n2){
  dp <- lapply(strsplit(ds, ''), function(x){ tibble(base = x, n = 1:length(x)) })
  dp <- bind_rows(mapply(function(x, n){ x$read <- n; x}, dp, 1:length(dp), SIMPLIFY = FALSE))
  dp$base <- factor(dp$base, levels = c('A', 'T', 'C', 'G', 'N'))
  
  ggplot(dp, aes(n, read, fill = base)) + theme_bw() + geom_tile() +
    scale_fill_manual(values =  c('red', 'green', 'blue', 'gold', 'gray50')) +
    #scale_x_continuous(limits = c(0, w), expand = c(0, 0)) +
    scale_y_continuous(label=comma, limits = c(0, n2), expand = c(0, 0)) +
    labs(x = 'Position', y = 'Reads') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}


r <- Reduce('append', lapply(split(d, d$sample), function(x){
  
       anchorReads <- unique(readDNAStringSet(x$file))
       adriftReads <- unique(readDNAStringSet(sub('anchor', 'adrift', x$file)))
       anchorReads <- anchorReads[width(anchorReads) >= w]
       adriftReads <- adriftReads[width(adriftReads) >= w]
       
       ids <- base::intersect(names(anchorReads), names(adriftReads))
       anchorReads <- anchorReads[names(anchorReads) %in% ids]
       adriftReads <- adriftReads[names(adriftReads) %in% ids]
      
       n2 <- ifelse(length(anchorReads) > n, n, length(adriftReads))
   
       if(n2 < 100){
         o <- DNAStringSet('NNN')
         names(o) <- x$sample[1]
         return(o)
       }
       message(paste0(x$sample[1], ' - ', n2, ' reads'))
       
       set.seed(1)
       i <- sample(1:length(anchorReads), n2)
      
       anchorReads <- anchorReads[i]
       adriftReads <- adriftReads[i]
       
       # Sort reads
       #i <- order(as.character(adriftReads))
       i <- order(as.character(anchorReads))
       
       anchorReads <- anchorReads[i]
       adriftReads <- adriftReads[i]
       
       message(all(names(anchorReads) == names(adriftReads)))
       
       ### if(x$sample == 'tumor1_1') browser()

       p1 <- createPlot(as.character(anchorReads), n2)
       p1 <- p1 + scale_x_continuous(limits = c(0, w), expand = c(0, 0)) +
             ggtitle(paste0('ITR reads: ', x$sample[1], ' - ', n2, ' reads'))
       
       p2 <- createPlot(as.character(adriftReads), n2)
       p2 <- p2 + ggtitle(paste0('Break reads: ', x$sample[1], ' - ', n2, ' reads')) +
             scale_x_reverse(limits = c(w, 0), expand = c(0, 0))
       
       png(filename = paste0(x$sample[1], '.png'))
       grid.arrange(p1, p2, ncol=1)
       dev.off()
       
       #o <- DNAString(gsub('[^A^T^C^G]', 'N', Biostrings::consensusString(anchorReads)))
       o <- DNAStringSet(Biostrings::consensusString(anchorReads))
       names(o) <- x$sample[1]
       o
}))

