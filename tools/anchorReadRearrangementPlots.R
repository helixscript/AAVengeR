library(tidyverse)

# Read-in / process data replicate and sample level data and create objects (r & s).
# Data should include a factor named label for labeling x-axis. 
# This may require reading in an external label file.
# Seed and label columns needs to be a factor.
#-------------------------------------------------------------------------------
r <- readRDS('output/anchorReadRearrangements2/result_replicates.rds')
r[r$sampleReplicate == 'HelaTAK981~HelaTAK981_set2_TRUF11_TAK_5~GTSP6682~2',]$sampleReplicate <- 'HelaTAK981~HelaTAK981_set2_TRUF11_TAK_5~GTSP6681~2'
r[r$sampleReplicate == 'HelaTAK981~HelaTAK981_set2_TRUF11_TAK_5~GTSP6683~3',]$sampleReplicate <- 'HelaTAK981~HelaTAK981_set2_TRUF11_TAK_5~GTSP6681~3'
r[r$sampleReplicate == 'HelaTAK981~HelaTAK981_set2_TRUF11_TAK_5~GTSP6684~4',]$sampleReplicate <- 'HelaTAK981~HelaTAK981_set2_TRUF11_TAK_5~GTSP6681~4'
r$sampleReplicate <- gsub('GSTP', 'GTSP', r$sampleReplicate)
r$gtsp <- stringr::str_extract(r$sampleReplicate, 'GTSP\\d+')
r$gtsp <- ifelse(grepl('PositiveControl_1', r$sample, ignore.case = TRUE), 'PositiveControl_1', r$gtsp)
r$gtsp <- ifelse(grepl('PositiveControl_2', r$sample, ignore.case = TRUE), 'PositiveControl_2', r$gtsp)
m <- read_tsv('metaData.tsv', show_col_types = FALSE)
r <- left_join(r, m, by = 'gtsp')
r$label <- factor(r$label)
r$seed <- factor(r$seed)


s <- readRDS('output/anchorReadRearrangements2/result_samples.rds')
s <- s[! grepl('TRUF11_TAK_5', s$sample),] # Replicate level corrections can not be applied in a sample level context.
s$sampleReplicate <- gsub('GSTP', 'GTSP', s$sampleReplicate)
s$gtsp <- stringr::str_extract(s$sample, 'GTSP\\d+')
s$gtsp <- ifelse(grepl('PositiveControl_1', s$sample, ignore.case = TRUE), 'PositiveControl_1', s$gtsp)
s$gtsp <- ifelse(grepl('PositiveControl_2', s$sample, ignore.case = TRUE), 'PositiveControl_2', s$gtsp)
m <- read_tsv('metaData.tsv', show_col_types = FALSE)
s <- left_join(s, m, by = 'gtsp')
s$label <- factor(s$label)
s$seed <- factor(s$seed)



# Replicate level plots
#-------------------------------------------------------------------------------

replicates_percentAltReadPlots <- lapply(split(r, paste(r$windowType, r$rarifactionLevel)), function(x){
  s <- sort(x$altStructs)
  m <- s[ceiling(length(s)*0.95):length(s)]
  x[x$altStructs %in% m,]$altStructs <- min(m)
  
  o <- list(ggplot() + theme_bw() +
    scale_fill_gradientn(name = 'altStructs', colors = c("red3", "orange", "yellow", "greenyellow", "darkgreen"),
                         breaks = c(1, floor(max(x$altStructs)/2), max(x$altStructs)),
                         labels = c(1, floor(max(x$altStructs)/2), paste0('>', max(x$altStructs)))) +
    geom_boxplot(data = x, aes(label, percentAltStructReads), color = 'black', outliers = FALSE, alpha = 0) +
    geom_jitter(data = x,  aes(label, percentAltStructReads, fill = altStructs), 
                shape = 21, size = 1.75, width = 0.35, height = 0) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(color = "black"),
          axis.text.x = element_text(angle = -45, vjust = 1.2, hjust=0),
          text = element_text(size = 14)) +
    labs(x = 'Sample', y = 'Percent alt struct reads') + scale_x_discrete(drop=FALSE))
  
  names(o) <- paste0('percentAltStructs_Window_', x$windowType[1], '_Seed_', x$seed[1], '_rarefactionLevel_', x$rarifactionLevel[1]) 
  o
})


replicates_rarifiedAltStructsPerKBplots <- lapply(split(r, paste(r$windowType, r$rarifactionLevel)), function(x){
  s <- sort(x$altStructs)
  m <- s[ceiling(length(s)*0.95):length(s)]
  x[x$altStructs %in% m,]$altStructs <- min(m)
  
  o <- list(ggplot() + theme_bw() +
              scale_shape_manual(values = c(21, 22, 23)) +
              scale_fill_gradientn(name = 'altStructs', colors = c("red3", "orange", "yellow", "greenyellow", "darkgreen"),
                                   breaks = c(1, floor(max(x$altStructs)/2), max(x$altStructs)),
                                   labels = c(1, floor(max(x$altStructs)/2), paste0('>', max(x$altStructs)))) +
              geom_boxplot(data = x, aes(label, altStructsPerKB), color = 'black', outliers = FALSE, alpha = 0) +
              geom_jitter(data = x,  aes(label, altStructsPerKB, fill = altStructs, shape = seed), 
                          size = 1.75, width = 0.35, height = 0) + 
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(color = "black"),
                    axis.text.x = element_text(angle = -45, vjust = 1.2, hjust=0),
                    text = element_text(size = 14)) +
              labs(x = 'Sample', y = paste0('Alt structs / KB')) + scale_x_discrete(drop=FALSE))
  
  names(o) <- paste0('rarifiedAltStructsPerKB_Window_', x$windowType[1], '_Seed_', x$seed[1], '_rarefactionLevel_', x$rarifactionLevel[1]) 
  o
})



# Sample level plots
#-------------------------------------------------------------------------------

samples_percentAltReadPlots <- lapply(split(s, paste(s$windowType, s$rarifactionLevel)), function(x){
  s <- sort(x$altStructs)
  m <- s[ceiling(length(s)*0.95):length(s)]
  x[x$altStructs %in% m,]$altStructs <- min(m)
  
  o <- list(ggplot() + theme_bw() +
              scale_fill_gradientn(name = 'altStructs', colors = c("red3", "orange", "yellow", "greenyellow", "darkgreen"),
                                   breaks = c(1, floor(max(x$altStructs)/2), max(x$altStructs)),
                                   labels = c(1, floor(max(x$altStructs)/2), paste0('>', max(x$altStructs)))) +
              geom_boxplot(data = x, aes(label, percentAltStructReads), color = 'black', outliers = FALSE, alpha = 0) +
              geom_jitter(data = x,  aes(label, percentAltStructReads, fill = altStructs), 
                          shape = 21, size = 1.75, width = 0.35, height = 0) +
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(color = "black"),
                    axis.text.x = element_text(angle = -45, vjust = 1.2, hjust=0),
                    text = element_text(size = 14)) +
              labs(x = 'Sample', y = 'Percent alt struct reads') + scale_x_discrete(drop=FALSE))
  
  names(o) <- paste0('percentAltStructs_Window_', x$windowType[1], '_Seed_', x$seed[1], '_rarefactionLevel_', x$rarifactionLevel[1]) 
  o
})


samples_rarifiedAltStructsPerKBplots <- lapply(split(s, paste(s$windowType, s$rarifactionLevel)), function(x){
  s <- sort(x$altStructs)
  m <- s[ceiling(length(s)*0.95):length(s)]
  x[x$altStructs %in% m,]$altStructs <- min(m)
  
  o <- list(ggplot() + theme_bw() +
              scale_shape_manual(values = c(21, 22, 23)) +
              scale_fill_gradientn(name = 'altStructs', colors = c("red3", "orange", "yellow", "greenyellow", "darkgreen"),
                                   breaks = c(1, floor(max(x$altStructs)/2), max(x$altStructs)),
                                   labels = c(1, floor(max(x$altStructs)/2), paste0('>', max(x$altStructs)))) +
              geom_boxplot(data = x, aes(label, altStructsPerKB), color = 'black', outliers = FALSE, alpha = 0) +
              geom_jitter(data = x,  aes(label, altStructsPerKB, fill = altStructs, shape = seed), 
                          size = 1.75, width = 0.35, height = 0) +
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(color = "black"),
                    axis.text.x = element_text(angle = -45, vjust = 1.2, hjust=0),
                    text = element_text(size = 14)) +
              labs(x = 'Sample', y = paste0('Alt structs / KB')) + scale_x_discrete(drop=FALSE))
  
  names(o) <- paste0('rarifiedAltStructsPerKB_Window_', x$windowType[1], '_Seed_', x$seed[1], '_rarefactionLevel_', x$rarifactionLevel[1]) 
  o
})




