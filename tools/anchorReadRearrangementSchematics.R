library(ggplot2)
library(yaml)
library(scales)
library(grid)

nSpecies <- 50
# vector <- 'Sabatino_PMC7855056_singleChain-plasmid.fasta'
vector <- 'Sabatino_PMC7855056_singleChain-plasmid.fasta'
expectedSeqs <- c('1..300[1944+2244]', '1..300[6815-7115]')
configFile <- '~/scratch/AAVengeR_3.0.7/config.yml'
rearrangementOutputDir <- '~/scratch/Sabatino/canine_AAV_factorVIII/inward/240501_M03249_0020_000000000-L4FK6/output/anchorReadRearrangements/'

# DogHemo~pTBG_HC_VCN20~pTBG_HC

opt <- yaml::read_yaml(configFile)

v <- opt$vectors[which(unlist(lapply(opt$vectors, function(x) names(x))) == vector)]

build_color_map <- function(feature_df) {
  vector_length <- unique(feature_df$vectorLength)
  
  color_map <- data.frame(
    pos = 1:vector_length,
    color = "#d3fbf9",
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(feature_df))) {
    row <- feature_df[i, ]
    start <- as.integer(row$start)
    end <- as.integer(row$end)
    vector_length <- as.integer(row$vectorLength)
    
    # Skip invalid ranges
    if (is.na(start) || is.na(end) || start > end || end > vector_length || start < 1) {
      warning(paste("Skipping invalid feature:", row$featureName, "start=", start, "end=", end))
      next
    }
    
    feature_positions <- start:end
    n_pos <- length(feature_positions)
    
    gradient_colors <- gradient_n_pal(c(row$startColor, row$endColor))(seq(0, 1, length.out = n_pos))
    color_map$color[feature_positions] <- gradient_colors
  }
  
  return(setNames(color_map$color, color_map$pos))
}

feature_df <- bind_rows(lapply(v, function(x){
  bind_rows(lapply(x[[1]]$features, function(xx){
    feature <- names(xx)
    tibble(vector = names(x), featureName = names(xx), start = xx[[1]]$start,  end = xx[[1]]$end,
           startColor = xx[[1]]$startColor, endColor = xx[[1]]$endColor, vectorLength = x[[1]]$length)
  }))
}))

color_map <- build_color_map(feature_df)

rearrangement_table <- readRDS(file.path(rearrangementOutputDir, 'result_samples_altStructs.rds')) %>%
  select(sample, seed, vector, rarifactionLevel, count, seqModel) %>%
  filter(vector == vector, rarifactionLevel == 0, seed == 1) %>%
  mutate(rearrangements = stringr::str_count(seqModel, ';')) %>%
  arrange(desc(count), desc(rearrangements)) %>%
  select(sample, count, seqModel) %>%
  dplyr::rename(shorthand = seqModel) %>%
  filter(grepl('^[12]\\.\\.', shorthand))

stats_table <- readRDS(file.path(rearrangementOutputDir, 'result_samples.rds')) %>%
  filter(vector == vector, rarifactionLevel == 0, seed == 1) %>%
  select(sample, totalReads, expectedStructReads)

rearrangement_table <- rearrangement_table[grepl('M50|M66', rearrangement_table$sample),]
stats_table <- stats_table[grepl('M50|M66', stats_table$sample),]

rearrangement_table <- left_join(rearrangement_table, stats_table, by = 'sample')
rearrangement_table$percent <- rearrangement_table$count / rearrangement_table$totalReads
rearrangement_table <- select(rearrangement_table, -totalReads, -expectedStructReads)


fullVector_table <- tibble(sample = 'vector', count = 1, 
                           shorthand = paste0('1..', v[[1]][[1]]$length, '[1+', v[[1]][[1]]$length, ']'),
                           totalReads = 1, percent = 1)

expectedSeq_table <- bind_rows(lapply(split(stats_table, stats_table$sample), function(x){
  tibble(sample = x$sample, count = NA, 
         shorthand = expectedSeqs, percent = x$expectedStructReads / x$totalReads)
}))


plots <- lapply(split(rearrangement_table, rearrangement_table$sample), function(tbl){
  
  y <- 1
  
  tbl <- arrange(tbl, desc(percent), desc(count))
  tbl <- tbl[1:nSpecies,]
  tbl <- bind_rows(subset(expectedSeq_table, sample == tbl$sample[1]), tbl)
  tbl <- arrange(tbl, desc(row_number()))
  
  o <- lapply(split(tbl, 1:nrow(tbl)), function(x){
    o <- unlist(strsplit(x$shorthand, ';'))
    a <- stringr::str_extract(o, '\\d+\\.\\.\\d+')
    a.start <-  stringr::str_extract(a, '^\\d+')
    a.end <-  stringr::str_extract(a, '\\d+$')
    
    b <- sub('\\[', '', stringr::str_extract(o, '\\[[^\\]]+'))
    b.start <-  vector()
    b.end <- vector()
    
    invisible(lapply(strsplit(b, '[\\+\\-]'), function(xx){
      if(! grepl('x', xx[1], ignore.case = TRUE)){
        b.start <<- append(b.start, xx[1])
        b.end <<- append(b.end, xx[2])
      } else {
        b.start <<- append(b.start, NA)
        b.end <<- append(b.end, NA)
      }
    }))
    
    k <- tibble(x1 = as.integer(a.start), x2 = as.integer(a.end), colorStart = as.integer(b.start))
    
    if(nrow(k) > 1){
      for(i in 2:nrow(k)){
        k[i,]$x1 <- k[i-1,]$x2 + 1
      }
    }
    
    k$strand <- stringr::str_extract(o, '[\\+\\-]')
    k$width <- k$x2 - k$x1 + 1
    
    m <- bind_rows(lapply(split(k, 1:nrow(k)), function(xx){
      if(! is.na(xx$strand)){
        colors <- unname(color_map[xx$colorStart:(xx$colorStart + xx$width - 1)])
      } else {
        colors <- 'white'
      }
      
      if(is.na(xx$strand)) xx$strand <- '+'
      if(xx$strand == '-') colors <- rev(colors)
      
      xx <- tidyr::uncount(xx, width)
      xx$x1 <- seq(from = xx$x1[1], to = xx$x2[1], by = 1)
      xx$x2 <- xx$x1 + 1
      xx$colors <- colors
      xx
    }))
    
    p <- tibble(p = sprintf("%.2f%%", x$percent * 100))
    
    p$y <- k$y1 <- m$y1 <- y
    k$y2 <- m$y2 <- y + 1
    y <<- y + 1
    
    list(k, m, p)
  })
  
  
  boxes  <- bind_rows(lapply(o, '[[', 1))
  cells  <- bind_rows(lapply(o, '[[', 2))
  labels <- bind_rows(lapply(o, '[[', 3))
  

  #browser()
  
  # Remove one of the top labels since they are duplicates 
  labels <- labels[1:(nrow(labels) - 1),]  
  
  # Increase y for the remaining top label.
  labels[nrow(labels),]$y <- labels[nrow(labels),]$y + 1.5
  
  # Raise the top two schematics.
  cells[cells$y1 %in% rev(sort(unique(cells$y1)))[1:2],]$y1 <- cells[cells$y1 %in% rev(sort(unique(cells$y1)))[1:2],]$y1 + 1
  cells[cells$y2 %in% rev(sort(unique(cells$y2)))[1:2],]$y1 <- cells[cells$y2 %in% rev(sort(unique(cells$y2)))[1:2],]$y2 + 1
  
  boxes[boxes$y1 %in% rev(sort(unique(boxes$y1)))[1:2],]$y1 <- boxes[boxes$y1 %in% rev(sort(unique(boxes$y1)))[1:2],]$y1 + 1
  boxes[boxes$y2 %in% rev(sort(unique(boxes$y2)))[1:2],]$y1 <- boxes[boxes$y2 %in% rev(sort(unique(boxes$y2)))[1:2],]$y2 + 1
  
 
  
  ggplot() +
    theme_void() +
    geom_rect(data = cells, mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = colors)) +
    geom_rect(data = boxes, mapping = aes(xmin = x1, xmax = x2 + 1, ymin = y1, ymax = y2), color = 'black', alpha = 0) +
    geom_text(data = labels, aes(x = -10, y = y+0.5, label = p), size = 2.8) +
    scale_fill_identity() +
    guides(fill="none") +
    ggtitle(tbl$sample[1])
  
})

ggsave('plot.pdf', plots[[1]], units = 'in', height = 8, width = 8, dpi = 300)




# 
# # Code for identifying position of elements within vectors to help build annotations
# #--------------------------------------------------------------------------------------------------------------
# library(Biostrings)
# 
# v <- readDNAStringSet('~/scratch/AAVengeR_3.0.7/data/vectors/Sabatino_PMC7855056_singleChain-plasmid.fasta')
# 
# vmatchPattern('CCCCTAAAATGGGCAAACATTGCAAGCAGCAAACAGCAAACACACAGCCCTCCCTGCCTGCTGACCTTGGAGCTGGGGCAGAGGTCAGAGACCTCTCTGGGCCCATGCCACCTCCAACATCCACTCGACCCCTTGGAATTTCGGTGGAGAGGAGCAGAGGTTGTCCTGGCGTGGTTTAGGTAGTGTGAGAGGGGAATGACTCCTTTCGGTAAGTGCAGTGGAAGCTGTACACTGCCCAGGCAAAGCGTCCGGGCAGCGTAGGCGGGCGACTCAGATCCCAGCCAGTGGACTTAGCCCCTGTTTGCTCCTCCGATAACTGGGGTGACCTTGGTTAATATTCACCAGCAGCCTCCCCCGTTGCCCCTCTGGATCCACTGCTTAAATACGGACGAGGACAGGGCCCTGTCTCCTCAGCTTCAGGCACCACCACTGACCTGGGACAGTGAATCCGGACTCTAAGGTAAATATAAAATTTTTAAGTGTATAATGTGTTAAACTACTGATTCTAATTGTTTCTCTCTTTTAG',
#               v, max.mismatch = 1)
# 
# 
# # Code for selecting rearranged species that have one or more pieces overlapping with a specific vector element.
# #----------------------------------------------------------------------------------------------------------------
# testFeature <- makeGRangesFromDataFrame(tibble(seqnames = 'X', start = 2482, end = 2685, strand = '*'))
# 
# rearrangement_table <- bind_rows(lapply(split(rearrangement_table, 1:nrow(rearrangement_table)), function(x){
#   g <- unlist(GRangesList(lapply(unlist(strsplit(x$shorthand, ';')), function(xx){
#     b <- sub('\\[', '', stringr::str_extract(xx, '\\[[^\\]]+'))
#     
#     if(! grepl('X', b, ignore.case = TRUE)){
#       o <- unlist(strsplit(b, '[\\+\\-]'))
#       return(makeGRangesFromDataFrame(tibble(seqnames = 'X', start = as.integer(o[1]), end = as.integer(o[2]), strand = '*')))
#     } else {
#       return(GRanges())
#     }
#   })))
#   
#   o <- findOverlaps(g, testFeature)
#   if(length(o) > 0){
#     return(x)
#   } else{
#     return(tibble())
#   }
# }))
# 
