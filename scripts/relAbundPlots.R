library(tidyverse)
library(RMySQL)
library(RColorBrewer)

args <- list()
args$numClones <- 10

sites <- readRDS('/data/project/Persaud_HIV/sites.rds')

invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 'select * from gtsp where Trial="Persaud HIV Integration"')
samples <- dplyr::rename(samples, timePoint = Timepoint, cellType = CellType)

sites <- left_join(sites, select(samples, SpecimenAccNum, timePoint, cellType), by = c("sample" = "SpecimenAccNum"))

sites <- group_by(sites, subject, sample) %>%
  mutate(relAbund = estAbund/sum(estAbund)) %>%
  ungroup()

ppNum <- function (n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

expandTimePoints <- function(tps){
  d <- tibble::tibble(tp = sub('_', '.', tps))
  d$n <- 1:nrow(d)
  
  d$timePointType <- stringr::str_match(base::toupper(d$tp), '[DMY]')
  d$timePointType[which(is.na(d$timePointType))] <- 'X'
  
  d <- dplyr::bind_rows(lapply(split(d, d$timePointType), function(x){
    n <- as.numeric(stringr::str_match(x$tp, '[\\d\\.]+')) * ifelse(grepl('\\-', x$tp), -1, 1)
    
    if(x$timePointType[1] == 'D'){
      x$timePointMonths <- base::round(n / 30.4167, digits = 0)
      x$timePointDays   <- base::round(n, digits = 0)
    } else if(x$timePointType[1] == 'M'){
      x$timePointMonths <- base::round(n, digits = 0)
      x$timePointDays   <- base::round(n * 30.4167, digits = 0)
    } else if(x$timePointType[1] == 'Y'){
      x$timePointMonths <- base::round(n * 12, digits = 0)
      x$timePointDays   <- base::round(n * 365, digits = 0)
    } else {
      message('Warning - could not determine date unit for: ', paste0(unique(x$timePoint), collapse = ', '))
      x$timePointMonths <- n
      x$timePointDays   <- n 
    }
    x
  }))
  
  data.frame(dplyr::arrange(d, n) %>% dplyr::select(timePointMonths, timePointDays))
}


sites$posidLabel <- sites$posid
sites$timePointDays <- expandTimePoints(sites$timePoint)$timePointDays

sites <- sites[! is.na(sites$timePointDays),]

# Create data frame needed to generate relative abundance plots.
abundantClones <- bind_rows(lapply(split(sites, sites$subject), function(s){
  
    o <- bind_rows(lapply(split(s, s$cellType), function(x){
      
     # Adjust the number of clones to return based on the number of sites per cell type.
     if(nrow(x) < args$numClones) args$numClones <- nrow(x)
  
     # Sort nearest genes by abundance.
     x <- x[order(x$estAbund, decreasing = TRUE),]
  
     # Select clones to report.
     topClones <-  unique(x$posidLabel)[1:args$numClones]
  
     # For each time point, create a data frame for relative abundance plots
     bind_rows(lapply(split(x, x$timePoint), function(x2){
    
     lowAbundData <- dplyr::mutate(x2, posidLabel = 'LowAbund',
                                   totalCells = sum(estAbund),
                                   relAbund   = 1) %>%
                    dplyr::slice(1) %>% 
                    dplyr::select(cellType, timePoint, posidLabel, totalCells, timePointDays, relAbund)
    
      x3 <- subset(x2, posidLabel %in% topClones)
      if(nrow(x3) == 0) return(lowAbundData)
      x3$totalCells <- sum(x2$estAbund)
    
      lowAbundData$relAbund <- 1 - sum(x3$relAbund)
     bind_rows(lowAbundData,  dplyr::select(x3, cellType, timePoint, posidLabel, totalCells, timePointDays, relAbund))
    }))
  }))
  
    o$subject <- s$subject[1]
    o
}))    







abundantClonesPlots <- 
  lapply(split(abundantClones, abundantClones$subject), function(s){
    
    # Create named color vector for unique clones.
    cloneColorsVector <- setNames(c('#eeeeee', colorRampPalette(brewer.pal(12, "Paired"))(n_distinct(s$posidLabel))),  c('LowAbund', unique(s$posidLabel)))
    
    r <- lapply(split(s, s$cellType), function(x){
           o <- subset(x, posidLabel != 'LowAbund')
           o <- o[order(o$relAbund, decreasing = TRUE),]
  
           x$posidLabel <- factor(x$posidLabel, levels = c('LowAbund', unique(o$posidLabel)))
           x <- x[order(x$timePointDays),]
           x$timePoint  <- factor(x$timePoint, levels = unique(subset(abundantClones, subject == x$subject[1])$timePoint))   # Unique subject level time points
  
           totalCellLabel <- unname(unlist(lapply(split(x, x$timePoint), function(x) ppNum(x$totalCells[1]))))
  
           ggplot(x) +
           theme_bw() +
           scale_x_discrete(drop=FALSE) + 
           geom_bar(aes(timePoint, relAbund, fill=posidLabel), stat='identity', color = 'black', size = 0.20) + 
           scale_fill_manual(name = 'Clones', values = cloneColorsVector) +
           scale_shape_manual(values = c(16, 17, 15), drop = FALSE) +
           labs(x = 'Timepoint', y = 'Relative Sonic Abundance') +
           ggtitle(x$cellType[1]) +
           guides(fill=guide_legend(title.position = "top", ncol=1, keyheight=0.35, default.unit="inch")) +
           scale_y_continuous(labels = scales::percent) + 
           annotate('text', x=1:length(totalCellLabel), y=1.04, label=totalCellLabel, size=2.7, angle=45, hjust=0.5) +
           theme(axis.text.x = element_text(angle = 315, hjust = 0))
          })
    
    names(r) <- unique(s$cellType)
    r
})

