library(dplyr)
library(lubridate)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)

source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$callNearestGenes_outputDir))

sites <- readRDS(file.path(opt$outputDir, opt$callNearestGenes_inputFile))

sites <- distinct(bind_rows(lapply(split(sites, sites$refGenome), function(x){

           if(! x$refGenome[1] %in% names(opt$callNearestGenes_boundaries)){
             write(c(paste(now(), paste0('   Error - ', x$refGenome[1], ' is not defined beneath callNearestGenes_boundaries in the config file.'))), file = file.path(opt$outputDir, 'log'), append = TRUE)
             q(save = 'no', status = 1, runLast = FALSE) 
           }
  
           if(! 'TUs' %in% names(opt$callNearestGenes_boundaries[[x$refGenome[1]]])){
             write(c(paste(now(), paste0('   Error - TUs file not defined for ', x$refGenome[1], ' in the config file.'))), file = file.path(opt$outputDir, 'log'), append = TRUE)
             q(save = 'no', status = 1, runLast = FALSE) 
           }
  
           if(! 'exons' %in% names(opt$callNearestGenes_boundaries[[x$refGenome[1]]])){
             write(c(paste(now(), paste0('   Error - exons file not defined for ', x$refGenome[1], ' in the config file.'))), file = file.path(opt$outputDir, 'log'), append = TRUE)
             q(save = 'no', status = 1, runLast = FALSE) 
           }
  
           if(! file.exists(file.path(opt$softwareDir, 'data', 'genomeAnnotations', opt$callNearestGenes_boundaries[[x$refGenome[1]]][['TUs']]))){
             write(c(paste(now(), paste0('   Error - TUs file set for ', x$refGenome[1], ' could not be found.'))), file = file.path(opt$outputDir, 'log'), append = TRUE)
             q(save = 'no', status = 1, runLast = FALSE) 
           }
  
           if(! file.exists(file.path(opt$softwareDir, 'data', 'genomeAnnotations', opt$callNearestGenes_boundaries[[x$refGenome[1]]][['exons']]))){
             write(c(paste(now(), paste0('   Error - exons file set for ', x$refGenome[1], ' could not be found.'))), file = file.path(opt$outputDir, 'log'), append = TRUE)
             q(save = 'no', status = 1, runLast = FALSE) 
           }
  
  
           genes <- readRDS(file.path(opt$softwareDir, 'data', 'genomeAnnotations', opt$callNearestGenes_boundaries[[x$refGenome[1]]][['TUs']]))
           exons <- readRDS(file.path(opt$softwareDir, 'data', 'genomeAnnotations', opt$callNearestGenes_boundaries[[x$refGenome[1]]][['exons']]))
  
           # Collapse TUs to single positions if requested. 
           # Remove exons since they would not be relevant.
           if(opt$callNearestGenes_TU_position != 'either'){
             options(warn=-1)
             if (opt$callNearestGenes_TU_position == "start")  genes <- GenomicRanges::flank(genes, width = -1)
             if (opt$callNearestGenes_TU_position == "end")    genes <- GenomicRanges::flank(genes, width = -1, start = FALSE)
             if (opt$callNearestGenes_TU_position == "center") ranges(genes) <- IRanges(mid(ranges(genes)), width = 1)
             options(warn=0)
             exons <- GenomicRanges::GRanges()
           }
          
          posids <- x$posid
          posids <- unique(sub('\\.\\d+$', '', posids))
          
          message('Calling nearestGene(), genes: ', length(genes), ', exons: ', length(exons))
          
          n <- nearestGene(posids, genes, exons, CPUs = opt$callNearestGenes_CPUs)
 
          message('nearestGene() done')
          
          n$posid2 <- paste0(n$chromosome, n$strand, n$position)
          n <- n[!duplicated(n$posid2),]
          x$posid2 <- sub('\\.\\d+', '', x$posid)
          
          len <- length(x)
          n <- select(n, nearestGene, nearestGeneStrand, nearestGeneDist, inGene, inExon, beforeNearestGene, posid2)
          
          x <- left_join(x, n, by = 'posid2') %>% select(-posid2)
          
          dplyr::relocate(x, names(n)[! grepl('posid', names(n))], .after = 'repLeaderSeq')
       })))

saveRDS(sites, file.path(opt$outputDir, opt$callNearestGenes_outputDir, opt$callNearestGenes_outputFile))
openxlsx::write.xlsx(sites, file.path(opt$outputDir, opt$callNearestGenes_outputDir, 'sites.xlsx'))

q(save = 'no', status = 0, runLast = FALSE) 
