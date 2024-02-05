suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(lubridate))

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')

opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

createOuputDir()
dir.create(file.path(opt$outputDir, opt$callNearestGenes_outputDir))

# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, opt$callNearestGenes_outputDir, 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)

updateLog('Starting callNearestGenes.')

if(! file.exists(file.path(opt$outputDir, opt$callNearestGenes_inputFile))){
  updateLog('Error - input file does not exist.')
  q(save = 'no', status = 1, runLast = FALSE)
}

sites <- readRDS(file.path(opt$outputDir, opt$callNearestGenes_inputFile))

if(nrow(sites) == 0){
  updateLog('Error - sites file was read and contains zero rows of data.')
  q(save = 'no', status = 1, runLast = FALSE)
}

if(! opt$callNearestGenes_addAfter %in% names(sites)){
  updateLog(paste0('Error - ', opt$callNearestGenes_addAfter, ' is not a column in your input data frame.'))
  q(save = 'no', status = 1, runLast = FALSE) 
}

sites <- distinct(bind_rows(lapply(split(sites, sites$refGenome), function(x){
           updateLog(paste0('Starting analysis of sites found in ', x$refGenome[1]))
  
           if(! x$refGenome[1] %in% names(opt$callNearestGenes_boundaries)){
             updateLog(paste0('   Error - ', x$refGenome[1], ' is not defined beneath callNearestGenes_boundaries in the config file.'))
             q(save = 'no', status = 1, runLast = FALSE) 
           }
  
           if(! 'TUs' %in% names(opt$callNearestGenes_boundaries[[x$refGenome[1]]])){
             updateLog(paste0('Error - TUs file not defined for ', x$refGenome[1], ' in the config file.'))
             q(save = 'no', status = 1, runLast = FALSE) 
           }
  
           if(! 'exons' %in% names(opt$callNearestGenes_boundaries[[x$refGenome[1]]])){
             updateLog(paste0('   Error - exons file not defined for ', x$refGenome[1], ' in the config file.'))
             q(save = 'no', status = 1, runLast = FALSE) 
           }
  
           if(! file.exists(file.path(opt$softwareDir, 'data', 'genomeAnnotations', opt$callNearestGenes_boundaries[[x$refGenome[1]]][['TUs']]))){
             updateLog(paste0('Error - TUs file set for ', x$refGenome[1], ' could not be found.'))
             q(save = 'no', status = 1, runLast = FALSE) 
           }
  
           if(! file.exists(file.path(opt$softwareDir, 'data', 'genomeAnnotations', opt$callNearestGenes_boundaries[[x$refGenome[1]]][['exons']]))){
             updateLog(paste0('   Error - exons file set for ', x$refGenome[1], ' could not be found.'))
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
          
          updateLog(paste0('Calling nearestGene() for ', n_distinct(posids), ' sites.'))
          n <- nearestGene(posids, genes, exons, CPUs = opt$callNearestGenes_CPUs)
          
          n$posid2 <- paste0(n$chromosome, n$strand, n$position)
          n <- n[!duplicated(n$posid2),]
          x$posid2 <- sub('\\.\\d+', '', x$posid)
          
          len <- length(x)
          n <- select(n, nearestGene, nearestGeneStrand, nearestGeneDist, inGene, inExon, beforeNearestGene, posid2)
          
          x <- left_join(x, n, by = 'posid2') %>% select(-posid2)
          
          dplyr::relocate(x, names(n)[! grepl('posid', names(n))], .after = opt$callNearestGenes_addAfter)
       })))

sites <- arrange(sites, desc(sonicLengths))

saveRDS(sites, file.path(opt$outputDir, opt$callNearestGenes_outputDir, 'sites.rds'))
openxlsx::write.xlsx(sites, file.path(opt$outputDir, opt$callNearestGenes_outputDir, 'sites.xlsx'))
readr::write_tsv(sites, file.path(opt$outputDir, opt$callNearestGenes_outputDir, 'sites.tsv.gz'))

updateLog('callNearestGenes completed.')

q(save = 'no', status = 0, runLast = FALSE) 
