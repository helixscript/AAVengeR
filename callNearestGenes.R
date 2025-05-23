#!/usr/bin/Rscript

# AAVengeR/buildStdFragments.R
# John K. Everett, Ph.D.
# 
# This script accepts input from buildSites or modules following it and annotates
# nearest genomic features. Details such as distance to nearest features and their
# orientation are included.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(GenomicRanges))

# Read in the configuration file and perform basic sanity checks.
args <- commandArgs(trailingOnly=TRUE)

if(length(args) == 0) stop('Expected at least one command line argument')
source(file.path(yaml::read_yaml(args[1])$softwareDir, 'lib.R'))
opt <- startModule(args)

createOuputDir()
dir.create(file.path(opt$outputDir, opt$callNearestGenes_outputDir), showWarnings = FALSE)

# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, opt$callNearestGenes_outputDir, 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)
write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)

quitOnErorr <- function(msg){
  updateLog(msg)
  updateLog(paste0('See log for more details: ', opt$defaultLogFile))
  q(save = 'no', status = 1, runLast = FALSE) 
}

updateLog(paste0('Starting callNearestGenes using ', opt$callNearestGenes_CPUs, ' CPUs.'))
updateMasterLog()

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
           updateLog(paste0('Starting analysis of sites found in ', x$refGenome[1], '.'))
           updateMasterLog()
  
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
           
           if(tolower(opt$callNearestGenes_geneList) != 'none'){
             updateLog(paste0('Filtering gene ranges with gene name list: ', opt$callNearestGenes_geneList)) 
             geneList <- toupper(readLines(opt$callNearestGenes_geneList))
             
             if(! file.exists(opt$callNearestGenes_geneList)){
               quitOnErorr(paste0('Gene filter list (', opt$callNearestGenes_geneList, ') could not be found.'))
             }
             
             genes <- genes[toupper(genes$name2) %in% geneList]
             exons <- exons[toupper(exons$name2) %in% geneList]
           }
           
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
          
          updateLog(paste0('Calling nearestGene() for ', ppNum(n_distinct(posids)), ' sites.'))
          updateMasterLog()
          n <- nearestGene(posids, genes, exons, CPUs = opt$callNearestGenes_CPUs)
          
          n$posid2 <- paste0(n$chromosome, n$strand, n$position)
          n <- n[!duplicated(n$posid2),]
          x$posid2 <- sub('\\.\\d+$', '', x$posid)
          
          len <- length(x)
          n <- select(n, nearestGene, nearestGeneStrand, nearestGeneDist, inGene, inExon, beforeNearestGene, posid2)
          
          if(tolower(opt$callNearestGenes_columnPrefix) != 'none'){
            names(n) <- paste0(opt$callNearestGenes_columnPrefix, names(n))
            x <- left_join(x, n, by = c('posid2' = paste0(opt$callNearestGenes_columnPrefix, 'posid2')), relationship = 'many-to-many') %>% select(-posid2)
          } else {
            x <- left_join(x, n, by = 'posid2', relationship = 'many-to-many') %>% select(-posid2)
          }
          
          dplyr::relocate(x, names(n)[! grepl('posid', names(n))], .after = opt$callNearestGenes_addAfter)
       })))

sites <- arrange(sites, desc(sonicLengths))

if(tolower(opt$database_configGroup) != 'none'){
  suppressPackageStartupMessages(library(RMariaDB))
  uploadSitesToDB(sites)
}

saveRDS(sites, file.path(opt$outputDir, opt$callNearestGenes_outputDir, 'sites.rds'))
openxlsx::write.xlsx(sites, file.path(opt$outputDir, opt$callNearestGenes_outputDir, 'sites.xlsx'))
readr::write_tsv(sites, file.path(opt$outputDir, opt$callNearestGenes_outputDir, 'sites.tsv.gz'))

updateLog('callNearestGenes completed.')
updateMasterLog()

q(save = 'no', status = 0, runLast = FALSE) 
