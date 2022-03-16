library(dplyr)
library(lubridate)

# configFile <- 'previous.yml'
configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)

source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$callNearestGenes_outputDir))

sites <- readRDS(file.path(opt$outputDir, opt$callNearestGenes_inputFile))

if(! 'posid' %in% names(sites)) stop('Error -- posid not found in sites data.')

samples <- loadSamples()

# Need to add trial ids to sites.
samples$n <- paste0(samples$trial, '~', samples$subject, '~', samples$sample)
sites$n <- paste0(samples$trial, '~', sites$subject, '~', sites$sample)

sites <- left_join(sites, select(samples, n, refGenome.id), by = 'n')

sites <- distinct(bind_rows(lapply(split(sites, sites$refGenome.id), function(x){
  
           ### if(x$refGenome.id[1] == 'hg38') return(tibble())
  
           genes <- file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x$refGenome.id[1], '.TUs.rds'))
           if(! file.exists(genes)) stop(paste0('Error -- could not file ', genes))
       
           exons <- file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x$refGenome.id[1], '.exons.rds'))
           if(! file.exists(exons)) stop(paste0('Error -- could not file ', exons))
  
           if(opt$callNearestGenes_use_humanXeno){
             genes <- file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x$refGenome.id[1], '.humanXeno.TUs.rds'))
             if(grepl('hg\\d+', x$refGenome.id[1])) genes <- sub('\\.humanXeno', '', genes) # Undo adding humanXeno to a human genome.
         
             if(! file.exists(genes)) stop(paste0('Error -- could not file ', genes))

             exons <- file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x$refGenome.id[1], '.humanXeno.exons.rds'))
             if(grepl('hg\\d+', x$refGenome.id[1])) exons <- sub('\\.humanXeno', '', exons)  # Undo adding humanXeno to a human genome.
         
             if(! file.exists(exons)) stop(paste0('Error -- could not file ', exons))
           }
  
           message('loading ', genes); genes <- readRDS(genes)
           message('loading ', exons); exons <- readRDS(exons)
  
           # Collapse TUs to single positions if requested. 
           # Remove exons since they would not be relevant.
           if(opt$callNearestGenes_TU_position != 'either'){
             options(warn=-1)
             if (subjectSide == "start") genes <- GenomicRanges::flank(genes, width = -1)
             if (subjectSide == "end")   genes <- GenomicRanges::flank(genes, width = -1, start = FALSE)
             if (subjectSide == "center") ranges(genes) <- IRanges(mid(ranges(genes)), width = 1)
             options(warn=0)
             exons <- GenomicRanges::GRanges()
           }
  
          message(x$refGenome.id[1])
          n <- nearestGene(unique(sub('\\.\\d+$', '', x$posid)), genes, exons, CPUs = opt$callNearestGenes_CPUs)
 
          n$posid2 <- paste0(n$chromosome, n$strand, n$position)
          n <- n[!duplicated(n$posid2),]
          x$posid2 <- sub('\\.\\d+', '', x$posid)
          
          left_join(x, select(n, nearestGene, nearestGeneStrand, nearestGeneDist, inGene, inExon, beforeNearestGene, posid2), by = 'posid2') %>% select(-n, -posid2)
       })))

saveRDS(sites, file.path(opt$outputDir, opt$callNearestGenes_outputDir, opt$callNearestGenes_outputFile))

q(save = 'no', status = 0, runLast = FALSE) 


