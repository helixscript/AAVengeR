library(dplyr)
library(lubridate)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)

source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$callNearestGenes_outputDir))

sites <- readRDS(file.path(opt$outputDir, opt$callNearestGenes_inputFile))

# Add refGenome to sites.
sites <- distinct(left_join(sites, select(samples, uniqueSample, refGenome.id), by = 'uniqueSample'))

sites <- distinct(bind_rows(lapply(split(sites, sites$refGenome.id), function(x){
  
           genes <- file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x$refGenome.id[1], '.TUs.rds'))
           if(! file.exists(genes)){
             write(c(paste(now(), "Error - could not find data file: ", genes)), file = file.path(opt$outputDir, 'log'), append = TRUE)
             q(save = 'no', status = 1, runLast = FALSE) 
           }
           genes <- readRDS(genes)
       
           exons <- file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x$refGenome.id[1], '.exons.rds'))
           if(! file.exists(exons)){
             write(c(paste(now(), "Error - could not find data file: ", exons)), file = file.path(opt$outputDir, 'log'), append = TRUE)
             q(save = 'no', status = 1, runLast = FALSE) 
           }
           exons <- readRDS(exons)
  
           if(opt$callNearestGenes_use_humanXeno){
             genes <- file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x$refGenome.id[1], '.humanXeno.TUs.rds'))
             if(grepl('hg\\d+', x$refGenome.id[1])) genes <- sub('\\.humanXeno', '', genes) # Undo adding humanXeno to a human genome.
         
             if(! file.exists(genes)){
               write(c(paste(now(), "Error - could not find data file: ", genes)), file = file.path(opt$outputDir, 'log'), append = TRUE)
               q(save = 'no', status = 1, runLast = FALSE) 
             }
             
             genes <- readRDS(genes)

             exons <- file.path(opt$softwareDir, 'data', 'genomeAnnotations', paste0(x$refGenome.id[1], '.humanXeno.exons.rds'))
             if(grepl('hg\\d+', x$refGenome.id[1])) exons <- sub('\\.humanXeno', '', exons)  # Undo adding humanXeno to a human genome.
         
             if(! file.exists(exons)){
               write(c(paste(now(), "Error - could not find data file: ", exons)), file = file.path(opt$outputDir, 'log'), append = TRUE)
               q(save = 'no', status = 1, runLast = FALSE) 
             }
             
             exons <- readRDS(exons)
           }
  
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
          
          posids <- x$posid
          posids <- unique(sub('\\.\\d+$', '', posids))
          
          n <- nearestGene(posids, genes, exons, CPUs = opt$callNearestGenes_CPUs)
 
          n$posid2 <- paste0(n$chromosome, n$strand, n$position)
          n <- n[!duplicated(n$posid2),]
          x$posid2 <- sub('\\.\\d+', '', x$posid)
          
          left_join(x, select(n, nearestGene, nearestGeneStrand, nearestGeneDist, inGene, inExon, beforeNearestGene, posid2), by = 'posid2') %>% select(-posid2)
       })))

saveRDS(select(sites, -uniqueSample, refGenome.id), file.path(opt$outputDir, opt$callNearestGenes_outputDir, opt$callNearestGenes_outputFile))
openxlsx::write.xlsx(select(sites, -uniqueSample, refGenome.id), file.path(opt$outputDir, opt$callNearestGenes_outputDir, 'sites.xlsx'))

q(save = 'no', status = 0, runLast = FALSE) 


