library(dplyr)
library(parallel)
library(pander)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

# Create directories and start logging.
dir.create(opt$outputDir)
dir.create(file.path(opt$outputDir, 'core'))
dir.create(file.path(opt$outputDir, 'core', 'demultiplex'))
dir.create(file.path(opt$outputDir, 'core', 'replicate_analyses'))
dir.create(file.path(opt$outputDir, 'core', 'subject_analyses'))           
          
write(paste0(date(), ' - start'), file = file.path(opt$outputDir, 'core', 'log'))

# Gather select sections from the configuration file.
o <- opt[grepl('^mode|^compressDataFiles|^Rscript|^softwareDir|^outputDir|^databaseGroup|^core|^demultiplex', names(opt))]
o$outputDir <- file.path(opt$outputDir, 'core')

# Run demultiplex module.
o$demultiplex_CPUs <- opt$core_CPUs
yaml::write_yaml(o, file.path(opt$outputDir, 'core',  'demultiplex', 'config.yml'))

write(c('#!/usr/bin/sh', 
        paste(opt$Rscript, file.path(opt$softwareDir, 'demultiplex.R'), file.path(opt$outputDir, 'core',  'demultiplex', 'config.yml'))), 
      file = file.path(opt$outputDir, 'core',  'demultiplex', 'run.sh'))

system(paste('chmod 755', file.path(opt$outputDir, 'core',  'demultiplex', 'run.sh')))
system(file.path(opt$outputDir, 'core',  'demultiplex', 'run.sh'), wait = TRUE, show.output.on.console = FALSE)

# Read in demultiplex result.
reads <- readRDS(file.path(opt$outputDir, 'core', 'demultiplex', 'reads.rds'))

jobStatus <- function(searchPattern = 'sites.done', outputFileName = 'jobTable.txt'){
  f <- list.files(file.path(opt$outputDir, 'core'), pattern = searchPattern, full = TRUE, recursive = TRUE)
  
  if(length(f) > 0){
    invisible(sapply(f, function(x){
      o <- unlist(strsplit(x, '/'))
      id <- o[length(o)-2]
      
      jobTable[jobTable$id == id,]$done     <<- TRUE
      jobTable[jobTable$id == id,]$active   <<- NA
      jobTable[jobTable$id == id,]$endTime  <<- as.character(lubridate::now()) 
      jobTable[jobTable$id == id,]$duration <<- paste0(as.integer(as.character(difftime(lubridate::now(), lubridate::ymd_hms(jobTable[jobTable$id == o[length(o)-2],]$startTime), units="mins"))), ' minutes')
      
      CPUs_used <<- CPUs_used - jobTable[jobTable$id == id,]$CPUs
      invisible(file.remove(f))
    }))
  }
  
  # Redefine CPUs.
  o <- jobTable[jobTable$done == FALSE & jobTable$active == FALSE,]
  
  if(nrow(o) > 0){
    # Increase core_maxPercentCPUs near the end of the job table.
    core_maxPercentCPUs <- opt$core_maxPercentCPUs
    if(100/nrow(o) > core_maxPercentCPUs) core_maxPercentCPUs <- ceiling(100/nrow(o)) 
  
    # The job with the most reads will evaluate to 1 * core_maxPercentCPUs and recieve core_maxPercentCPUs cores.
    o$CPUs <- as.integer(ceiling((((o$reads / max(o$reads)) * core_maxPercentCPUs) / 100) * opt$core_CPUs))
  
    o$CPUs <- ifelse(o$CPUs == 1, 2, o$CPUs) # Each job must get at least 2 CPUs.
  
    jobTable[match(o$id, jobTable$id),]$CPUs <- o$CPUs
  }
  
  tab <- pandoc.table.return(arrange(jobTable, desc(done), desc(active)), style = "simple", split.tables = Inf, plain.ascii = TRUE)
  write(c(date(), paste0('Available CPUs: ', (opt$core_CPUs - CPUs_used)), tab), file = file.path(opt$outputDir, 'core', outputFileName), append = FALSE)
}

# Build a jobs table for prepReads, alignReads, and buildFragments.
jobTable <- data.frame(table(reads$uniqueSample)) %>% dplyr::rename(id = Var1) %>% dplyr::rename(reads = Freq)
jobTable$CPUs <- as.integer(ceiling((((jobTable$reads / max(jobTable$reads)) * opt$core_maxPercentCPUs) / 100) * opt$core_CPUs))
jobTable$CPUs <- ifelse(jobTable$CPUs == 1, 2, jobTable$CPUs)
jobTable$active <- FALSE
jobTable$startTime <- NA
jobTable$endTime <- NA
jobTable$done <- FALSE
jobTable$duration <- NA

CPUs_used <- 0

# Run prepReads, alignReads, and buildFragments.
while(! all(jobTable$done == TRUE)){
  
  CPUs_available <- opt$core_CPUs - CPUs_used
  tab <- jobTable[jobTable$CPUs <= CPUs_available & jobTable$active == FALSE & jobTable$done == FALSE,] %>% dplyr::slice_max(CPUs, with_ties = FALSE)
  
  if(nrow(tab) == 0){
    jobStatus(searchPattern = 'fragments.done', outputFileName = 'replicateJobTable')
    Sys.sleep(5)
    next
  }
  
  # Gather select sections from the configuration file and set CPU values.
  o <- opt[grepl('^mode|^compressDataFiles|^Rscript|^softwareDir|^outputDir|^databaseGroup|^core|^prepReads|^alignReads|^buildFragments', names(opt))]
  o$prepReads_CPUs <- as.integer(tab$CPUs)
  o$alignReads_CPUs <- as.integer(tab$CPUs)
  o$buildFragments_CPUs <- as.integer(tab$CPUs)
  o$outputDir <- file.path(opt$outputDir, 'core', 'replicate_analyses', tab$id)
  
  # Instruct the pipeline to create a buildFragments/fragments.done file for jobs that failed.
  o$core_createFauxFragDoneFiles <- TRUE
  
  # Create the output directory for this replicate and save a copy of the replicate reads within it.
  dir.create(o$outputDir)
  saveRDS(subset(reads, uniqueSample == tab$id), file.path(o$outputDir, paste0(tab$id, '.CORE_TMP.rds')), compress = opt$compressDataFiles)
  
  # Update configuration to point to the subset of replicate reads.
  o$prepReads_readsTable <- paste0(tab$id, '.CORE_TMP.rds')
  
  # Define modules to run in replicate level configuration file and write out file.
  o$modules <- list()
  o[['modules']] <- c('prepReads', 'alignReads', 'buildFragments')
  yaml::write_yaml(o, file.path(opt$outputDir, 'core',  'replicate_analyses', tab$id, 'config.yml'))
  
  # Create AAVengeR launching script since system(x, wait = TRUE) does not wait for commands w/ arguments.
  write(c('#!/usr/bin/sh', 
          paste(opt$Rscript, file.path(opt$softwareDir, 'aavenger.R'), file.path(opt$outputDir, 'core',  'replicate_analyses', tab$id, 'config.yml'))), 
        file = file.path(opt$outputDir, 'core',  'replicate_analyses', tab$id, 'run.sh'))
  
  system(paste('chmod 755', file.path(opt$outputDir, 'core',  'replicate_analyses', tab$id, 'run.sh')))
  system(file.path(opt$outputDir, 'core', 'replicate_analyses', tab$id, 'run.sh'), wait = FALSE, show.output.on.console = FALSE)
  
  CPUs_used <- CPUs_used + tab$CPUs
  
  jobTable[which(jobTable$id == tab$id),]$active <- TRUE
  jobTable[which(jobTable$id == tab$id),]$startTime <- as.character(lubridate::now())

  jobStatus(searchPattern = 'fragments.done', outputFileName = 'replicateJobTable')
  Sys.sleep(5)
}

invisible(file.remove(list.files(file.path(opt$outputDir, 'core', 'replicate_analyses'), recursive = TRUE, pattern = 'CORE_TMP', full.names = TRUE)))

# Bundle together fragment output files.
f <- list.files(file.path(opt$outputDir, 'core'), pattern = 'fragments.rds', recursive = TRUE, full.names = TRUE)
frags <- bind_rows(lapply(f, readRDS))

# Create a subject level to-do table.
jobTable <- tibble(u = unique(frags$uniqueSample))
jobTable$id <- sapply(jobTable$u, function(x) paste0(unlist(strsplit(x, '~'))[1:2], collapse = '~'))

# Collate logs.
write('Replicate level modules:\n', file = file.path(opt$outputDir, 'core', 'log'), append = TRUE)

for(x in list.files(file.path(opt$outputDir, 'core', 'replicate_analyses'))){
  write(paste0('\n', x), file = file.path(opt$outputDir, 'core', 'log'), append = TRUE)
  
  if(file.exists(file.path(opt$outputDir, 'core', 'replicate_analyses', x, 'prepReads', 'log'))){
   write(paste0(' \n[prepReads]'), file = file.path(opt$outputDir, 'core', 'log'), append = TRUE)
   write(paste0(' ',readLines(file.path(opt$outputDir, 'core', 'replicate_analyses', x, 'prepReads', 'log'))), file = file.path(opt$outputDir, 'core', 'log'), append = TRUE)
  }
  
  if(file.exists(file.path(opt$outputDir, 'core', 'replicate_analyses', x, 'alignReads', 'log'))){
   write(paste0(' \n[alignReads]'), file = file.path(opt$outputDir, 'core', 'log'), append = TRUE)
   write(paste0(' ',readLines(file.path(opt$outputDir, 'core', 'replicate_analyses', x, 'alignReads', 'log'))), file = file.path(opt$outputDir, 'core', 'log'), append = TRUE)
  }
  
  if(file.exists(file.path(opt$outputDir, 'core', 'replicate_analyses', x, 'buildFragments', 'log'))){
  # write(paste0(' \n[buildFragments]'), file = file.path(opt$outputDir, 'core', 'log'), append = TRUE)
  # write(paste0(' ',readLines(file.path(opt$outputDir, 'core', 'replicate_analyses', x, 'buildFragments', 'log'))), file = file.path(opt$outputDir, 'core', 'log'), append = TRUE)
  }
}

# Create a subject-level splitting vector for fragments records.
frags <- left_join(frags, jobTable, by = c('uniqueSample' = 'u'))

jobTable <- group_by(distinct(frags), id) %>% summarise(reads = n()) %>% ungroup() %>% arrange(desc(reads))

jobTable$CPUs <- as.integer(ceiling((((jobTable$reads / max(jobTable$reads)) * opt$core_maxPercentCPUs) / 100) * opt$core_CPUs))
jobTable$CPUs <- ifelse(jobTable$CPUs == 1, 2, jobTable$CPUs)
jobTable$active <- FALSE
jobTable$startTime <- NA
jobTable$endTime <- NA
jobTable$done <- FALSE
jobTable$duration <- NA

### frags$trialSubject <- sub('~[^~]+~\\d+$', '', frags$uniqueSample)

d <- tibble(uniqueSample = unique(frags$uniqueSample)) %>% tidyr::separate(uniqueSample, c('trial', 'subject', 'sample', 'replicate'), sep = '~', remove = FALSE)
d$trialSubject <- paste0(d$trial, '~', d$subject)
frags <- left_join(frags, distinct(dplyr::select(d, uniqueSample, trialSubject)), by = 'uniqueSample')


CPUs_used <- 0

while(! all(jobTable$done == TRUE)){
  
  CPUs_available <- opt$core_CPUs - CPUs_used
  tab <- jobTable[jobTable$CPUs <= CPUs_available & jobTable$active == FALSE & jobTable$done == FALSE,] %>% dplyr::slice_max(CPUs, with_ties = FALSE)
  
  if(nrow(tab) == 0){
    jobStatus(searchPattern = 'sites.done', outputFileName = 'subjectJobTable')
    Sys.sleep(5)
    next
  }

  # Gather select sections from the configuration file and set CPU values.
  o <- opt[grepl('^mode|^compressDataFiles|^Rscript|^softwareDir|^outputDir|^databaseGroup|^core|^buildStdFragments|^buildSites', names(opt))]
  o$buildStdFragments_CPUs <- as.integer(tab$CPUs)
  o$buildSites_CPUs <- as.integer(tab$CPUs)
  o$outputDir <- file.path(opt$outputDir, 'core', 'subject_analyses', tab$id)
  
  # Instruct the pipeline to create a buildFragments/fragments.done file for jobs that failed.
  o$core_createFauxSiteDoneFiles <- TRUE
  
  # Create the output directory for this replicate and save a copy of the replicate reads within it.
  dir.create(o$outputDir)
  
  saveRDS(subset(frags, trialSubject == tab$id), file.path(o$outputDir, paste0(tab$id, '.CORE_TMP.rds')), compress = opt$compressDataFiles)
  
  # Update configuration to point to the subset of replicate reads.
  o$buildStdFragments_inputFile <- file.path(paste0(tab$id, '.CORE_TMP.rds'))
  
  # Define modules to run in replicate level configuration file and write out file.
  o$modules <- list()
  o[['modules']] <- c('buildStdFragments', 'buildSites')
  
  yaml::write_yaml(o, file.path(opt$outputDir, 'core', 'subject_analyses', tab$id, 'config.yml'))
  
  # Create AAVengeR launching script since system(x, wait = TRUE) does not wait for commands w/ arguments.
  write(c('#!/usr/bin/sh', 
          paste(opt$Rscript, file.path(opt$softwareDir, 'aavenger.R'), file.path(opt$outputDir, 'core', 'subject_analyses',  tab$id, 'config.yml'))), 
        file = file.path(opt$outputDir, 'core', 'subject_analyses',  tab$id, 'run.sh'))
  
  system(paste('chmod 755', file.path(opt$outputDir, 'core', 'subject_analyses', tab$id, 'run.sh')))
  system(file.path(opt$outputDir, 'core', 'subject_analyses', tab$id, 'run.sh'), wait = FALSE, show.output.on.console = FALSE)
  
  CPUs_used <- CPUs_used + tab$CPUs
  
  jobTable[which(jobTable$id == tab$id),]$active <- TRUE
  jobTable[which(jobTable$id == tab$id),]$startTime <- as.character(lubridate::now())
  
  jobStatus(searchPattern = 'sites.done', outputFileName = 'subjectJobTable')
  Sys.sleep(5)
}

# Collate logs.
write('Subject level modules:\n', file = file.path(opt$outputDir, 'core', 'log'), append = TRUE)

for(x in list.files(file.path(opt$outputDir, 'core', 'subject_analyses'))){
  write(paste0('\n', x), file = file.path(opt$outputDir, 'core', 'log'), append = TRUE)
  
  if(file.exists(file.path(opt$outputDir, 'core', 'subject_analyses', x, 'buildStdFragments', 'log'))){
    write(paste0(' \n[buildStdFragments]'), file = file.path(opt$outputDir, 'core', 'log'), append = TRUE)
    write(paste0(' ',readLines(file.path(opt$outputDir, 'core', 'subject_analyses', x, 'buildStdFragments', 'log'))), file = file.path(opt$outputDir, 'core', 'log'), append = TRUE)
  }
  
  if(file.exists(file.path(opt$outputDir, 'core', 'subject_analyses', x, 'buildSites', 'log'))){
    write(paste0(' \n[buildSites]'), file = file.path(opt$outputDir, 'core', 'log'), append = TRUE)
    write(paste0(' ',readLines(file.path(opt$outputDir, 'core', 'subject_analyses', x, 'buildSites', 'log'))), file = file.path(opt$outputDir, 'core', 'log'), append = TRUE)
  }
}

invisible(file.remove(list.files(file.path(opt$outputDir, 'core', 'subject_analyses'), recursive = TRUE, pattern = 'CORE_TMP', full.names = TRUE)))

# Bundle together multi-hit cluster output files.
f <- list.files(file.path(opt$outputDir, 'core'), pattern = 'multiHitClusters.rds', recursive = TRUE, full.names = TRUE)
multiHitClusters <- bind_rows(lapply(f, readRDS))
saveRDS(multiHitClusters, file = file.path(opt$outputDir, 'core', 'multiHitClusters.rds'), compress = opt$compressDataFiles)

# Bundle together site output files.
f <- list.files(file.path(opt$outputDir, 'core'), pattern = 'sites.rds', recursive = TRUE, full.names = TRUE)
sites <- bind_rows(lapply(f, readRDS))
saveRDS(sites, file = file.path(opt$outputDir, 'core', 'sites.rds'), compress = opt$compressDataFiles)
openxlsx::write.xlsx(sites, file.path(opt$outputDir, 'core', 'sites.xlsx'))
readr::write_tsv(sites, file.path(opt$outputDir, 'core', 'sites.tsv.gz'))

write(paste0(date(), ' - end'), file = file.path(opt$outputDir, 'core', 'log'), append = TRUE)
