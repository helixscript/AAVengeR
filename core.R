# AAVengeR/core.R
# John K. Everett, Ph.D.
# 
# This script runs AAVengeR's six core modules:
# demultiplex, prepReads, alignReads, buildFragments, buildStdFragments, and buildSites.
# It is advantageous to use this module rather than calling the core modules serially
# for large data sets since the core module dynamically sets CPU thresholds for each 
# sample replicate according to the number of associated reads.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))

set.seed(1)

source(file.path(yaml::read_yaml(commandArgs(trailingOnly=TRUE)[1])$softwareDir, 'lib.R'))
opt <- loadConfig()
optionsSanityCheck()

# Create directories and start logging.
createOuputDir()
dir.create(file.path(opt$outputDir, 'core'))
dir.create(file.path(opt$outputDir, 'core', 'demultiplex'))
dir.create(file.path(opt$outputDir, 'core', 'replicate_analyses'))
dir.create(file.path(opt$outputDir, 'core', 'subject_analyses')) 

# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, 'core', 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)
write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)

quitOnErorr <- function(msg){
  updateLog(msg)
  message(msg)
  message(paste0('See log for more details: ', opt$defaultLogFile))
  q(save = 'no', status = 1, runLast = FALSE) 
}

opt$calledFromCore <- TRUE

updateLog('Starting core module.')
updateLog('This module will call the following modules:')
updateLog('demultiplex, prepReads, alignReads, buildFragments, buildStdFragments, and buildSites.')


# Gather select sections from the configuration file.
o <- opt
o$outputDir <- file.path(opt$outputDir, 'core')

# Run demultiplex module.
o$demultiplex_CPUs <- opt$core_CPUs
yaml::write_yaml(o, file.path(opt$outputDir, 'core',  'demultiplex', 'config.yml'))


# Create a shell script to start demultiplex module.
write(c('#!/usr/bin/sh', 
        paste(opt$Rscript, file.path(opt$softwareDir, 'demultiplex.R'), file.path(opt$outputDir, 'core',  'demultiplex', 'config.yml')),
        'echo $?'), 
        file = file.path(opt$outputDir, 'core',  'demultiplex', 'run.sh'))

updateLog('Running the demultiplex module.')

# All core CPUs allocated to demultiplex.
system(paste('chmod 755', file.path(opt$outputDir, 'core',  'demultiplex', 'run.sh')))

r <- system(file.path(opt$outputDir, 'core',  'demultiplex', 'run.sh'), wait = TRUE, intern = TRUE)
#r <- 0


if(r != 0){
  message('The demultiplex call from the demultiplex module failed.')
  q(save = 'no', status = 1, runLast = FALSE) 
}

updateLog('Demultiplex completed.')

# Read in demultiplex result.
reads <- readRDS(file.path(opt$outputDir, 'core', 'demultiplex', 'reads.rds'))

if(nrow(reads) > 0){
  updateLog(paste0(ppNum(n_distinct(reads$readID)), ' reads demultiplexed across ', ppNum(n_distinct(reads$uniqueSample)), ' sample replicates.'))
} else {
  updateLog('Error -- no reads were demultiplexed.')
  message('Error -- no reads were demultiplexed.')
  q(save = 'no', status = 1, runLast = FALSE)
}


jobStatus <- function(searchPattern = 'sites.done', outputFileName = 'jobTable.txt'){
  f <- list.files(file.path(opt$outputDir, 'core'), pattern = searchPattern, full = TRUE, recursive = TRUE)
  
  if(length(f) > 0){
    invisible(sapply(f, function(x){
      o <- unlist(strsplit(x, '/'))
      id <- o[length(o)-2]
      
      jobTable[jobTable$id == id,]$done       <<- TRUE
      jobTable[jobTable$id == id,]$active     <<- NA
      jobTable[jobTable$id == id,]$endTimeDsp <<- base::format(Sys.time(), "%m.%d.%Y %l:%M%P")
      jobTable[jobTable$id == id,]$endTime    <<- as.character(lubridate::now()) 
      jobTable[jobTable$id == id,]$duration   <<- paste0(as.integer(as.character(difftime(lubridate::now(), lubridate::ymd_hms(jobTable[jobTable$id == o[length(o)-2],]$startTime), units="mins"))), ' minutes')
      
      CPUs_used <<- CPUs_used - jobTable[jobTable$id == id,]$CPUs
      invisible(suppressWarnings(file.remove(f)))
    }))
  }
  
  # Redefine CPUs.
  o <- jobTable[jobTable$done == FALSE & jobTable$active == FALSE,]
  
  # End game CPU assignments.
  if(nrow(o) > 0 & nrow(o) <= 3){
    # a <- floor((1 / nrow(o)) * opt$core_CPUs) - 1
    a <- floor((1 / nrow(o)) * (opt$core_CPUs - CPUs_used))
    
    a[a < 1] <- 1
    
    a[a > opt$core_processMaxCPUs] <- opt$core_processMaxCPUs
    
    o$CPUs <- a
    
    jobTable[match(o$id, jobTable$id),]$CPUs <<- o$CPUs
  }
  
  dspJobTable <- jobTable
  dspJobTable$startTime <- dspJobTable$startTimeDsp;  dspJobTable$startTimeDsp <- NULL
  dspJobTable$endTime   <- dspJobTable$endTimeDsp;    dspJobTable$endTimeDsp   <- NULL
  
  tab <- pandoc.table.return(arrange(dspJobTable, desc(done), desc(active)), style = "simple", split.tables = Inf, plain.ascii = TRUE)
  
  write(c(date(), paste0('Available CPUs: ', (opt$core_CPUs - CPUs_used)), tab), file = file.path(opt$outputDir, 'core', outputFileName), append = FALSE)
}

# Build a jobs table for prepReads, alignReads, and buildFragments.
jobTable <- data.frame(table(reads$uniqueSample)) %>% 
            dplyr::rename(id = Var1) %>%
            dplyr::rename(reads = Freq)

a <- ceiling((jobTable$reads / max(jobTable$reads)) * (opt$core_CPUs * (opt$core_processMaxPercentCPU / 100)))

a[a < 1] <- 1
a[a > opt$core_processMaxCPUs] <- opt$core_processMaxCPUs

jobTable$CPUs         <- a
jobTable$active       <- FALSE
jobTable$startTime    <- NA
jobTable$endTime      <- NA
jobTable$startTimeDsp <- NA
jobTable$endTimeDsp   <- NA
jobTable$done         <- FALSE
jobTable$duration     <- NA

CPUs_used <- 0

updateLog('Starting replicate level jobs.')

# !!!!!!!!!
#jobTable$done <- TRUE

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
  o <- opt
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
  
  updateLog(paste0('Starting ',  tab$id, '.'))
  
  waitForMemory(stepDesc = 'Core module, replicate level jobs', minMem = opt$system_minMemThreshold, maxWaitSecs = opt$system_minMemWaitTime, sleepSecs = opt$system_minMemSleepTime)
  
  system(paste('chmod 755', file.path(opt$outputDir, 'core',  'replicate_analyses', tab$id, 'run.sh')))
  system(file.path(opt$outputDir, 'core', 'replicate_analyses', tab$id, 'run.sh'), wait = FALSE)
  
  CPUs_used <- CPUs_used + tab$CPUs
  
  jobTable[which(jobTable$id == tab$id),]$active <- TRUE
  jobTable[which(jobTable$id == tab$id),]$startTimeDsp <- base::format(Sys.time(), "%m.%d.%Y %l:%M%P")
  jobTable[which(jobTable$id == tab$id),]$startTime <- as.character(lubridate::now())
  
  jobStatus(searchPattern = 'fragments.done', outputFileName = 'replicateJobTable')
  Sys.sleep(5)
}

updateLog('Replicate level jobs completed.')

if(! opt$core_keepIntermediateFiles) invisible(file.remove(list.files(file.path(opt$outputDir, 'core', 'replicate_analyses'), recursive = TRUE, pattern = 'CORE_TMP', full.names = TRUE)))

updateLog('Bundling together replicate level fragment records.')

# Bundle together fragment output files.
f <- list.files(file.path(opt$outputDir, 'core'), pattern = 'fragments.rds', recursive = TRUE, full.names = TRUE)
frags <- bind_rows(lapply(f, readRDS))

# # Identify unique trial / patient combinations from returned fragments.
# u <- tidyr::separate(tibble(u = unique(frags$uniqueSample)), u, c('trial', 'subject', 'sample', 'replicate'), sep = '~') %>% select(-sample, -replicate) %>% distinct()
# 
# opt$core_maxCPUsPerProcess <- floor(opt$core_CPUs / nrow(u))
# if(opt$core_maxCPUsPerProcess == 0) opt$core_maxCPUsPerProcess <- 1
# 
# updateLog(paste0('Setting max. process CPU limit to ',opt$core_maxCPUsPerProcess, ' CPUs for subject level jobs.'))
# 
# opt$core_maxPercentCPUs <- floor((opt$core_maxCPUsPerProcess / opt$core_CPUs) * 100)

# Create a subject level to-do table.
jobTable <- tibble(u = unique(frags$uniqueSample))
jobTable$id <- sapply(jobTable$u, function(x) paste0(unlist(strsplit(x, '~'))[1:2], collapse = '~'))

# Create a subject-level splitting vector for fragments records.
frags <- left_join(frags, jobTable, by = c('uniqueSample' = 'u'))

# n is trial ~ subject
# Fragments are on the read level.
jobTable <- group_by(distinct(frags), id) %>% 
  summarise(reads = n_distinct(readID)) %>% 
  ungroup() %>% 
  arrange(desc(reads))

# Let the process with the greatest number of reads have 25% of all cores, scale other processes. 
a <- ceiling((jobTable$reads / max(jobTable$reads)) * (opt$core_CPUs * (opt$core_processMaxPercentCPU / 100)))

a[a < 1] <- 1
a[a > opt$core_processMaxCPUs] <- opt$core_processMaxCPUs

jobTable$CPUs         <- a
jobTable$active       <- FALSE
jobTable$startTime    <- NA
jobTable$endTime      <- NA
jobTable$startTimeDsp <- NA
jobTable$endTimeDsp   <- NA
jobTable$done         <- FALSE
jobTable$duration     <- NA

d <- tibble(uniqueSample = unique(frags$uniqueSample)) %>% tidyr::separate(uniqueSample, c('trial', 'subject', 'sample', 'replicate'), sep = '~', remove = FALSE)
d$trialSubject <- paste0(d$trial, '~', d$subject)
frags <- left_join(frags, distinct(dplyr::select(d, uniqueSample, trialSubject)), by = 'uniqueSample')

updateLog('Starting subject level jobs.')

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
  o <- opt
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
  
  updateLog(paste0('Starting ',  tab$id, '.'))
  
  waitForMemory(stepDesc = 'Core module, subject level jobs', minMem = opt$system_minMemThreshold, maxWaitSecs = opt$system_minMemWaitTime, sleepSecs = opt$system_minMemSleepTime)
  
  system(paste('chmod 755', file.path(opt$outputDir, 'core', 'subject_analyses', tab$id, 'run.sh')))
  system(file.path(opt$outputDir, 'core', 'subject_analyses', tab$id, 'run.sh'), wait = FALSE)
  
  CPUs_used <- CPUs_used + tab$CPUs
  
  jobTable[which(jobTable$id == tab$id),]$active <- TRUE
  jobTable[which(jobTable$id == tab$id),]$startTimeDsp <- base::format(Sys.time(), "%m.%d.%Y %l:%M%P")
  jobTable[which(jobTable$id == tab$id),]$startTime <- as.character(lubridate::now())
  
  jobStatus(searchPattern = 'sites.done', outputFileName = 'subjectJobTable')
  Sys.sleep(5)
}

updateLog('Subject level jobs completed.')

if(! opt$core_keepIntermediateFiles) invisible(file.remove(list.files(file.path(opt$outputDir, 'core', 'subject_analyses'), recursive = TRUE, pattern = 'CORE_TMP', full.names = TRUE)))

# Bundle together multi-hit cluster output files.
f <- list.files(file.path(opt$outputDir, 'core'), pattern = 'multiHitClusters.rds', recursive = TRUE, full.names = TRUE)

if(length(f) > 0){
  updateLog('Bundling multi-hit data objects.')
  multiHitClusters <- bind_rows(lapply(f, readRDS))
  saveRDS(multiHitClusters, file = file.path(opt$outputDir, 'core', 'multiHitClusters.rds'), compress = opt$compressDataFiles)
}

# Bundle together site output files.
f <- list.files(file.path(opt$outputDir, 'core'), pattern = 'sites.rds', recursive = TRUE, full.names = TRUE)

if(length(f) > 0){
  updateLog('Bundling integration site data objects.')
  sites <- bind_rows(lapply(f, readRDS))
  saveRDS(sites, file = file.path(opt$outputDir, 'core', 'sites.rds'), compress = opt$compressDataFiles)
  openxlsx::write.xlsx(sites, file.path(opt$outputDir, 'core', 'sites.xlsx'))
  readr::write_tsv(sites, file.path(opt$outputDir, 'core', 'sites.tsv.gz'))
} else {
  updateLog('Error -- no integration sites were recovered from any sample replicates.')
  message('Error -- no integration sites were recovered from any sample replicates.')
  q(save = 'no', status = 1, runLast = FALSE)
}

updateLog('Core module completed.')

q(save = 'no', status = 0, runLast = FALSE) 
