#!/usr/bin/Rscript
options(scipen = 999, useFancyQuotes = FALSE)

# AAVengeR/core.R
# John K. Everett, Ph.D.
# 
# This script runs AAVengeR's six core modules:
# demultiplex, prepReads, alignReads, buildFragments, buildStdFragments, and buildSites.
# It is advantageous to use this module rather than calling the core modules serially
# for large data sets since the core module dynamically sets CPU thresholds for each 
# sample replicate according to the number of associated reads.

for (p in c('dplyr', 'pander')) suppressPackageStartupMessages(library(p, character.only = TRUE))

# Read in the configuration file and perform basic sanity checks.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) stop('Expected at least one command line argument')
source(file.path(yaml::read_yaml(args[1])$softwareDir, 'lib', 'main.R'))
opt <- startModule(args)


# Create directories and start logging.
createOuputDir()
dir.create(file.path(opt$outputDir, 'core'), showWarnings = FALSE)
dir.create(file.path(opt$outputDir, 'core', 'demultiplex'), showWarnings = FALSE)
dir.create(file.path(opt$outputDir, 'core', 'replicate_analyses'), showWarnings = FALSE)
dir.create(file.path(opt$outputDir, 'core', 'subject_analyses'), showWarnings = FALSE) 

# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, 'core', 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)
write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)

quitOnErorr <- function(msg){
  updateLog(msg)
  updateLog(paste0('See log for more details: ', opt$defaultLogFile))
  closeAllConnections()
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
        paste(opt$Rscript, file.path(opt$softwareDir, 'modules', 'demultiplex.R'), file.path(opt$outputDir, 'core',  'demultiplex', 'config.yml')),
        'echo $?'), 
        file = file.path(opt$outputDir, 'core',  'demultiplex', 'run.sh'))

updateLog('Running the demultiplex module.')

# All core CPUs allocated to demultiplex.
system(paste('chmod 755', file.path(opt$outputDir, 'core',  'demultiplex', 'run.sh')))


# Run demultiplex module.
if(opt$core_skipDemultiplexing == FALSE){
  r <- system(file.path(opt$outputDir, 'core',  'demultiplex', 'run.sh'), wait = TRUE, intern = TRUE)
} else {
  r <- 0
}

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
  closeAllConnections()
  q(save = 'no', status = 1, runLast = FALSE)
}


distribute_integer <- function(total, n) {
  base_value <- total %/% n
  remainder <- total %% n
  
  distribution <- rep(base_value, n)
  if (remainder > 0) {
    distribution[1:remainder] <- distribution[1:remainder] + 1
  }
  
  return(distribution)
}


allocate_CPUs <- function(totalCPUs, jobReads, minCPUsPerJob = 1, maxCPUsPerJob = totalCPUs, maxReadsPerCPU = 5000, readsDownStep = 1){
  o <- round(jobReads / maxReadsPerCPU)             # Set an initial number of CPUs for each job.
  o <- ifelse(o > opt$core_CPUs, opt$core_CPUs, o)  # Make sure no job has more than the total number of available CPUs.
  o <- ifelse(o < minCPUsPerJob, minCPUsPerJob, o)  # Make sure each job has at least the min. number of CPUs.
  
  o_min <- o
    
  # Allocate more CPUs if available.
  if(sum(o) < totalCPUs){
    while(sum(o) < totalCPUs){
        maxReadsPerCPU <- maxReadsPerCPU - readsDownStep
        o <- ifelse(round(jobReads / maxReadsPerCPU) == 0, 1, round(jobReads / maxReadsPerCPU))
    }
  
    # Handle rounding issues that occasionally lead to high CPU allocations.
    maxReadsPerCPU <- maxReadsPerCPU + readsDownStep
    o <- ifelse(round(jobReads / maxReadsPerCPU) == 0, 1, round(jobReads / maxReadsPerCPU))
  }
  
  o <- ifelse(is.infinite(o), minCPUsPerJob, o)
  
  ifelse(o < o_min, o_min, o)
}


jobStatus <- function(searchPattern = 'sites.done', outputFileName = 'jobTable.txt'){
  f <- list.files(file.path(opt$outputDir, 'core'), pattern = searchPattern, full = TRUE, recursive = TRUE)
  
  # Mark completed jobs, update CPUs_used.
  if(length(f) > 0){
    invisible(sapply(f, function(x){
      o <- unlist(strsplit(x, '/'))
      id <- o[length(o)-2]
      
      jobTable[jobTable$id == id,]$done       <<- TRUE
      jobTable[jobTable$id == id,]$active     <<- NA
      jobTable[jobTable$id == id,]$endTimeDsp <<- base::format(Sys.time(), "%m.%d.%Y %l:%M%P")
      jobTable[jobTable$id == id,]$endTime    <<- as.character(lubridate::now()) 
      jobTable[jobTable$id == id,]$duration   <<- paste0(as.integer(as.character(difftime(lubridate::now(), lubridate::ymd_hms(jobTable[jobTable$id == o[length(o)-2],]$startTime), units="mins"))), ' minutes')
      
      # Return the CPUs used for this completed job.
      CPUs_used <<- CPUs_used - jobTable[jobTable$id == id,]$CPUs
      invisible(suppressWarnings(file.remove(f)))
    }))
  }
  
  free_CPUs <- opt$core_CPUs - CPUs_used
  
  # Determine which jobs still need to be completed.
  o <- jobTable[jobTable$done == FALSE & jobTable$active == FALSE,]
   
  if(nrow(o) > 0 & free_CPUs > 0){
    o$CPUs <- allocate_CPUs(free_CPUs, o$reads, maxReadsPerCPU = opt$core_readsPerCPU, minCPUsPerJob = opt$core_minCPUsPerJob, readsDownStep = opt$core_readsDownSte)
    jobTable[match(o$id, jobTable$id),]$CPUs <<- o$CPUs
   }
  
  jobTable <<- arrange(jobTable, desc(reads))
  
  dspJobTable <- jobTable
  dspJobTable$startTime <- dspJobTable$startTimeDsp;  dspJobTable$startTimeDsp <- NULL
  dspJobTable$endTime   <- dspJobTable$endTimeDsp;    dspJobTable$endTimeDsp   <- NULL
  
  tab <- pandoc.table.return(arrange(dspJobTable, desc(done), desc(active)), style = "simple", split.tables = Inf, plain.ascii = TRUE)
  
  freeMem <- system("free -h | awk '/^Mem:/ {print $4}'", intern =  TRUE)
 
  cmd <- "free -m | awk '/^Mem:/ {printf \"%.0f%% of total memory used.\\n\", ($3/$2)*100}'"
  percentMemUsed <- system(cmd, intern = TRUE)
  
  write(c(date(), paste0('Available CPUs: ', free_CPUs, '\nAvailable memory: ', freeMem, '\n', percentMemUsed), tab), file = file.path(opt$outputDir, 'core', outputFileName), append = FALSE)
  updateMasterLog()
}

# Build an initial jobs table for prepReads, alignReads, and buildFragments.
jobTable <- data.frame(table(reads$uniqueSample)) %>% 
            dplyr::rename(id = Var1) %>%
            dplyr::rename(reads = Freq)

# Assign a job an initial number of cores based on the number of its reads.
jobTable$CPUs         <- allocate_CPUs(opt$core_CPUs, jobTable$reads, maxReadsPerCPU = opt$core_readsPerCPU, minCPUsPerJob = opt$core_minCPUsPerJob, readsDownStep = opt$core_readsDownStep)
jobTable$active       <- FALSE
jobTable$startTime    <- NA
jobTable$endTime      <- NA
jobTable$startTimeDsp <- NA
jobTable$endTimeDsp   <- NA
jobTable$done         <- FALSE
jobTable$duration     <- NA

jobTable <- arrange(jobTable, desc(reads))

CPUs_used <- 0

updateLog('Starting replicate level jobs.')

if(opt$core_skipReplicateLevelJobs) jobTable$done <- TRUE

# Run prepReads, alignReads, and buildFragments.
while(! all(jobTable$done == TRUE)){
  
  # Determine the number of free CPUs.
  CPUs_available <- opt$core_CPUs - CPUs_used
  
  # Select a single job requiring the most CPUs that can be done now.
  tab <- jobTable[jobTable$CPUs <= CPUs_available & jobTable$active == FALSE & jobTable$done == FALSE,] %>% dplyr::slice_max(CPUs, with_ties = FALSE)
  
  # If no jobs can be started now, search for completed jobs then loop.
  if(nrow(tab) == 0){
    Sys.sleep(5)
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
  dir.create(o$outputDir, showWarnings = FALSE)
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
  
  # Record the number of CPUs used.
  CPUs_used <- CPUs_used + tab$CPUs
  
  # Update Job table.
  jobTable[which(jobTable$id == tab$id),]$active <- TRUE
  jobTable[which(jobTable$id == tab$id),]$startTimeDsp <- base::format(Sys.time(), "%m.%d.%Y %l:%M%P")
  jobTable[which(jobTable$id == tab$id),]$startTime <- as.character(lubridate::now())
  
  # Search for completed jobs.
  jobStatus(searchPattern = 'fragments.done', outputFileName = 'replicateJobTable')
  Sys.sleep(5)
}

updateLog('Replicate level jobs completed.')
updateMasterLog()

if(! opt$core_keepIntermediateFiles) invisible(file.remove(list.files(file.path(opt$outputDir, 'core', 'replicate_analyses'), recursive = TRUE, pattern = 'CORE_TMP', full.names = TRUE)))

updateLog('Bundling together replicate level fragment records.')

# Bundle together fragment output files.
f <- list.files(file.path(opt$outputDir, 'core'), pattern = 'fragments.rds', recursive = TRUE, full.names = TRUE)

if(length(f) == 0){
  updateLog('Error -- no fragments were generated for any of the sample replicates.')
  message('Error -- no fragments were generated for any of the sample replicates.')
  closeAllConnections()
  q(save = 'no', status = 1, runLast = FALSE)
}

frags <- bind_rows(lapply(f, readRDS))

# Identify unique trial / patient combinations from returned fragments.
u <- tidyr::separate(tibble(u = unique(frags$uniqueSample)), u, c('trial', 'subject', 'sample', 'replicate'), sep = '~') %>% select(-sample, -replicate) %>% distinct()

opt$core_maxCPUsPerProcess <- floor(opt$core_CPUs / nrow(u))
if(opt$core_maxCPUsPerProcess == 0) opt$core_maxCPUsPerProcess <- 1

updateLog(paste0('Setting max. process CPU limit to ',opt$core_maxCPUsPerProcess, ' CPUs for subject level jobs.'))

opt$core_maxPercentCPUs <- floor((opt$core_maxCPUsPerProcess / opt$core_CPUs) * 100)

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

jobTable$CPUs         <- allocate_CPUs(opt$core_CPUs, jobTable$reads, maxReadsPerCPU = opt$core_readsPerCPU, minCPUsPerJob = opt$core_minCPUsPerJob, readsDownStep = opt$core_readsDownSte)
jobTable$active       <- FALSE
jobTable$startTime    <- NA
jobTable$endTime      <- NA
jobTable$startTimeDsp <- NA
jobTable$endTimeDsp   <- NA
jobTable$done         <- FALSE
jobTable$duration     <- NA

jobTable <- arrange(jobTable, desc(reads))

d <- tibble(uniqueSample = unique(frags$uniqueSample)) %>% tidyr::separate(uniqueSample, c('trial', 'subject', 'sample', 'replicate'), sep = '~', remove = FALSE)
d$trialSubject <- paste0(d$trial, '~', d$subject)
frags <- left_join(frags, distinct(dplyr::select(d, uniqueSample, trialSubject)), by = 'uniqueSample')

updateLog('Starting subject level jobs.')
updateMasterLog()

CPUs_used <- 0

while(! all(jobTable$done == TRUE)){
  
  # Determine how many CPUs are available to start jobs.
  CPUs_available <- opt$core_CPUs - CPUs_used
  
  # Select the job requiring the most CPUs that can be started now.
  tab <- jobTable[jobTable$CPUs <= CPUs_available & jobTable$active == FALSE & jobTable$done == FALSE,] %>% dplyr::slice_max(CPUs, with_ties = FALSE)
  
  # If no job can be started, look for completed jobs the loop again.
  if(nrow(tab) == 0){
    Sys.sleep(5)
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
  dir.create(o$outputDir, showWarnings = FALSE)
  
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
updateMasterLog()

if(! opt$core_keepIntermediateFiles) invisible(file.remove(list.files(file.path(opt$outputDir, 'core', 'subject_analyses'), recursive = TRUE, pattern = 'CORE_TMP', full.names = TRUE)))

# Bundle together multi-hit cluster output files.
f <- list.files(file.path(opt$outputDir, 'core'), pattern = 'multiHitClusters.rds', recursive = TRUE, full.names = TRUE)

if(length(f) > 0){
  updateLog('Bundling multi-hit data objects.')
  multiHitClusters <- bind_rows(lapply(f, readRDS))
  saveRDS(multiHitClusters, file = file.path(opt$outputDir, 'core', 'multiHitClusters.rds'), compress = opt$compressDataFiles)
  readr::write_tsv(multiHitClusters,  file.path(opt$outputDir, 'core', 'multiHitClusters.tsv'))
}


# Bundle together anchor read clusters files.
f <- list.files(file.path(opt$outputDir, 'core'), pattern = 'anchorReadClusterDecisionTable.rds', recursive = TRUE, full.names = TRUE)

if(length(f) > 0){
  updateLog('Bundling anchorReadClusterTable data objects.')
  anchorReadClusterTable <- bind_rows(lapply(f, readRDS))
  saveRDS(anchorReadClusterTable, file = file.path(opt$outputDir, 'core', 'anchorReadClusterDecisionTable.rds'), compress = opt$compressDataFiles)
  readr::write_tsv(anchorReadClusterTable, file = file.path(opt$outputDir, 'core', 'anchorReadClusterDecisionTable.tsv'))
}


# Collect and write out read attrition values.
updateLog('Compiling attrition logs.')
atl <- bind_rows(lapply(list.files(opt$outputDir, pattern = 'attritionLog.tsv', recursive = TRUE, full.names = TRUE), function(x){
  d <- readr::read_tsv(x, show_col_type = FALSE)
  #browser()
  d$sampleRep <- rev(unlist(strsplit(x, '/')))[3]
  d
}))

if(nrow(atl) > 0){
  updateLog(paste0(nrow(atl), ' rows found in compiled attrition logs.'))
  atlOutput <- vector()
  
  # Protection against duplicate rows from re-running analyses.
  i <- duplicated(paste(atl$label, atl$sampleRep))
  if(length(i) > 0) atl <- atl[!i,]
  
  atl <- tidyr::pivot_wider(distinct(atl), names_from = label, values_from = value)

  cols <- c("sampleRep", "PRD1", "PRD2", "PRD3", "PRD4", "PRD5", "ALR1", "ALR2", "ALR3", "ALR4", "ALR5", "BFR1", "BSF1")
  cols <- cols[cols %in% names(atl)]
  atl <- atl[, cols]
  
  updateLog('Writing out tsv formatted attrition log.')
  readr::write_tsv(atl, file = file.path(opt$outputDir, 'core', 'readAttrition.tsv'))
  
  updateLog('Compiling readAttritionReport.txt.')
  
  invisible(lapply(split(atl, atl$sampleRep), function(x){
    sampleRep <- x$sampleRep[1]
    x <- tidyr::pivot_longer(x, cols[cols != 'sampleRep'])[,2:3]
    names(x) <- c('label', 'value')
    
    reads <- as.integer(x[1,]$value)    # Get top row value for divisor.
    if(! is.integer(reads)) reads <- NA # Safe for division.
    x$value <- as.integer(x$value)      # Ensure integer type.
    x$value <- x$value / reads
   
    x <- x[! is.na(x$value),]
    atlOutput <<- c(atlOutput, sampleRep, paste0(ppNum(reads), ' demultiplexed reads'), asciiPercentBarChart(x), '\n')
  }))
  
  write(c(readLines(file.path(opt$softwareDir, 'figures', 'readAttritionDesc.txt')), atlOutput), file.path(opt$outputDir, 'core', 'readAttritionReport.txt'))
}


# Bundle together site output files.
f <- list.files(file.path(opt$outputDir, 'core'), pattern = 'sites.rds', recursive = TRUE, full.names = TRUE)

if(length(f) > 0){
  updateLog('Bundling integration site data objects.')
  sites <- bind_rows(lapply(f, readRDS))
  saveRDS(sites, file = file.path(opt$outputDir, 'core', 'sites.rds'), compress = opt$compressDataFiles)
  openxlsx::write.xlsx(sites, file.path(opt$outputDir, 'core', 'sites.xlsx'))
  readr::write_tsv(sites, file.path(opt$outputDir, 'core', 'sites.tsv'))
} else {
  updateLog('Error -- no integration sites were recovered from any sample replicates.')
  message('Error -- no integration sites were recovered from any sample replicates.')
  closeAllConnections()
  q(save = 'no', status = 1, runLast = FALSE)
}

updateLog('Core module completed.')
updateMasterLog()
closeAllConnections()

q(save = 'no', status = 0, runLast = FALSE) 
