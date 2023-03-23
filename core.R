library(dplyr)
library(parallel)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)
opt$runDir <- getwd()
source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$runDir, file.path(opt$outputDir)))
dir.create(file.path(opt$runDir, opt$outputDir, 'core'))
write(paste0(date(), ' - start'), file = file.path(opt$runDir, opt$outputDir, 'core', 'log'))

opt$demultiplex_CPUs <- opt$core_CPUs

# Gather select sections from the configuration file.
o <- opt[grepl('^Rscript|^softwareDir|^outputDir|^runDir|^databaseGroup|^core|^demultiplex', names(opt))]
o$outputDir <- file.path(opt$runDir, opt$outputDir, 'core')

# Run demultiplex.
yaml::write_yaml(o, file.path(opt$runDir, opt$outputDir, 'core',  'demultiplex_config.yml'))
system(paste(opt$Rscript, file.path(opt$softwareDir, 'demultiplex.R'), file.path(opt$runDir, opt$outputDir, 'core', 'demultiplex_config.yml')))

reads <- readRDS(file.path(opt$runDir, opt$outputDir, 'core', 'demultiplex', 'reads.rds'))

jobStatus <- function(searchPattern = 'sites.done', sleepSeconds = 5){
  f <- list.files(file.path(opt$runDir, opt$outputDir, 'core'), pattern = searchPattern, full = TRUE, recursive = TRUE)
  
  if(length(f) > 0){
    # For each found done file, find the matching entry in the jobTable.
    i <- jobTable$id %in% sapply(f, function(x){
                            o <- unlist(strsplit(x, '/'))
                            o[length(o)-2]
                          }, USE.NAMES = FALSE)
    
    if(any(i)) jobTable[i,]$done <<- TRUE
    nActiveJobs <<- nActiveJobs - length(f)
    nJobsRemaining <<- nJobsRemaining - length(f)
    if(nJobsRemaining < maxSimultaneousJobs) maxSimultaneousJobs <<- nJobsRemaining   # This will increase the number of CPUs allowed per job
    invisible(file.remove(f))
  }
  
  print(jobTable)
  message('waiting - nActiveJobs: ', nActiveJobs, ', maxSimultaneousJobs: ', maxSimultaneousJobs, ', nJobsRemaining: ', nJobsRemaining - nActiveJobs)
  Sys.sleep(sleepSeconds)
}

jobTable <- data.frame(table(reads$uniqueSample)) %>% dplyr::rename(id = Var1) %>% arrange(Freq) %>% select(-Freq)
jobTable$done <- FALSE

maxSimultaneousJobs <- ceiling(opt$core_CPUs / opt$core_minCPUs)
nActiveJobs <- 0
nJobsRemaining <- nrow(jobTable)

# Run prepReads, alignReads, and buildFragments.
for(j in jobTable$id){
  # Gather select sections from the configuration file and set CPU values.
  o <- opt[grepl('^Rscript|^softwareDir|^outputDir|^runDir|^databaseGroup|^core|^prepReads|^alignReads|^buildFragments', names(opt))]
  o$prepReads_CPUs <- as.integer(ceiling(opt$core_CPUs / maxSimultaneousJobs))
  o$alignReads_CPUs <- as.integer(ceiling(opt$core_CPUs / maxSimultaneousJobs))
  o$buildFragments_CPUs <- as.integer(ceiling(opt$core_CPUs / maxSimultaneousJobs))
  o$outputDir <- file.path(opt$runDir, opt$outputDir, 'core', j)
  
  # Instruct the pipeline to create a buildFragments/fragments.done file for jobs that failed.
  o$core_createFauxFragDoneFiles <- TRUE
  
  # Create the output directory for this replicate and save a copy of the replicate reads within it.
  dir.create(o$outputDir)
  saveRDS(subset(reads, uniqueSample == j), file.path(o$outputDir, paste0(j, '.rds')))
  
  # Update configuration to point to the subset of replicate reads.
  o$prepReads_readsTable <- paste0(j, '.rds')
  
  # Define modules to run in replicate level configuration file and write out file.
  o$modules <- list()
  o[['modules']] <- c('prepReads', 'alignReads', 'buildFragments')
  yaml::write_yaml(o, file.path(opt$runDir, opt$outputDir, 'core',  paste0(j,'.yml')))
  
  # Create AAVengeR launching script since system(x, wait = TRUE) does not wait for commands w/ arguments.
  write(c('#!/usr/bin/sh', 
          paste(opt$Rscript, file.path(opt$softwareDir, 'aavenger.R'), file.path(opt$runDir, opt$outputDir, 'core',  paste0(j,'.yml')))), 
        file = file.path(opt$runDir, opt$outputDir, 'core',  paste0(j,'.sh')))
  
  system(paste('chmod 755', file.path(opt$runDir, opt$outputDir, 'core',  paste0(j,'.sh'))))
  system(file.path(opt$runDir, opt$outputDir, 'core',  paste0(j,'.sh')), wait = FALSE, show.output.on.console = FALSE)
  nActiveJobs <- nActiveJobs + 1
  
  # Wait before starting more jobs.
  while(nActiveJobs >= maxSimultaneousJobs){
    jobStatus(searchPattern = 'fragments.done')
    if(all(jobTable$done == TRUE)) break
  }  
}

# Catch remaining jobs.
while(TRUE){
  jobStatus(searchPattern = 'fragments.done')
  if(all(jobTable$done == TRUE)) break
  sleep(5)
}

# Bundle together fragment output files.
f <- list.files(file.path(opt$outputDir, 'core'), pattern = 'fragments.rds', recursive = TRUE, full.names = TRUE)
frags <- bind_rows(lapply(f, readRDS))
saveRDS(frags, file = file.path(opt$outputDir, 'core', 'fragments.rds'))

# Create a subject level to-do table.
jobTable <- tibble(u = unique(frags$uniqueSample))
jobTable$id <- sapply(jobTable$u, function(x) paste0(unlist(strsplit(x, '~'))[1:2], collapse = '~'))

# Create a subject-level splitting vector for fragments records.
frags <- left_join(frags, jobTable, by = c('uniqueSample' = 'u'))

jobTable <- group_by(distinct(frags), id) %>% summarise(n = n()) %>% ungroup() %>% arrange(n) %>% mutate(done = FALSE) %>% select(-n)

maxSimultaneousJobs <- ceiling(opt$core_CPUs / opt$core_minCPUs)
nActiveJobs <- 0
nJobsRemaining <- nrow(jobTable)
if(nJobsRemaining < maxSimultaneousJobs) maxSimultaneousJobs <- nJobsRemaining  

# Run buildStdFragments and buildSites.
for(j in jobTable$id){
  # Gather select sections from the configuration file and set CPU values.
  o <- opt[grepl('^Rscript|^softwareDir|^outputDir|^runDir|^databaseGroup|^core|^buildStdFragments|^buildSites', names(opt))]
  o$buildStdFragments_CPUs <- as.integer(ceiling(opt$core_CPUs / maxSimultaneousJobs))
  o$buildSites_CPUs <- as.integer(ceiling(opt$core_CPUs / maxSimultaneousJobs))
  o$outputDir <- file.path(opt$runDir, opt$outputDir, 'core', j)
  
  # Instruct the pipeline to create a buildSites/sites.done file for jobs that failed.
  o$core_createFauxSiteDoneFiles <- TRUE

  # Create the output directory for this replicate and save a copy of the replicate reads within it.
  dir.create(o$outputDir)
  saveRDS(subset(frags, id == j), file.path(o$outputDir, paste0('stdFragments.rds')))
 
  # Update configuration to point to the subset of replicate reads.
  o$buildStdFragments_inputFile <- 'stdFragments.rds'
  
  # Define modules to run in replicate level configuration file and write out file.
  o$modules <- list()
  o[['modules']] <- c('buildStdFragments', 'buildSites')
  yaml::write_yaml(o, file.path(opt$runDir, opt$outputDir, 'core',  paste0(j,'.yml')))
  
  # Create AAVengeR launching script since system(x, wait = TRUE) does not wait for commands w/ arguments.
  write(c('#!/usr/bin/sh', 
          paste(opt$Rscript, file.path(opt$softwareDir, 'aavenger.R'), 
          file.path(opt$runDir, opt$outputDir, 'core',  paste0(j,'.yml')))), 
        file = file.path(opt$runDir, opt$outputDir, 'core',  paste0(j,'.sh')))
  
  system(paste('chmod 755', file.path(opt$runDir, opt$outputDir, 'core',  paste0(j,'.sh'))))
  system(file.path(opt$runDir, opt$outputDir, 'core',  paste0(j, '.sh')), wait = FALSE, show.output.on.console = FALSE)
  nActiveJobs <- nActiveJobs + 1
  
  while(nActiveJobs >= maxSimultaneousJobs){
    jobStatus(searchPattern = 'sites.done')
    if(all(jobTable$done == TRUE)) break
  }  
}

# Catch remaining jobs.
while(TRUE){
  jobStatus(searchPattern = 'sites.done')
  if(all(jobTable$done == TRUE)) break
  sleep(5)
}

# Bundle together site output files.
f <- list.files(file.path(opt$outputDir, 'core'), pattern = 'sites.rds', recursive = TRUE, full.names = TRUE)
sites <- bind_rows(lapply(f, readRDS))
saveRDS(sites, file = file.path(opt$outputDir, 'core', 'sites.rds'))
openxlsx::write.xlsx(sites, file.path(opt$outputDir, 'core', 'sites.xlsx'))
readr::write_tsv(sites, file.path(opt$outputDir, 'core', 'sites.tsv.gz'))

write(paste0(date(), ' - end'), file = file.path(opt$runDir, opt$outputDir, 'core', 'log'), append = TRUE)
