#!/usr/bin/Rscript

# AAVengeR/pullDatabaseFragments.R
# John K. Everett, Ph.D.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RMariaDB))

# Read in the configuration file and perform basic sanity checks.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) stop('Expected at least one command line argument')
source(file.path(yaml::read_yaml(args[1])$softwareDir, 'lib.R'))
opt <- startModule(args)

createOuputDir()
if(! dir.exists(file.path(opt$outputDir, opt$pullDatabaseFragments_outputDir))) dir.create(file.path(opt$outputDir, opt$pullDatabaseFragments_outputDir), showWarnings = FALSE)
if(! dir.exists(file.path(opt$outputDir, opt$pullDatabaseFragments_outputDir, 'tmp')))  dir.create(file.path(opt$outputDir, opt$pullDatabaseFragments_outputDir, 'tmp'), showWarnings = FALSE)

opt$defaultLogFile <- file.path(opt$outputDir, opt$pullDatabaseFragments_outputDir, 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)
write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)

quitOnErorr <- function(msg){
  updateLog(msg)
  updateLog(paste0('See log for more details: ', opt$defaultLogFile))
  updateMasterLog()
  q(save = 'no', status = 1, runLast = FALSE) 
}

dbConn <- createDBconnection()

frags <- bind_rows(lapply(unlist(strsplit(opt$pullDatabaseFragments_trialSubjectSamples, '\\s*\\|\\s*')), function(x){
           x <- unlist(strsplit(x, '\\s*,\\s*'))
           query <- paste0("select data from fragments where trial = '", x[1], "' and subject = '", x[2], "' and sample = '", x[3], "' and refGenome = '", x[4], "'")
           updateLog(paste0('Submitting query: ', query))
           o <- reconstructDBtable(dbGetQuery(dbConn, query), tmpDirPath = file.path(opt$outputDir, opt$pullDatabaseFragments_outputDir, 'tmp'))
           updateLog(paste0(ppNum(nrow(o)), ' rows returned from the database.'))
           o
         }))

dbDisconnect(dbConn)

saveRDS(frags, file.path(opt$outputDir, opt$pullDatabaseFragments_outputDir, 'fragments.rds'))

updateLog('pullDatabaseFragments completed.')

q(save = 'no', status = 0, runLast = FALSE) 
