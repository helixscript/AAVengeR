# AAVengeR/aavenger.R
# John K. Everett, Ph.D.
# 
# This is the launch script for the AAVengeR pipeline.
# Its primary function is to check for configuration errors and execute
# the module list found in the configuration file.

suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(pander))

configFile = commandArgs(trailingOnly=TRUE)

if(nchar(configFile) == 0){
  stop('Error -- the config file paramater was not provided.')
}

if(! file.exists(configFile)){
  stop('Error -- the config file could not be found.')
}

opt <- read_yaml(configFile)
source(file.path(opt$softwareDir, 'lib.R'))

createOuputDir()

dir.create(file.path(opt$outputDir, 'src'))

invisible(sapply(list.files(opt$softwareDir, pattern = '*.R', full.names = TRUE), file.copy, to = file.path(opt$outputDir, 'src')))

# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)
write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)

# Check for third party software.
checkSoftware()

invisible(file.copy(opt$sampleConfigFile, file.path(opt$outputDir, 'src', 'sampleData.tsv')))
invisible(file.copy(configFile, file.path(opt$outputDir, 'src', 'config.yml')))

updateLog("Starting AAVengeR's module chain.")
updateLog("Each module will create its own output folder in the main output folder.")
updateLog("The log files for each module are stored within their output folders.")

# Execute each module defined in opt$modules.
invisible(lapply(opt$modules, function(m){
  updateLog(paste0('Starting ', m, '.'))
  
  r <- system(paste(opt$Rscript, file.path(opt$softwareDir, paste0(m, '.R')), configFile))
  if(r != 0){
    updateLog(paste0('Module ', m, ' failed. Please see logs for details.'))
    q(save = 'no', status = 1, runLast = FALSE)
  }
}))

updateLog("Done.")
q(save = 'no', status = 0, runLast = FALSE)