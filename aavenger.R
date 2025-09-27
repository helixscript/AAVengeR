#!/usr/bin/Rscript

# AAVengeR/aavenger.R
# John K. Everett, Ph.D.
# 
# This is the launch script for the AAVengeR pipeline.
# Its primary function is to check for configuration errors and execute
# the module list found in the configuration file.

for (p in c('yaml', 'dplyr', 'lubridate', 'pander', 'optparse')) suppressPackageStartupMessages(library(p, character.only = TRUE))

args <- commandArgs(trailingOnly=TRUE)

knownParameters <- c('list_available_genomes', 'list_installed_genomes', 'install_genome', 'remove_genome')

if(length(args) == 0){
  message('usage: ')
  message('  aavenger.R <config file>')
  message('  aavenger.R list_available_genomes')
  message('  aavenger.R list_installed_genomes')
  message('  aavenger.R install_genome <genome id>')
  message('  aavenger.R remove_genome <genome id>')
  quit(save = "no", status = 0, runLast = FALSE)
}

parseParameters <- function(){
  if(args[1] == 'list_available_genomes'){
    o <- yaml::read_yaml('data/defaults.yml')
    f <- system(paste0('wget -qO- ', o$remoteDataURL, '/genomePackages/fileList'), intern = TRUE)
    f <- f[2:length(f)]
    f <- f[! grepl('fileList', f)]
    sizes <-  unlist(lapply(strsplit(f, '\\s+'), '[[', 5))
    sizes <- sub('G', ' GB', sizes)
    sizes <- sub('M', ' MB', sizes)
    genomes <- unlist(lapply(strsplit(f, '\\s+'), '[[', 9))
    genomes <- sub('\\.tar$', '', genomes)
    message('Available genomes:', pandoc.table.return(tibble(refGenome = genomes, dataSize = sizes), style = "simple", split.tables = Inf, plain.ascii = TRUE))
    closeAllConnections()
    quit(save = "no", status = 0, runLast = FALSE)
  }
  
  if(args[1] == 'list_installed_genomes'){
    f <- list.files('data/referenceGenomes/blat', pattern = "*.2bit", full.names = TRUE)
    
    if(length(f) == 0){
      message('No installed genomes found.\n')
      closeAllConnections()
      quit(save = "no", status = 0, runLast = FALSE)
    }
    
    g <- unique(sub('\\.2bit', '', stringr::str_extract(f, '[a-zA-Z\\d]+.2bit')))
    
    tab <- bind_rows(lapply(g, function(x){
      f <- list.files('data', recursive = TRUE, pattern = paste0('^', x, '.'), full.names = TRUE)
      tibble(refGenome = x, dataFiles = length(f), spaceUsed = paste0(round(sum(file.size(f)) / 1024 / 1024 / 1024, 2), ' GB'))
    }))
    
    message(pandoc.table.return(tab, style = "simple", split.tables = Inf, plain.ascii = TRUE))
    closeAllConnections()
    quit(save = "no", status = 0, runLast = FALSE)
  }
  
  if(args[1] == 'install_genome'){
    if(length(args) == 1) stop('The install_genome parameter requires a genome identifier.')
    
    o <- yaml::read_yaml('data/defaults.yml')
    f <- system(paste0('wget -qO- ', o$remoteDataURL, '/genomePackages/fileList'), intern = TRUE)
    f <- f[2:length(f)]
    f <- f[! grepl('fileList', f)]
    genomes <- unlist(lapply(strsplit(f, '\\s+'), '[[', 9))
    genomes <- sub('\\.tar$', '', genomes)
    
    if(! args[2] %in% genomes){
      message("Error - ", args[2], " is not in the list of available genomes.\n")
      closeAllConnections()
      quit(save = "no", status = 0, runLast = FALSE)
    }
    
    ff <- list.files('data/referenceGenomes/blat', pattern = "*.2bit", full.names = TRUE)
    gg <- unique(sub('\\.2bit', '', stringr::str_extract(ff, '[a-zA-Z\\d]+.2bit')))
    
    if(args[2] %in% gg){
      message('Error - ', args[2], ' is already installed.')
      closeAllConnections()
      quit(save = "no", status = 0, runLast = FALSE)
    }
    
    tmpDir <- paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')
    dir.create(tmpDir, showWarnings = FALSE)
    
    message('Downloading genome data.')
    system(paste0('wget -q -O ', tmpDir, '/archive.tar ', o$remoteDataURL, '/genomePackages/', args[2], '.tar'))
    
    message('Extracting data.')
    system(paste0('tar xvf ', tmpDir, '/archive.tar -C ', tmpDir))
    
    message('Installing data.')
    system(paste0('rsync -a ', tmpDir, '/data/ data/'))
    unlink(tmpDir, recursive = TRUE)
    
    message('done.\n')
    closeAllConnections()
    quit(save = "no", status = 0, runLast = FALSE)
  }
  
  if(args[1] == 'remove_genome'){
    if(length(args) == 1) stop('The remove_genome parameter requires a genome identifier.')
    
    f <- list.files('data', pattern = paste0('^', args[2], '\\.'), recursive = TRUE, full.names = TRUE)
    
    if(length(f) == 0){
      message('No genome files were found matching ', args[2], '.\n')
      closeAllConnections()
      quit(save = "no", status = 0, runLast = FALSE)
    }
    
    print(f)
    message("Remove these files? (y/n)?")
    line <- readLines("stdin", n=1)
    
    if(line == 'y'){
      invisible(file.remove(f)) 
    } 
    
    message(args[2], ' removed.')
    closeAllConnections()
    quit(save = "no", status = 0, runLast = FALSE)
  }
}

# If a flagged parameter was provided, provide the requested action then exit.
if(args[1] %in% knownParameters) parseParameters()
 
# Exit if config file does not exist.
if(! file.exists(args[1])) stop(paste0("Config file '", args[1], "' does not exist."))

# Read in the configuration file and perform basic sanity checks.
set.seed(1)
source(file.path(yaml::read_yaml(args[1])$softwareDir, 'lib', 'main.R'))
opt <- startModule(args)

# Show AAVengeR banner.
if(! opt$calledFromCore) invisible(sapply(readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt')), message))

# Test database connection if databasing is turned on.
if('database_configGroup' %in% names(opt)){
  if(opt$database_configGroup != 'none'){
    suppressPackageStartupMessages(library(RMariaDB))

    if(! file.exists(opt$database_configFile)){
      stop('Error - database_configGroup is not set to none (databasing turned on) and database_configFile does not point to a file.')
    }
    
    dbConn <- tryCatch({
      dbConnect(RMariaDB::MariaDB(), group = opt$database_configGroup, default.file = opt$database_configFile)
    },
    error=function(cond) {
      stop(paste0('Error - databasing was requested but AAVengeR could not connect to the database defined in the ', 
                  opt$database_configGroup, ' block of ', opt$database_configFile, '.'))
    })

    dbDisconnect(dbConn)
  }
}

# Create main output folder and make a copy of the source code and configurations.
createOuputDir()

# Check for third party software.
checkSoftware()

# Make a copy of all R and config files if aavenger.R is not called from the core module.
if(! opt$calledFromCore){
  if(! dir.exists(file.path(opt$outputDir, 'src'))){
    dir.create(file.path(opt$outputDir, 'src'), showWarnings = FALSE)
    invisible(sapply(list.files(file.path(opt$softwareDir, 'modules'), pattern = '*.R', full.names = TRUE), file.copy, to = file.path(opt$outputDir, 'src')))
    invisible(file.copy(opt$demultiplex_sampleDataFile, file.path(opt$outputDir, 'src', 'sampleData.tsv')))
    invisible(file.copy(opt$configFile, file.path(opt$outputDir, 'src', 'config.yml')))
    dir.create(file.path(opt$outputDir, 'src', 'version'), showWarnings = FALSE)
    invisible(file.copy(list.files(file.path(opt$softwareDir, 'version'), full.names = TRUE), file.path(opt$outputDir, 'src', 'version'), recursive = TRUE))
  }
}

# Start log.
updateMasterLog()

# Execute each module defined in opt$modules.
invisible(lapply(opt$modules, function(m){
  m <- unlist(strsplit(m, '\\s+'))
  args <- ifelse(length(m) == 2, m[2], '')

  r <- system(paste(opt$Rscript, file.path(opt$softwareDir, 'modules', paste0(m[1], '.R')), opt$configFile, args))

  if(r != 0){
    closeAllConnections()
    quit(save = "no", status = 0, runLast = FALSE)
  }
}))

closeAllConnections()
quit(save = "no", status = 0, runLast = FALSE)
