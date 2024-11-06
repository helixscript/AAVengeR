#!/usr/bin/Rscript

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
suppressPackageStartupMessages(library(optparse))

args <- commandArgs(trailingOnly=TRUE)

if(any(grepl('\\-\\-', args))){
  option_list = list(
    make_option(c("--list.available.genomes"), type="character", default=NULL, help="Display a list of available genomes.", action="store_true"),
    make_option(c("--list.installed.genomes"), type="character", default=NULL, help="Display a list of installed genomes.", action="store_true"),
    make_option(c("--install.genome"), type="character", default=NULL, help="install an available genome.", action="store"),
    make_option(c("--remove.genome"), type="character", default=NULL, help="Remove a genome from your installation.", action="store"))
  
  opt_parser = OptionParser(option_list=option_list, add_help_option = FALSE)
  opt = parse_args(opt_parser)
  
  if(length(opt) > 1) stop('Error - more than one argument was provided.')
 
  if('list.available.genomes' %in% names(opt)){
    o <- yaml::read_yaml('defaults.yml')
    f <- system(paste0('wget -qO- ', o$remoteDataURL, '/fileList'), intern = TRUE)
    f <- f[2:length(f)]
    f <- f[! grepl('fileList', f)]
    sizes <-  unlist(lapply(strsplit(f, '\\s+'), '[[', 5))
    sizes <- sub('G', ' GB', sizes)
    sizes <- sub('M', ' MB', sizes)
    genomes <- unlist(lapply(strsplit(f, '\\s+'), '[[', 9))
    genomes <- sub('\\.tar$', '', genomes)
    message('Available genomes:', pandoc.table.return(tibble(refGenome = genomes, dataSize = sizes), style = "simple", split.tables = Inf, plain.ascii = TRUE))
    q()
  }
  
  if('list.installed.genomes' %in% names(opt)){
    f <- list.files('data/referenceGenomes/blat', pattern = "*.2bit", full.names = TRUE)
    
    if(length(f) == 0){
      message('No installed genomes found.\n')
      q()
    }
    
    g <- unique(sub('\\.2bit', '', stringr::str_extract(f, '[a-zA-Z\\d]+.2bit')))
    
    tab <- bind_rows(lapply(g, function(x){
      f <- list.files('data', recursive = TRUE, pattern = paste0('^', x, '.'), full.names = TRUE)
      tibble(refGenome = x, dataFiles = length(f), spaceUsed = paste0(round(sum(file.size(f)) / 1024 / 1024 / 1024, 2), ' GB'))
    }))
    
    message(pandoc.table.return(tab, style = "simple", split.tables = Inf, plain.ascii = TRUE))
    
    q()
  }
  
  if('install.genome' %in% names(opt)){
    o <- yaml::read_yaml('defaults.yml')
    f <- system(paste0('wget -qO- ', o$remoteDataURL, '/fileList'), intern = TRUE)
    f <- f[2:length(f)]
    f <- f[! grepl('fileList', f)]
    genomes <- unlist(lapply(strsplit(f, '\\s+'), '[[', 9))
    genomes <- sub('\\.tar$', '', genomes)
    
    if(! opt$install.genome %in% genomes){
      message("Error - ", opt$install.genome, " is not in the list of available genomes.\n")
      q()
    }
    
    ff <- list.files('data/referenceGenomes/blat', pattern = "*.2bit", full.names = TRUE)
    gg <- unique(sub('\\.2bit', '', stringr::str_extract(ff, '[a-zA-Z\\d]+.2bit')))
    
    if(opt$install.genome %in% gg){
      message('Error - ', opt$install.genome, ' is already installed.')
      q()
    }
    
    tmpDir <- paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')
    dir.create(tmpDir)
    
    message('Downloading genome data.')
    system(paste0('wget -q -O ', tmpDir, '/archive.tar ', o$remoteDataURL, '/', opt$install.genome, '.tar'))
    
    message('Extracting data.')
    system(paste0('tar xvf ', tmpDir, '/archive.tar -C ', tmpDir))
    
    message('Installing data.')
    system(paste0('rsync -a ', tmpDir, '/data/ data/'))
    unlink(tmpDir, recursive = TRUE)
    
    message('done.\n')
    q()
  }
  
  if('remove.genome' %in% names(opt)){
    
    f <- list.files('data', pattern = paste0('^', opt$remove.genome, '\\.'), recursive = TRUE, full.names = TRUE)
    
    if(length(f) == 0){
      message('No genome files were found matching ', opt$remove.genome, '.\n')
      q()
    }
    
    print(f)
    message("Remove these files? (y/n)?")
    line <- readLines("stdin", n=1)
    
    if(line == 'y'){
      invisible(file.remove(f)) 
    } 
    
    message(opt$remove.genome, ' removed.')
    q()
  }
}

source(file.path(yaml::read_yaml(commandArgs(trailingOnly=TRUE)[1])$softwareDir, 'lib.R'))
opt <- loadConfig()
optionsSanityCheck()

# Create main output folder and make a copy of the source code and configurations.
createOuputDir()

# Check for third party software.
checkSoftware()

if(! opt$calledFromCore){
  dir.create(file.path(opt$outputDir, 'src'))
  invisible(sapply(list.files(opt$softwareDir, pattern = '*.R', full.names = TRUE), file.copy, to = file.path(opt$outputDir, 'src')))
  invisible(file.copy(opt$sampleConfigFile, file.path(opt$outputDir, 'src', 'sampleData.tsv')))
  invisible(file.copy(opt$configFile, file.path(opt$outputDir, 'src', 'config.yml')))
  dir.create(file.path(opt$outputDir, 'src', 'version'))
  invisible(file.copy(list.files(file.path(opt$softwareDir, 'version'), full.names = TRUE), file.path(opt$outputDir, 'src', 'version'), recursive = TRUE))
}


# Start log.
opt$defaultLogFile <- file.path(opt$outputDir, 'log')
logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
write(logo, opt$defaultLogFile, append = FALSE)
write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), opt$defaultLogFile, append = TRUE)

updateLog("Starting AAVengeR's module chain.")
updateLog("Each module will create its own output folder in the main output folder.")
updateLog("The log files for each module are stored within their output folders.")

# Execute each module defined in opt$modules.
invisible(lapply(opt$modules, function(m){
  m <- unlist(strsplit(m, '\\s+'))
  args <- ifelse(length(m) == 2, m[2], '')
  
  updateLog(paste0('Starting ', m, '.'))

  r <- system(paste(opt$Rscript, file.path(opt$softwareDir, paste0(m[1], '.R')), opt$configFile, args))
  
  if(r != 0){
    updateLog(paste0('Module ', m[1], ' failed. Please see logs for details.'))
    if(length(m) == 2) updateLog(paste0('Module arguments: ', args))
    q(save = 'no', status = 1, runLast = FALSE)
  }
}))

updateLog("Done.")
q(save = 'no', status = 0, runLast = FALSE)
