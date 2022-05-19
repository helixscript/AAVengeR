library(yaml)
library(lubridate)

configFile = commandArgs(trailingOnly=TRUE)

if(nchar(configFile) == 0){
  stop('Error -- the config file paramater was not provided.')
}

if(! file.exists(configFile)){
  stop('Error -- the config file could not be found.')
}

opt <- read_yaml(configFile)

if(! dir.exists(opt$outputDir)){
  dir.create(file.path(opt$outputDir))
  
  if(! dir.exists(opt$outputDir)){
    message('Error: Can not create the output directory.')
    q(save = 'no', status = 1, runLast = FALSE)
  }
}


dir.create(file.path(opt$outputDir, 'tmp'))
if(! dir.exists(file.path(opt$outputDir, 'tmp'))){
  message('Error: Can not create the tmp directory.')
  q(save = 'no', status = 1, runLast = FALSE)
}

write(c(date(), floor(as.numeric(now()))), file.path(opt$outputDir, 'log'), append = FALSE)

# Execute each module defined in opt$modules.
invisible(lapply(opt$modules, function(m){
  write(c(paste(now(), 'Starting', m)), file = file.path(opt$outputDir, 'log'), append = TRUE)
  r <- system(paste(opt$Rscript, file.path(opt$softwareDir, m), configFile))
  if(r != 0){
    write(c(paste(now(), 'module', m, 'failed, please see log for details.')), file = file.path(opt$outputDir, 'log'), append = TRUE) 
    q(save = 'no', status = 1, runLast = FALSE)
  }
}))

write(c(paste(now(), 'Done.')), file = file.path(opt$outputDir, 'log'), append = TRUE) 
q(save = 'no', status = 0, runLast = FALSE)