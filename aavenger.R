library(yaml)
library(lubridate)

opt <- read_yaml('config.yml')


if(! dir.exists(opt$outputDir)){
  dir.create(file.path(opt$outputDir, opt$demultiplex_outputDir))
}

if(! dir.exists(opt$outputDir)){
  message('Error: Can not create the output directory.')
  q(save = 'no', status = 1, runLast = FALSE)
}

dir.create(file.path(opt$outputDir, 'tmp'))
if(! dir.exists(file.path(opt$outputDir, 'tmp'))){
  message('Error: Can not create the tmp directory.')
  q(save = 'no', status = 1, runLast = FALSE)
}

write(c(date(), floor(as.numeric(now()))), file.path(opt$outputDir, 'log'), append = FALSE)

# Execute each module defined in opt$modules.
invisible(lapply(opt$modules, function(m){
  system(paste(opt$Rscript, file.path(opt$softwareDir, m)))
  }))
