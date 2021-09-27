opt <- yaml::read_yaml('config.yml')

# Check for an existing output directory and then create it.
#if(dir.exists(opt$outputDir)) stop('Output directory already exists.')
dir.create(opt$outputDir)
dir.create(file.path(opt$outputDir, 'tmp'))

# Execute each module defined in opt$modules.
invisible(lapply(opt$modules, function(m){
  system(paste(opt$Rscript, file.path(opt$softwareDir, m)))
  }))
