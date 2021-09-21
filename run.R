Rscript <- 'xxx'

modules <- c('demultiplex.R',    'createUniqueFasta.R',    'vectorFilter.R', 'mapLeaderSequences.R',    'alignReads.R', 
             'buildFragments.R', 'standardizeFragments.R', 'buildSites.R',   'mapSiteLeaderSequences.R')


opt <- yaml::read_yaml('config.yml')

if(dir.exists(opt$outputDir)) stop('Output directory already exists.')
dir.create(opt$outputDir)
dir.create(file.path(opt$outputDir, 'tmp'))



invisible(lapply(modules, function(m) system(paste(Rscript, file.path(opt$softwareDir, m)))))