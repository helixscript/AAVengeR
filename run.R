
opt <- yaml::read_yaml('config.yml')

dir.create(opt$outputDir)
if(! dir.exists(opt$outputDir)) stop('Error - could not create the output directory.')
dir.create(file.path(opt$outputDir, 'tmp'))

source(file.path(opt$softwareDir, 'demultiplex.R'))
source(file.path(opt$softwareDir, 'createUniqueSampleFasta.R'))
source(file.path(opt$softwareDir, 'vectorFilter.R'))




source(file.path(opt$softwareDir, 'mapLeaderSequences.R'))



c('demultiplex.R', 'createUniqueSampleFasta.R', 'vectorFilter.R', 'mapLeaderSequences.R', 'alignReads.R', 'buildFragments.R', 'standardizeFragments.R', 'buildSites.R')

source(file.path(opt$softwareDir, 'demultiplex.R'))
source(file.path(opt$softwareDir, 'createUniqueSampleFasta.R'))
source(file.path(opt$softwareDir, 'alignReads.R'))
source(file.path(opt$softwareDir, 'buildFragments.R'))
source(file.path(opt$softwareDir, 'standardizeFragments.R'))
source(file.path(opt$softwareDir, 'mapLeaderSequences.R'))
source(file.path(opt$softwareDir, 'buildSites.R'))
