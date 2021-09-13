library(yaml)

opt <- read_yaml('config.yml')

dir.create(opt$outputDir)
if(! dir.exists(opt$outputDir)) stop('Error - could not create the output directory.')

dir.create(file.path(opt$outputDir, 'tmp'))

source(file.path(opt$softwareDir, 'demultiplex.R'))
source(file.path(opt$softwareDir, 'createUniqueSampleFasta.R'))





Rscript <- '/home/opt/R-3.4.0/lib/R/bin/Rscript'

dir.create('out')
dir.create('out/demultiplexed_seqs')

system(paste(Rscript, 'demuliplex.R', 
             '--softwareDir /home/everett/projects/AAVengeR2',
             '--sampleConfigFile /home/everett/projects/AAVengeR_runs/configs/sampleConfig.Sabatino.tsv',
             '--index1seq /home/everett/projects/AAVengeR/data/testData/I1.test.fastq.gz',
             '--adriftReadsFile /home/everett/projects/AAVengeR/data/testData/R1.test.fastq.gz',
             '--anchorReadsFile /home/everett/projects/AAVengeR/data/testData/R2.test.fastq.gz',
             '--softwareDir /home/everett/projects/AAVengeR2',
             '--outputDir /home/everett/projects/AAVengeR2/out/demultiplexed_seqs'))

