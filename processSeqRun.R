suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(RMySQL))

source(file.path(this.path::this.dir(), 'lib.R'))

option_list = list(
  make_option(c("--runConfigFile"), type="character", default=NULL, help="R1 fasta file", metavar="character"),
  make_option(c("--outputDir"), type="character", default=NULL, help="R2 fasta file", metavar="character"),
  make_option(c("--dbgroup"), type="character", default='specimen_management', metavar="character"),
  make_option(c("--seqRunID"), type="character", default=NULL, metavar="character"),
  make_option(c("--sampleDataFormat"), type="character", default='AAVengeR', metavar="character"),            
  make_option(c("--anchorReadsFile"), type="character", default=NULL, metavar="character"),
  make_option(c("--adriftReadsFile"), type="character", default=NULL, metavar="character"),
  make_option(c("--index1ReadsFile"), type="character", default=NULL, metavar="character"),
  make_option(c("--sampleDataFile"), type="character", default=NULL, metavar="character"),
  make_option(c("--CPUs"), type="integer", default=10, metavar="integer"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

opt$CPUs <- 30
opt$runConfigFile <- '/home/ubuntu/projects/CART/190525_M00281_0495_000000000-CGY8F/SampleSheet.csv'
opt$outputDir <- '/home/ubuntu/projects/CART/190525_M00281_0495_000000000-CGY8F'
opt$seqRunID  <- '190525_M00281_0495_000000000-CGY8F'
opt$sampleDataFormat <- 'intSiteCaller'
opt$anchorReadsFile <- '/home/ubuntu/projects/CART/190525_M00281_0495_000000000-CGY8F/Undetermined_S0_L001_R2_001.fastq.gz'
opt$adriftReadsFile <- '/home/ubuntu/projects/CART/190525_M00281_0495_000000000-CGY8F/Undetermined_S0_L001_R1_001.fastq.gz'
opt$index1ReadsFile <- '/home/ubuntu/projects/CART/190525_M00281_0495_000000000-CGY8F/Undetermined_S0_L001_I1_001.fastq.gz'

if(is.null(opt$runConfigFile))    stop('Error - the inputFile paramter was not provided.')
if(is.null(opt$outputDir))        stop('Error - the outputFile paramter was not provided.')
if(is.null(opt$seqRunID))         stop('Error - the seqRunID paramter was not provided.')
if(is.null(opt$sampleDataFormat)) stop('Error - the sampleFormat paramter was not provided.')
if(is.null(opt$anchorReadsFile))  stop('Error - the anchorReadsFile paramter was not provided.')
if(is.null(opt$adriftReadsFile))  stop('Error - the adriftReadsFile paramter was not provided.')
if(is.null(opt$index1ReadsFile))  stop('Error - the index1ReadsFile paramter was not provided.')

if(! dir.exists(opt$outputDir)) dir.create(opt$outputDir)
if(! dir.exists(opt$outputDir)) stop('Error - could not create the output directory.')

if(! file.exists(opt$runConfigFile))   stop('Error - run configuration file does not exist.')
if(! file.exists(opt$anchorReadsFile)) stop('Error - the anchor read file does not exist.')
if(! file.exists(opt$adriftReadsFile)) stop('Error - the adrift read file does not exist.')
if(! file.exists(opt$index1ReadsFile)) stop('Error - the index1 read file does not exist.')
if(! opt$sampleDataFormat %in% c('AAVengeR', 'intSiteCaller')) stop('Error - sampleDataFormat must be set to either AAVengeR or intSiteCaller.')

# Retrieve all sample data from the specimen database.
dbConn  <- dbConnect(MySQL(), group = opt$dbgroup)
if(! base::exists('dbConn')) stop('Error - Could not connect to the sample database.')

sampleData <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)

o <- readLines(opt$runConfigFile)
i <- which(grepl('^\\[metaData\\]', o))

if(length(i) != 1) stop('Error - can not locate the [metaData] block in the runConfigFile.')

if(opt$sampleDataFormat == 'intSiteCaller'){
  r <- intSiteCallerConversions()
  o <- read.table(textConnection(o[(i+1):length(o)]), header = TRUE, sep = ',')

  if(! all(o$vectorSeq %in% names(r))) stop('Error - all of the vector names in the SampleSheet are not defined in script.')

  s <- strsplit(o$alias, '\\-')

  d <- tibble(sample = unlist(lapply(s, '[[', 1)),
              replicate = unlist(lapply(s, '[[', 2)),
              adriftReadLinkerSeq = o$linkerSequence,
              index1Seq = o$bcSeq,
              refGenome = o$refGenome,
              vectorFastaFile = unlist(lapply(o$vectorSeq, function(x) r[[x]][['vector']])),
              anchorReadStartSeq = unlist(lapply(o$vectorSeq, function(x) r[[x]][['startSeq']])),
              leaderSeqHMM = unlist(lapply(o$vectorSeq, function(x) r[[x]][['hmm']])),
              flags = 'integrase') %>%
       left_join(select(sampleData, Patient, SpecimenAccNum, Trial) %>%  dplyr::rename('trial' = 'Trial', 'subject' = 'Patient'), by = c('sample' = 'SpecimenAccNum'))
       
  d[is.na(d$trial),]$trial <- 'Control'
  d[is.na(d$subject),]$subject <- opt$seqRunID
  
} else {
  d <- as_tibble(read.table(textConnection(o[(i+1):length(o)]), header = TRUE, sep = ','))
}

readr::write_tsv(d, file.path(opt$outputDir,'sampleData.tsv'))

config <- yaml::read_yaml(file.path(this.path::this.dir(), 'config.yml'))

config$softwareDir <- this.path::this.dir()
config$core_CPUs <- as.integer(opt$CPUs)
config$outputDir <- opt$outputDir
config$databaseConfigGroup <- 'AAVengeR'
config$pullDBsubjectFrags <- TRUE
config$demultiplex_sampleDataFile  <- file.path(opt$outputDir,'sampleData.tsv')
config$demultiplex_anchorReadsFile <- opt$anchorReadsFile
config$demultiplex_adriftReadsFile <- opt$adriftReadsFile
config$demultiplex_index1ReadsFile <- opt$index1ReadsFile
config$buildStdFragments_createMultiHitClusters <- FALSE
config$modules <- list('core')

yaml::write_yaml(config, file.path(opt$outputDir,'config.yml'))

system(paste('Rscript', file.path(this.path::this.dir(), 'aavenger.R'), file.path(opt$outputDir,'config.yml')))
