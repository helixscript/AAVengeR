calculate_total_sites_detected <- function(path_to_sites) {
  return(readRDS(path_to_sites) %>%
    dplyr::group_by(subject, sample) %>%
    dplyr::summarize(count = dplyr::n()))
}

get_aavenger_result_files <- function(results_path) {
  file_list <- list.files(file.path(results_path, 'core'), recursive = TRUE, full.names = TRUE, pattern  = '.rds')
  df_file_list <- data.frame(
    'file_path' = file_list,
    'basename' = sapply(file_list, function(x) {basename(x)}),
    'relative_path' = sapply(file_list, function(x) {
      stringr::str_replace(x, paste0(results_path, '/'), '')
    })
  )
  rownames(df_file_list) <- seq(1, nrow(df_file_list))
  
  return(df_file_list)
}


announce_run_steps <- function(run_type, filename, is_start) {
  if (is_start) {
    print(paste0('Running ', filename, ' in ', run_type))
  } else {
    print(paste0('Finished running ', filename, ' in ', run_type))
  }
}

parse_inputs <- function() {
  parser <- argparse::ArgumentParser()
  
  # Define command-line arguments
  parser$add_argument("--run_type", type = "character", help = "Specify the run type (build or test)", choices = c('build', 'test'))
  parser$add_argument("--aavenger_dir_path", type = "character", help = "Specify the AAVengeR directory path")
  parser$add_argument("--output_dir", type = "character", help = "Specify the output directory path")
  parser$add_argument("--overwrite", action = 'store_true', default = FALSE, help = "Overwrite files when in build mode.")
  parser$add_argument("--number_of_cpus", type = 'integer', default = 8, help = "Number of CPUs to use for processing.")
  
  args <- parser$parse_args()
  
  return(args)
}


test_equal_order_by_readID <- function(path_to_rds_reference, path_to_rds_test) {
  testthat::expect_equal(
    readRDS(file.path(path_to_rds_reference)) %>%
      dplyr::arrange(readID),
    readRDS(file.path(path_to_rds_test)) %>%
      dplyr::arrange(readID),
    info = 'Check that output file pairs are identical. Failure means file pairs are not identical.')
}


make_core_module_integration_test_config_df <- function(aavenger_dir_path, number_of_cpus, output_dir) {

  test_config <- data.frame('mode' = 'integrase')
  
  # CONFIG START
  test_config$mode <- 'integrase'
  
  test_config$Rscript <- '/usr/bin/Rscript'
  test_config$softwareDir <- aavenger_dir_path
  
  test_config$outputDir <- output_dir # output
  
  test_config$core_CPUs <- number_of_cpus# 15
  test_config$core_maxPercentCPUs <- 25 # 25                                                         # For 8 CPUs per chunk using 47 CPUs -- (8 / 47)*100 =~ 17
  test_config$core_keepIntermediateFiles <- TRUE # TRUE
  
  test_config$demultiplex_seqRunID <- 'test_small' # xxx
  test_config$demultiplex_anchorReadsFile <- file.path(aavenger_dir_path, 'tests', 'test_data', 'small_R2.fastq.gz') # data/Undetermined_S0_R2_001.fastq.gz
  test_config$demultiplex_adriftReadsFile <- file.path(aavenger_dir_path, 'tests', 'test_data', 'small_R1.fastq.gz') # data/Undetermined_S0_R1_001.fastq.gz
  test_config$demultiplex_index1ReadsFile <- file.path(aavenger_dir_path, 'tests', 'test_data', 'small_I1.fastq.gz') # data/Undetermined_S0_I1_001.fastq.gz
  test_config$demultiplex_sampleDataFile <- file.path(aavenger_dir_path, 'tests', 'test_data', 'small_sample_data.tsv') #  data/sampleData.tsv
  test_config$demultiplex_level <- 'unique' # unique
  test_config$demultiplex_quickAlignFilter <- FALSE # FALSE                                              # (!) Designed for AAV work - will be slow for integrase work.
  test_config$demultiplex_useAdriftReadUniqueLinkers <- TRUE # TRUE
  test_config$demultiplex_requirePostUmiLinker <- TRUE # TRUE
  test_config$demultiplex_requirePostUmiLinker_maxMismatch <- 1 # 1
  
  test_config$processAdriftReadLinkerUMIs <- TRUE # TRUE
  
  test_config$modules <- 'core'
  
  test_config$compressDataFiles <-  FALSE
  
  test_config$demultiplex_CPUs <- number_of_cpus                                                            # [integer]   Number of CPUs to use for demultiplexing.
  test_config$demultiplex_outputDir <- 'demultiplex'                                              # [character] Demultiplexing ouput directory name that will be created in the output directory.
  test_config$demultiplex_RC_I1_barcodes <- FALSE                                               # [boolean]   Use the reverse compliment of I1 barcode sequences found in [sampleConfigFile].
  test_config$demultiplex_RC_I1_barcodes_auto <- TRUE                                           # [boolean]   Attempt to automatically determine best setting for demultiplex_RC_I1_barcodes. Will override [demultiplex_RC_I1_barcodes].
  test_config$demultiplex_qualtrim_code <- '5'                                                  # [character] ASCII Q-score threshold to fail sliding window read quality trimming.
  test_config$demultiplex_qualtrim_halfWidth <- 5                                               # [integer]   Half width of sliding window used to evaluate base call scores.
  test_config$demultiplex_qualtrim_events <- 2                                                  # [integer]   Number of events within a sliding that need to fall below the min. score to trigger trimming.
  test_config$demultiplex_qualtrim_minLength <- 50                                              # [integer]   Minimum length of anchor and adrift reads post quality trimming.
  test_config$demultiplex_correctGolayIndexReads <- TRUE                                        # [boolean]   Correct Golay barcodes using Golay alogrithm. Requires Python2 and numPy module installed.
  test_config$demultiplex_sequenceChunkSize <- 300000                                           # [integer]   Chunk size to use when breaking sequencing data into chunks for parallel processing.
  test_config$demultiplex_index1ReadMaxMismatch <- 1                                            # [integer]   Max. number of mismatches allowed in I1 barcodes.
  test_config$demultiplex_adriftReadLinkerBarcodeMaxMismatch <- 1                               # [integer]   (evaluated if demultiplex_useAdriftReadUniqueLinkers = TRUE) Max. number of mismatches allowed in R1 unique linker sequences.
  test_config$demultiplex_anchorReadStartSeq.maxMisMatch <- 1                                   # [integer]   (evaluated if anchorReadStartSeq column provided in sampleData file) Max. number of mismatches allowed in anchor read start sequences.
  test_config$demultiplex_exportFASTQ <- FALSE                                                  # [bollean]   Export demultiplexed reads as FASTQ to the demultiplex output folder.
  test_config$demultiplex_mergeSimilarReadPairsParams <- '-c 0.98 -d 0 -M 0 -G 1 -aS 0.98 -aL 0.98 -s 0.98 -n 10 -gap -1000'                                     # [string]    CD-HIT-EST parameters for clustering similiar read pairs.
  test_config$demultiplex_quickAlignFilter_minPercentID <- 97                                   # 
  test_config$demultiplex_quickAlignFilter_minMatches <- 30
  test_config$demultiplex_quickAlignFilter_minEstLeaderSeqLength <- 15
  
  
  test_config$prepReads_CPUs <- number_of_cpus                                                              # [integer]   Number of CPUs to use when preparing reads for alignment to a reference genome.
  test_config$prepReads_outputDir <- 'prepReads'                                                  # [character] prepReads ouput directory name that will be created in the output directory.
  test_config$prepReads_readsTable <-  'demultiplex/reads.rds'                                  # [file path] (@) Path to input reads data object created by demultiplex.R.
  test_config$prepReads_mergeSimilarReadPairs <- TRUE
  test_config$prepReads_excludeAnchorReadVectorHits <- TRUE                                     # [boolean]   Exclude anchor reads whoes ends align strongly to the vector sequence.
  test_config$prepReads_excludeAdriftReadVectorHits <- FALSE                                    # [boolean]   Exclude adrift reads whoes ends align strongly to the vector sequence.
  test_config$prepReads_vectorAlignmentTestLength <- 25                                         # [integer]   Length of NTs to from ends of reads to align to vectors to test for vector alignments.
  test_config$prepReads_vectorAlignmentTest_minPercentCoverage <- 90                            # [float]     Min. percent of {prepReads_vectorAlignmentTestLength} required for a vector test alignment to be considered. 
  test_config$prepReads_vectorAlignmentTest_minPercentID <- 90                                  # [float]     Min. percent seq id required for a vector test alignment to be considered. 
  test_config$prepReads_minAnchorReadLength <- 25                                               # [integer]   Min. length of anchor read post trim. Shorter reads would likely not align well to reference genomes.
  test_config$prepReads_minAdriftReadLength <- 25                                               # [integer]   Min. length of adrift read post trim. Shorter reads would likely not align well to reference genomes.
  test_config$prepReads_buildReadMaps_blastReconstruction <- TRUE                               # [boolean]   Use local alignments or reads to vector sequences to build leader sequence models. 
  test_config$prepReads_mapLeaderSeqsMinAlignmentLength <- 15                                   # [integer]   (*) min. length of local alignment to vector used when building leader seq models.
  test_config$prepReads_mapLeaderSeqsMinPercentID <- 95                                         # [float]     (evaluated if repReads_buildReadMaps_blastReconstruction = TRUE) min. % seq id of local alignment to vector used when building leader seq models.
  test_config$prepReads_mapLeaderSeqsMaxGapBetweenAlignments <- 5                               # [integer]   The max. distanct between local alignments to vector sequences that can be merged into single repLeader sequences.
  test_config$prepReads_buildReadMaps_minMapStartPostion <- 5                                   # [integer]   (if prepReads_buildReadMaps_blastReconstruction = FALSE) Use the best alignment begining at a position <= this value as the leader seque model.
  test_config$prepReads_useDefaultHMMsettings <- TRUE                                           # [boolean]   Use the default HMM settings, overrides other HMM settings in this configuration file.
  test_config$prepReads_HMMsearchReadStartPos <- 1                                              # [integer]   (^) Read start position of region to be evaluated by HMM.
  test_config$prepReads_HMMsearchReadEndPos <-  100                                             # [integer]   (^) Read end position of region to be evaluated by HMM.
  test_config$prepReads_HMMmaxStartPos <- 5                                                     # [integer]   (^) Max position on anchor read that a HMM hit can begin.
  test_config$prepReads_HMMminFullBitScore <- 20                                                # [float]     (^) Min. HMM score for a leader sequence to be considered.
  test_config$prepReads_HMMmatchEnd <- TRUE                                                     # [boolean]   (^) Require a match to the end of the HMM profile. 
  test_config$prepReads_HMMmatchTerminalSeq <- 'CA'                                              # [character] (^) Required NTs at the end of HMM hits. Enter 'none' to disable.
  test_config$prepReads_limitLeaderSeqsWithQuickAlignFilter <- TRUE                             # [boolean]   Always set to FALSE for non-AAV work. Reccomended for AAV work.
  
  
  test_config$alignReads_CPUs <- number_of_cpus                                                             # [integer]   Number of CPUs to use when aligning reads to a reference genome.
  test_config$alignReads_outputDir <- 'alignReads'                                                # [character] alignReads ouput directory name that will be created in the output directory.
  test_config$alignReads_inputFile <- 'prepReads/reads.rds'                                      # [file path] (@) Path to input reads RDS file.
  test_config$alignReads_aligner <- 'blat'
  test_config$alignReads_blat_fastMap <- TRUE
  test_config$alignReads_blatUseOocFile <- FALSE
  test_config$alignReads_genomeAlignment_anchorRead_maxStartPos <- 3                            # [integer]   Max. start position along trimmed anchor reads that an alignment can begin. (!) Set to a small value such as 3 when using HMMs, large value such as 250 for AAV analyses.
  test_config$alignReads_genomeAlignment_adriftRead_maxStartPos <- 3                            # [integer]   Max. start position along trimmed adrift reads that an alignment can begin.
  test_config$alignReads_genomeAlignment_minPercentID <- 98                                     # [float]     Min. % seq id required to accept an alignment to the reference genome.
  test_config$alignReads_genomeAlignment_blatStepSize <- 9                                      # [integer]   BLAT step size parameter. BLAT default 9.
  test_config$alignReads_genomeAlignment_blatTileSize <- 11                                     # [integer]   BLAT tile size parameter. BLAT default 11.
  test_config$alignReads_genomeAlignment_blatRepMatch <- 3000                                   # [integer]   BLAT repMatch parameter. Number of times a tile needs to be seen before masking.
  test_config$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned <- 5                        # [integer]   Max. number of NTs at the end of anchor reads not aligning to the genome.
  test_config$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned <- 5                        # [integer]   Max. number of NTs at the end of adrift reads not aligning to the genome.
  #             This value should be generous for AAV analysis, eg. 10, since ITR remnants can not always be fully modeled
  #             and trimmed before aligning to genomes. Unlaligned NTs may arise from over-reading into poorly ITR remnants.
  # (@) path relative to output directory
  test_config$buildFragments_CPUs <- number_of_cpus                                                          # [integer]   Number of CPUs to use when generating initial fragments.
  test_config$buildFragments_outputDir <- 'buildFragments'                                        # [character] buildFragments ouput directory name that will be created in the output directory.
  test_config$buildFragments_duplicateReadFile <- 'prepReads/duplicateReads.rds'                 # [file path] (@)(optional) Path to input duplicate reads RDS file. Typically created by prepReads module.
  test_config$buildFragments_anchorReadsAlignmentFile <- 'alignReads/anchorReadAlignments.rds'    # [file path] (@)Path to anchor reads alignment data RDS file.
  test_config$buildFragments_adriftReadsAlignmentFile <- 'alignReads/adriftReadAlignments.rds'    # [file path] (@)Path to adrift reads alignment data RDS file.
  test_config$buildFragments_idGroup_size <- 500                                                # [integer]   Number of read ids to group into initial fragments at one time. Grouping prevents overloading join functions.
  test_config$buildFragments_maxReadAlignments <- 1500                                          # [integer]   Max. number of alignments allowed for a given sequencing read. If both the anchor and adrift read in a pair exceed this value, the software will subsample alignments.
  test_config$buildFragments_salvageReadsBeyondMaxNumAlignments <- FALSE
  test_config$buildFragments_minFragLength <- 20                                                # [integer]   Min. length allowed for a fragment.
  test_config$buildFragments_maxFragLength <- 3000                                              # [integer]   Max. length allowed for a fragment.
  test_config$buildFragments_randomLinkerID_minReadCountToSegreagate <- 10                      # [integer]   When testing random linker ids for sample crossover, only consider resolving crossovers if a random id was seen this many times.
  test_config$buildFragments_randomLinkerID_minSingleSampleMajorityPercent <- 90                # [integer]   When testing random linker ids for sample crossover, a single sample must have this majority percentage for a correction to occur.
  
  
  test_config$buildStdFragments_CPUs <- number_of_cpus                                                       # [integer]   Number of CPUs to use when generating standardized fragments.
  test_config$buildStdFragments_inputFile <- 'buildFragments/fragments.rds'                       # [file path relative to output path] Path to initial fragment RDS file (ouput of buildFragments module).
  test_config$buildStdFragments_outputDir <- 'buildStdFragments'
  test_config$buildStdFragments_remnantClusterParams <- '-c 0.85 -d 0 -M 0 -G 1 -aS 0.90 -aL 0.90 -s 0.90 -n 4'
  test_config$buildStdFragments_autoPullTrialSamples <- FALSE                                   # [boolean]   Automatically pull additional trial samples from the database to co-standardize with data found in buildStdFragments_inputFile.
  test_config$buildStdFragments_randomIDdupReadMult <- 5                                        # [integer]   When evaluating UMI sequences, if the most common position id associated with anUMI is x times larger than the second most common position then all positions assoicated with an UMI are switched to the dominate position.
  test_config$buildStdFragments_randomIDdupAbundMult <- 3                                       # [integer]   When evaluating UMI sequences, if the most common position id associated with an UMI has a sonicLength abundance x times larger than the second most common position then all positions assoicated with an UMI are switched to the dominate position.
  test_config$buildStdFragments_intSiteStdWindowWidth <- 5                                      # [integer]   Half width of window centered on each intSite when looking for other sites to merge.
  test_config$buildStdFragments_breakPointStdWindowWidth <- 4                                   # [integer]   Half width of window centered on each sonic break when looking for other break points to merge.
  test_config$buildStdFragments_minReadsPerFrag <- 1                                            # [integer]   Min number of reads required to reatin a UMI specific fragment.
  test_config$buildStdFragments_createMultiHitClusters <- FALSE                                  # [boolean]   Create networks of sites for reads that align well to multiple locations in a reference genome.
  test_config$buildStdFragments_minUMIfreq <- 15
  
  test_config$buildSites_CPUs <- number_of_cpus 
  test_config$buildSites_inputFile <- 'buildStdFragments/stdFragments.rds'
  test_config$buildSites_outputDir <- 'buildSites'
  test_config$buildSites_dualDetectWidth <- 6
  
  test_config$predictPCRartifacts_CPUs <- number_of_cpus 
  test_config$predictPCRartifacts_inputFile <- 'core/sites.rds'
  test_config$predictPCRartifacts_outputDir <- 'predictPCRartifacts'
  test_config$predictPCRartifacts_adjacentSeqLength <- 14
  test_config$predictPCRartifacts_minReportHalfMatches <- 5
  test_config$predictPCRartifacts_minReportMatches <- 10
  test_config$predictPCRartifacts_maxAlnGaps <- 1
  test_config$predictPCRartifacts_wordSize <- 5
  test_config$predictPCRartifacts_gapOpeningPenalty <- 10
  test_config$predictPCRartifacts_gapExtensionPenalty <- 5
  test_config$predictPCRartifacts_addAfter <- 'posid'
  
  test_config$barcodeAssocLinkers_outputDir <- 'barcodeAssocLinkers'                              # [character] barcodeAssocLinkers ouput directory name that will be created in the output directory.
  test_config$barcodeAssocLinkers_sampleDataFile <-  file.path(aavenger_dir_path, 'tests', 'test_data', 'small_sample_data.tsv') #data/sampleData.tsv                        # [file path] (@) sampleData file path.
  test_config$barcodeAssocLinkers_adriftReadsFile <- file.path(aavenger_dir_path, 'tests', 'test_data', 'small_R1.fastq.gz') #data/Undetermined_S0_R1_001.fastq.gz       # [file path] (@) Adrift read FASTQ path.
  test_config$barcodeAssocLinkers_index1ReadsFile <- file.path(aavenger_dir_path, 'tests', 'test_data', 'small_I1.fastq.gz') #data/Undetermined_S0_I1_001.fastq.gz       # [file path] (@) I1 barcode read FASTQ path.
  test_config$barcodeAssocLinkers_nCodes <- 50                                                  # [integer]   Number of most abundant I1 barcodes to examine.
  test_config$barcodeAssocLinkers_excludePolyG_index1BarCodes <- TRUE                           # [boolean]   Exclude I1 barcodes that contain GGGGG (most likely PhiX).
  test_config$barcodeAssocLinkers_adriftReadUniqueLinkerLength <- 20                            # [integer]   First n NTs of adrift read NTs to consider as a unique linker sequence.
  
  test_config$mapSiteLeaderSequences_CPUs <- number_of_cpus
  test_config$mapSiteLeaderSequences_outputDir <- 'mapSiteLeaderSequences'
  test_config$mapSiteLeaderSequences_inputFile <- 'core/sites.rds'
  test_config$mapSiteLeaderSequences_addAfter <- 'repLeaderSeq'
  test_config$mapSiteLeaderSequences_minAllowableGap <- 5                                       # Min. space between blast hits allowed. Larger gaps will be filled with [x].
  test_config$mapSiteLeaderSequences_alignmentChunkSize <- 5000
  test_config$mapSiteLeaderSequences_minAlignmentPercentID <- 90
  test_config$mapSiteLeaderSequences_minAlignmentLength <- 15

  return(test_config)
  
}


