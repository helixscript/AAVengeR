#' Run AAVengeR in build or test mode. Build makes a new set of results to use as a reference.
#' Test mode tests the results of a new run against the reference results

rm(list=ls())
library(dplyr)
library(tidyr)
library(stringr)
library(yaml)
library(argparse)

aavenger_base_dir <- dirname(dirname(stringr::str_replace(commandArgs(trailingOnly = FALSE)[4], '--file=', '')))

source(file.path(aavenger_base_dir, 'tests', 'lib', 'helpers.R'))

input_args <- parse_inputs()

## ----- DEBUGGING PARAMS START----- ##
# input_args <- data.frame('run_type' = NA_character_, 'aavenger_dir_path' = NA_character_)
# input_args$run_type <- 'build'
# input_args$run_type <- 'test'
# aavenger_base_dir <- '/data/AAVengeR'
# input_args$output_dir <- '/data/fake_test2'
# input_args$number_of_cpus <- 4
# input_args$overwrite <- TRUE
# source(file.path(aavenger_base_dir, 'tests', 'lib', 'helpers.R'))

# sudo rm -R /data/fake_test2
# sudo cp -R /data/AAVengeR/tests/test_data/small_test_results /data/fake_test2
## ----- DEBUGGING PARAMS END----- ##

input_args$aavenger_dir_path <- aavenger_base_dir

if (input_args$run_type == 'build' & input_args$overwrite) {
  
  announce_run_steps(input_args$run_type, 'core_module_integration_test.R', TRUE)
  
  input_args$output_dir <- file.path(input_args$aavenger_dir_path, 'tests', 'test_data', 'small_test_results')
  
  # make config test file
  test_config <- make_core_module_integration_test_config_df(aavenger_dir_path = input_args$aavenger_dir_path, number_of_cpus = input_args$number_of_cpus, output_dir = input_args$output_dir)
  
  writeLines(yaml::as.yaml(as.list(test_config), column.major = TRUE), file.path(input_args$aavenger_dir_path, 'tests', 'test_data', 'small_config_temp.yaml'))
  
  # Run AAVengeR
  system(paste0('Rscript ', file.path(input_args$aavenger_dir_path, 'aavenger.R'), ' ', file.path(input_args$aavenger_dir_path, 'tests', 'test_data', 'small_config_temp.yaml')))
  
  # remove test config file
  file.remove(file.path(input_args$aavenger_dir_path, 'tests', 'test_data', 'small_config_temp.yaml'))
  
  announce_run_steps(input_args$run_type, 'core_module_integration_test.R', FALSE)
}

if (input_args$run_type == 'test') {
  
  announce_run_steps(input_args$run_type, 'core_module_integration_test.R', TRUE)
  
  # make config test file
  test_config <- make_core_module_integration_test_config_df(aavenger_dir_path = input_args$aavenger_dir_path, number_of_cpus = input_args$number_of_cpus, output_dir = input_args$output_dir)
  
  writeLines(yaml::as.yaml(as.list(test_config), column.major = TRUE), file.path(input_args$aavenger_dir_path, 'tests', 'test_data', 'small_config_temp.yaml'))
  
  system(paste0('Rscript ', file.path(input_args$aavenger_dir_path, 'aavenger.R'), ' ', file.path(input_args$aavenger_dir_path, 'tests', 'test_data', 'small_config_temp.yaml')))
  
  # Gather companion files from test run and reference run
  aavenger_reference_results_path <- file.path(input_args$aavenger_dir_path, 'tests', 'test_data', 'small_test_results')

  df_results <- rbind(
    get_aavenger_result_files(aavenger_reference_results_path) %>%
      dplyr::mutate('data_source' = 'reference'),
    get_aavenger_result_files(input_args$output_dir) %>%
      dplyr::mutate('data_source' = 'test')
    )
    
  tryCatch(
    {
    # Test that there is 1 file per file pair for each data source
    df_results_file_counts <- df_results %>%
      dplyr::group_by(data_source, relative_path) %>%
      dplyr::summarize(count = dplyr::n()) %>%
      tidyr::pivot_wider(names_from = 'data_source', values_from = 'count', values_fill = 0)
    
    lapply(split(df_results_file_counts, 1:nrow(df_results_file_counts)), function(row) {
      
      testthat::expect_equal(row$reference, row$test, info = 'Check that filenames and their numbers match. Failure means different number of files.')
      })
    
    rm(df_results_file_counts)
    
    # Test that the contents within each file pair are identical  
    df_results_relative_path <- df_results %>%
      tidyr::pivot_wider(names_from = 'data_source', values_from = 'file_path')  
    
    ## These files need to be ordered by readID before being compared since their order of addition to the table is not reproducible
    read_id_based_files <- df_results_relative_path$basename[grep("*CORE_TMP.rds", df_results_relative_path$basename , ignore.case = TRUE)]
    read_id_based_files <- c(read_id_based_files, 'adriftReadAlignments.rds', 'anchorReadAlignments.rds', 'fragments.rds', 'reads.rds')
    
    lapply(split(df_results_relative_path, 1:nrow(df_results_relative_path)), function(row) {
      
      # print(row$relative_path)
  
      if (row$basename %in% read_id_based_files) {
        test_equal_order_by_readID(file.path(row$reference), file.path(row$test))
      } else {
        testthat::expect_equal(readRDS(file.path(row$reference)), readRDS(file.path(row$test)), info = 'Check that output file pairs are identical. Failure means file pairs are not identical.')
      }
    })
    
    # Report total number of sites detected. Not necessary but useful for talking to user.
    df_results_relative_path <- df_results_relative_path %>%
      dplyr::filter(relative_path == 'core/sites.rds')
    
    print('total sites detected reference:')
    print(calculate_total_sites_detected(df_results_relative_path$reference))
    
    print('total sites detected test:')
    print(calculate_total_sites_detected(df_results_relative_path$test))
    
    print('Tests passed successfully!')
    },
    
    error = function(e) {
      print(paste('Test failed with error:', conditionMessage(e)))
      
    })  
  
  # remove test config file
  file.remove(file.path(input_args$aavenger_dir_path, 'tests', 'test_data', 'small_config_temp.yaml'))
  
  announce_run_steps(input_args$run_type, 'core_module_integration_test.R', FALSE)
}
