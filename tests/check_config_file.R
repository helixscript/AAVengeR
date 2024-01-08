#' Check whether variables in production yaml config file are the same as those in a table of known and validated variables.
#' Does not test values for variables EXCEPT module (should always be core and only core)

rm(list=ls())
library(dplyr)
library(tidyr)
library(stringr)
library(yaml)
library(testthat)
library(argparse)

aavenger_base_dir <- dirname(dirname(stringr::str_replace(commandArgs(trailingOnly = FALSE)[4], '--file=', '')))

source(file.path(aavenger_base_dir, 'tests', 'lib', 'helpers.R'))

input_args <- parse_inputs()

## ----- DEBUGGING PARAMS START----- ##
# input_args <- data.frame('run_type' = NA_character_, 'aavenger_dir_path' = NA_character_)
# input_args$run_type <- 'build'
# input_args$run_type <- 'test'
# aavenger_base_dir <- '/data/AAVengeR'
## ----- DEBUGGING PARAMS END----- ##

input_args$aavenger_dir_path <- aavenger_base_dir

if (input_args$run_type == 'build') {
  
  announce_run_steps(input_args$run_type, 'check_config_file.R', TRUE)
  
  production_config <- yaml::read_yaml(file.path(input_args$aavenger_dir_path, 'config.yml'))
  
  test_config_variables <- data.frame(
    'config_variable' = names(production_config)
  )
  
  tryCatch(
    {
    testthat::expect_equal(length(base::unique(test_config_variables$config_variable)), length(test_config_variables$config_variable), info = 'Check that there are no duplicate variables in yaml file. Error means there are duplicates.') 
    
    testthat::expect_equal(1, length(production_config$modules), info = 'Check that "core" is the only module listed. Error means that core is not the only module listed.') 
      
    testthat::expect_equal(TRUE, input_args$overwrite, info = 'Check that build mode is allowed to overwrite. Error means that build was not allowed to overwrite. Enable with --overwrite') 
    
    write.csv(test_config_variables, file.path(input_args$aavenger_dir_path, 'tests', 'test_data', 'test_config_variables.csv'), row.names = FALSE)
    
    print(paste0('New config variable reference table written to ', file.path(input_args$aavenger_dir_path, 'tests', 'test_data', 'test_config_variables.csv')))
    },
    error = function(e) {
      
      print(paste('Test failed with error:', conditionMessage(e)))
      
  })
  announce_run_steps(input_args$run_type, 'check_config_file.R', FALSE)
}

if (input_args$run_type == 'test') {
  
  announce_run_steps(input_args$run_type, 'check_config_file.R', TRUE)
  
  production_config <- yaml::read_yaml(file.path(input_args$aavenger_dir_path, 'config.yml'))
  
  test_config_variables <- read.csv(file.path(input_args$aavenger_dir_path, 'tests', 'test_data', 'test_config_variables.csv'))
  
  tryCatch(
    {
      testthat::expect_equal(0, length(setdiff(names(production_config), test_config_variables$config_variable)), info = 'Check that variables in in production yaml match known list of known variables. Failure means different number of variables.') 
      
      testthat::expect_equal(length(names(production_config)), length(test_config_variables$config_variable), info = 'Check no duplicates in production yaml variable names. Failure means there are duplicates.')
      
      print('Tests passed successfully!')
    },
    
    error = function(e) {
      print(paste('Test failed with error:', conditionMessage(e)))
      
    }
  )
  
  announce_run_steps(input_args$run_type, 'check_config_file.R', FALSE)
}




