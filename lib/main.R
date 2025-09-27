updateLog <- function(msg, logFile = NULL){
  if(is.null(logFile)) logFile <- opt$defaultLogFile
  msg <- paste0(base::format(Sys.time(), "%m.%d.%Y %l:%M%P"), "\t", msg)
  write(msg, file = logFile, append = TRUE)
}


remove_whitespace_in_quotes <- function(file_path) {
  # Read all lines
  lines <- readLines(file_path)
  
  # Collapse into one text blob
  text <- paste(lines, collapse = "\n")
  
  # Now, match quoted strings manually
  text_fixed <- stringr::str_replace_all(text, '"([^"]*)"', function(x) {
    # x comes in as one full match: "something inside"
    inner <- stringr::str_sub(x, 2, -2)  # Remove the outer quotes manually (first and last char)
    inner_nospace <- stringr::str_remove_all(inner, "\\s+")  # Remove all whitespace inside
    stringr::str_c('"', inner_nospace, '"')  # Re-add quotes
  })
  
  # Finally, split back into lines
  fixed_lines <- stringr::str_split(text_fixed, "\n")[[1]]
  
  return(fixed_lines)
}



updateMasterLog <- function(){
  lockFile <- file.path(opt$orgOutputDir, 'log.lock')
  
  if(! file.exists(lockFile)){
    
    write(date(), lockFile)
    
    f  <- file.path(opt$orgOutputDir, 'log.tmp')
    f2 <- file.path(opt$orgOutputDir, 'log')
    
    # Create log skeleton.
    orgOpt <- yaml::read_yaml(textConnection(remove_whitespace_in_quotes(file.path(opt$orgOutputDir, 'src', 'config.yml'))))
 
    logo <- readLines(file.path(opt$softwareDir, 'figures', 'ASCII_logo.txt'))
    o <- lapply(orgOpt$modules, function(x) return('none\n'))
    names(o) <- orgOpt$modules
    
    # Create core module subsections.
    if('core' %in% names(o)){
      o$core <- c("<core>", "{core demultiplex}", "none\n", "{core replicate job table}", "none\n", 
                  "{core subject job table}", "none\n")
      
      d <- readr::read_tsv(opt$demultiplex_sampleDataFile, show_col_types = FALSE)
      
      r <- unique(paste0('{core replicate ', d$trial, '~', d$subject, '~', d$sample, '~', d$replicate, '}'))
      r <- as.vector(rbind(r, 'none\n'))
      
      s <- unique(paste0('{core subject ', d$trial, '~', d$subject, '}'))
      s <- as.vector(rbind(s, 'none\n'))
      
      o$core <- c(o$core, r, s)
      
      z <- split(o$core, cumsum(grepl('^\\{', o$core)))
      names(z) <- unlist(lapply(z, function(x) gsub('\\{|\\}', '', x[1])))
      z <- lapply(z, function(x) x[-1])
      o$core <- z
    }
    
    if(any(grepl('outputDir', names(o)))){
      # Handle instances where outputDir is changed by parameter injection.
      names(o) <- mapply(function(x, y){
                    if(length(x) > 0){
                      return(sub('outputDir:', '', stringr::str_extract(x, 'outputDir\\:[^,^\\:^\\"]+')))
                    } else {
                      return(y)
                    }
                  }, stringr::str_extract_all(names(o), 'outputDir:\\S+'), names(o))
      
    }
    
    if(file.exists(file.path(opt$orgOutputDir, 'core', 'replicateJobTable'))){
      o$core$`core replicate job table` <-  readLines(file.path(opt$orgOutputDir, 'core', 'replicateJobTable'))
    }
    
    if(file.exists(file.path(opt$orgOutputDir, 'core', 'subjectJobTable'))){
      o$core$`core subject job table` <-  readLines(file.path(opt$orgOutputDir, 'core', 'subjectJobTable'))
    }
    
    # Clear out the master log file and rebuild it.
    write(logo, f, append = FALSE)
    write(paste0('version: ', readLines(file.path(opt$softwareDir, 'version', 'version')), "\n"), f, append = TRUE)
    
    # Cycle through the sections in skeleton. Take care to handle core subsections. 
    invisible(mapply(function(section, x){
      write(c('\n', logBanner(section)), f, append = TRUE)
      
      x <- x[-1]
      
      if(section == 'core' & is.list(x)){
        coreLogFiles <- list.files(file.path(opt$orgOutputDir, 'core'), recursive = TRUE, pattern = '^log$', full.names = TRUE)
        invisible(lapply(names(x), function(xx){
          
          if(grepl('demultiplex', xx)){
            if(length(coreLogFiles[grepl('demultiplex', coreLogFiles)]) == 1){
              if(unlist(x[xx])[1] == 'none\n') x[xx] <<- ''
              x$`core demultiplex` <<- readLines(coreLogFiles[grepl('demultiplex', coreLogFiles)])[-(1:11)]
            }
          }
          
          if(grepl('core replicate', xx)){
            logs <- coreLogFiles[grepl(paste0(stringr::str_extract(xx, '\\S+$'), '/'), coreLogFiles)]
            
            if(length(logs) > 0){
              if(unlist(x[xx])[1] == 'none\n') x[xx] <<- ''
              if(length(logs[grepl('prepReads', logs)]) == 1)  x[xx] <<- list(c(unlist(x[xx]), '\n# prepReads', paste0(c('#', rep('-', 99)), collapse = ''), readLines(logs[grepl('prepReads', logs)])[-(1:11)]))
              if(length(logs[grepl('alignReads', logs)]) == 1)  x[xx] <<- list(c(unlist(x[xx]), '\n# alignReads',  paste0(c('#', rep('-', 99)), collapse = ''), readLines(logs[grepl('alignReads', logs)])[-(1:11)]))
              if(length(logs[grepl('buildFragments', logs)]) == 1)  x[xx] <<- list(c(unlist(x[xx]), '\n# buildFragments',  paste0(c('#', rep('-', 99)), collapse = ''), readLines(logs[grepl('buildFragments', logs)])[-(1:11)]))
              
              attritionLogs <- list.files(opt$orgOutputDir, pattern = 'attritionLog.tsv', full.names = TRUE, recursive = TRUE)
              attritionLogs <- attritionLogs[grepl(paste0(stringr::str_extract(xx, '\\S+$'), '/'), attritionLogs)]
              
              if(length(attritionLogs) > 0){
                attrition <- bind_rows(lapply(attritionLogs, readr::read_tsv, show_col_types = FALSE))
                levs <- c("PRD1", "PRD2", "PRD3", "PRD4", "PRD5", "ALR1", "ALR2", "ALR3", "ALR4", "ALR5", "BFR1", "BSF1")
                attrition$label <- factor(attrition$label, levels = levs)
                attrition <- attrition[order(attrition$label),]
                totalReads <- ppNum(attrition[1,]$value)
                attrition <- attrition[! is.na(attrition$value),]
                attrition$value <- attrition$value / attrition[1,]$value
                title <- paste0('Read attrition of ', totalReads, ' demultiplexed reads.')
                x[xx] <<- list(c(unlist(x[xx]), '\n# readAttrition',  paste0(c('#', rep('-', 99)), collapse = ''), title, asciiPercentBarChart(attrition)))
              }
            }
          }
          
          if(grepl('core subject', xx)){
            logs <- coreLogFiles[grepl(paste0(stringr::str_extract(xx, '\\S+$'), '/'), coreLogFiles)]
            
            if(length(logs) > 0){
              if(unlist(x[xx])[1] == 'none\n') x[xx] <<- ''
              if(length(logs[grepl('buildStdFragments', logs)]) == 1)  x[xx] <<- list(c(unlist(x[xx]), '\n# buildStdFragments',  paste0(c('#', rep('-', 99)), collapse = ''), readLines(logs[grepl('buildStdFragments', logs)])[-(1:11)]))
              if(length(logs[grepl('buildSites', logs)]) == 1)  x[xx] <<- list(c(unlist(x[xx]), '\n# buildSites',  paste0(c('#', rep('-', 99)), collapse = ''), readLines(logs[grepl('buildSites', logs)])[-(1:11)]))
            }
          }
        }))
        
        invisible(mapply(function(section2, x2){
          write(c('\n', logBanner(section2, corner = '+', vert = '|', horiz = '-'),  x2), f, append = TRUE)
        }, names(x), x))
        
      } else {
        # Non-core sections.
        
        logFile <- list.files(opt$orgOutputDir, recursive = TRUE, pattern = 'log', full.names = TRUE)
        logFile <- logFile[grepl(section, logFile, fixed = TRUE)]
        
        if(length(logFile) == 1){
          x2 <- readLines(logFile)[-(1:11)]
        } else {
          x2 <- 'none'
        }
        
        write(c('\n', x2), f, append = TRUE)
      }
      
    }, names(o), o, SIMPLIFY = FALSE))
    
    invisible(system(paste('mv', f, f2)))
    invisible(file.remove(lockFile))
  }
}


logBanner <- function(text, width = NULL, corner = "#", horiz = "#", vert = "#") {
  if (is.null(width)) width <- nchar(text) + 10
  if (nchar(text) > (width - 2))  stop("Text is too long to fit in the specified width.")
  
  padding_total <- width - 2 - nchar(text)
  left_padding <- floor(padding_total / 2)
  right_padding <- ceiling(padding_total / 2)
  
  top_bottom <- paste0(corner, paste(rep(horiz, width - 2), collapse = ""), corner)
  middle <- paste0(vert, paste(rep(" ", left_padding), collapse = ""), text, paste(rep(" ", right_padding), collapse = ""), vert)
  
  return(c(top_bottom, middle, top_bottom))
}


asciiPercentBarChart <- function(df, max_width = 50, label_width = 5, bar_char = "#") {
  if (! all(c("label", "value") %in% names(df))) stop("Data frame must have columns named 'label' and 'value'.")
  
  r <- vector()
  for (i in seq_len(nrow(df))) {
    label <- format(df$label[i], width = label_width, justify = "left")
    percent_val <- round(df$value[i] * 100, 1)
    bar_len <- round(df$value[i] * max_width)
    bar <- paste(rep(bar_char, bar_len), collapse = "")
    r <- append(r, paste0(label, "|", bar, " ", percent_val, "%"))
  }
  
  r
}


startModule <- function(args){
  set.seed(1)
  configFile <- args[1]
  if(! file.exists(configFile)) stop('Error - the configuration file does not exists.')
  
  opt <- yaml::read_yaml(textConnection(remove_whitespace_in_quotes(configFile)))
  
  opt <- ensureDataTypes(opt)
  
  if(length(args) == 2){
    message('Injecting module...')
    invisible(lapply(strsplit(unlist(strsplit(args[2], '\\s*,\\s*')), '\\s*:\\s*'), function(x){
      val <- x[2]
      if(! suppressWarnings(is.na(as.numeric(val))) & grepl('\\.', val)) val <- as.numeric(val)
      if(! suppressWarnings(is.na(as.numeric(val))) & ! grepl('\\.', val)) val <- as.integer(val)
      message(paste0('Setting ', x[1], ' to ', val))
      opt[[x[1]]] <<- val
    })) 
  }
  
  opt$configFile <- configFile
  if(! 'orgOutputDir' %in% names(opt)) opt$orgOutputDir <- opt$outputDir
  
  optionsSanityCheck(opt)
}


optionsSanityCheck <- function(opt){
  userRequired <- c('mode', 'softwareDir', 'outputDir', 'database_configFile', 'database_configGroup',
                    'demultiplex_anchorReadsFile', 'demultiplex_adriftReadsFile', 
                    'demultiplex_index1ReadsFile', 'demultiplex_sampleDataFile', 'modules')
  
  if(! all(userRequired %in% names(opt))){
    msg <- paste0('Error, these user required parameters are missing from the configuration file: ', paste0(userRequired[! userRequired %in% names(opt)], collapse = ', '))
    message(msg)
    closeAllConnections()
    quit(save = "no", status = 0, runLast = FALSE) 
  }
  
  if(! opt$mode %in% c('integrase', 'AAV', 'transposase', 'manual')){
    message('Error, the mode parameter must be set to integrase, AAV, transposase, or manual.')
    closeAllConnections()
    quit(save = "no", status = 0, runLast = FALSE)
  }
    
  # Read in the default values and supplement user configs where appropriate.
  defaultConfig <- yaml::read_yaml(file.path(opt$softwareDir, 'data/defaults.yml'))
  
  allowedToBeAdded <- names(defaultConfig)[! names(defaultConfig) %in% userRequired]
  
  toAdd <- allowedToBeAdded[! allowedToBeAdded %in% names(opt)]
  
  if(length(toAdd) > 0){
    for(x in toAdd){
      opt[[x]] <- defaultConfig[[x]]
    }
  }
  
  # The core_CPU value can not exceed the number of system cores.
  cores <- parallel::detectCores()
  if(opt$core_CPUs > cores) opt$core_CPUs <- cores
  
  # Create list of all options containing 'CPUs' except core_CPUs.
  a <- names(opt)[grepl('_CPUs', names(opt))]
  a <- a[a != 'core_CPUs']
  
  # Do not allow any non-core module exceed the number of system cores.
  # Set their values to core_CPUs value if requested except for cases when the config
  # file is written by the core module since it balances the CPUs based on read counts.
  
  if(length(a) > 0){
    for(x in a){
      if(opt[[x]] > cores) opt[[x]] <- cores
      if(opt$core_applyCoreCPUsToAllModules & ! opt$calledFromCore) opt[[x]] <- opt$core_CPUs
    }
  }
  
  # Set mode specific values.
  if(grepl('integrase', opt$mode, ignore.case = TRUE)){
    opt$buildSites_integraseCorrectionDist <- 2
    opt$alignReads_genomeAlignment_anchorRead_maxStartPos <- 3
    opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned <- 5
    opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned <- 5
    opt$prepReads_HMMmatchEnd <- TRUE
    opt$prepReads_HMMmatchTerminalSeq <- 'CA'
    opt$prepReads_forceAnchorReadStartSeq <- FALSE
  } else if(grepl('AAV', opt$mode, ignore.case = TRUE)){
    opt$buildSites_integraseCorrectionDist <- 0
    opt$alignReads_genomeAlignment_anchorRead_maxStartPos <- 300
    opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned <- 5
    opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned <- 10 # 10 default
    opt$prepReads_forceAnchorReadStartSeq <- TRUE # Use provided start sequence if a leader seq model is not returned.
    opt$prepReads_HMMmatchEnd <- FALSE # Triggers leaderSeq extension in alignReads.R.
  } else if(grepl('transposase', opt$mode, ignore.case = TRUE)){
    opt$buildSites_integraseCorrectionDist <- 0
    opt$alignReads_genomeAlignment_anchorRead_maxStartPos <- 3
    opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned <- 5
    opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned <- 5
    opt$prepReads_HMMmatchEnd <- TRUE
    opt$prepReads_HMMmatchTerminalSeq <- 'TA'
  } else if(grepl('manual', opt$mode, ignore.case = TRUE)){
    # Do nothing.
  } else {
    stop('Error -- mode set to an unknown value.')
  }
  
  opt
}


ppNum <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)


createOuputDir <- function(){
  if(! dir.exists(opt$outputDir)){
    dir.create(file.path(opt$outputDir), showWarnings = FALSE)
    
    if(! dir.exists(opt$outputDir)){
      message('Error: Can not create the output directory.')
      closeAllConnections()
      quit(save = "no", status = 0, runLast = FALSE)
    }
  }
}


previousSampleDatabaseCheck <- function(samples){
  if(tolower(opt$database_configGroup) != 'none'){
      conn <- createDBconnection()
      
      r <- bind_rows(lapply(split(samples, 1:nrow(samples)), function(x){
             o <- dbGetQuery(conn, paste0("select * from samples where trial = '", x$trial, "' and subject = '", x$subject, "' and sample = '", x$sample, "' and replicate = '", x$replicate, "' and refGenome = '", x$refGenome, "'"))
             x <- select(x, trial, subject, sample, replicate, refGenome)
             x$inDatabase <- nrow(o) > 0
             x
           }))
    
      write(pander::pandoc.table.return(r, style = 'simple'), file.path(opt$outputDir, 'databaseSampleCheck.txt'))
    
      dbDisconnect(conn)
      
      if(any(r$inDatabase)){
        message('Error -- one or more sample replicates are already in the database.')
        message('See report located here: ', file.path(opt$outputDir, 'databaseSampleCheck.txt'))
        return(TRUE)
      } else {
        return(FALSE)
      }
  } else {
    return(FALSE)
  }
}


runSam2Psl <- function(samFile, outputFile, cl){
   o <- data.table(line = readLines(samFile))
   o$n <- ntile(1:nrow(o), CPUs)
   o <- split(o, o$n)
   gc()

   sam2psl <- function(x){
     system(paste0(file.path(opt$softwareDir, 'bin', 'sam2psl.py'), ' -i ',  x, ' -o ', paste0(x, '.out')))
     invisible(file.remove(x))
     x <- readLines(paste0(x, '.out'))
     invisible(file.remove(paste0(x, '.out')))
     x
   }

   f <- lapply(o, function(x){
     tmpFile <- paste0('sam2psl.', paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''))
     write(x$line, file.path(opt$outputDir, tmpFile))
     file.path(opt$outputDir, tmpFile)
   })

   write(unlist(parLapply(cl, f, sam2psl)), outputFile)
}


percentSysMemUsed <- function(){
  m <- as.integer(unlist(strsplit(system('free', intern = TRUE)[2], '\\s+'))[2:3])
  (m[2] / m[1]) * 100
}


lowMemoryException <- function(stepDesc = 'none', maxMem = 98){
  if(percentSysMemUsed() > maxMem){
    msg <- paste0('Error - the system has exceeded ', maxMem, '% memory usage -- stopping all processes. Step desc: ', stepDesc)
    updateLog(msg, logFile = file.path(opt$outputDir, 'Error.log'))
    updateMasterLog()
    closeAllConnections()
    quit(save = "no", status = 0, runLast = FALSE)
  }
}


waitForMemory <- function(stepDesc = 'none', minMem = 80, maxWaitSecs = 1200, sleepSecs = 10){
  gc()
  t = 0
  
  repeat{
    lowMemoryException(stepDesc = stepDesc)
    
    if(t > maxWaitSecs){
      msg <- paste0('Error - the software has waited for ', maxWaitSecs, 
                    ' seconds for at least ', 100 - minMem, '% of the system memory to be free.',
                    ' Stopping all processes. Step desc: ', stepDesc)
      updateLog(msg, logFile = file.path(opt$outputDir, 'Error.log'))
      updateMasterLog()
      closeAllConnections()
      q(save = 'no', status = 1, runLast = FALSE) 
    }
    
    m <- percentSysMemUsed()
    if(m <= minMem) break
    
    updateLog(paste0('Waiting for memory to be freed, total memory use currently ', sprintf("%.2f%%", m), 
                     '%, waiting for no more than ', minMem, '% usage - waited for ', t, ' seconds.'))
    
    Sys.sleep(sleepSecs)
    
    t <- t + sleepSecs
  }
}



tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }

# AAVengeR is provided with default configuration files for different types of analyses
# where all the expected options are provided. Removing options from the configuration file
# may break the software in some instances. Some modules, such as the core module, injects 
# non-public options and this function ensures that they are always set even when not needed
# so modules will not throw errors. 


# Helper function that creates the expected done files expected by the core module when modules fail.
# The core module would become stuck in a perpetual loop if the done files did not appear after a failure.

core_createFauxFragDoneFiles <- function(){
  if(! dir.exists(file.path(opt$outputDir, 'buildFragments'))) dir.create(file.path(opt$outputDir,'buildFragments'), showWarnings = FALSE)
  write(date(), file = file.path(opt$outputDir, 'buildFragments', 'fragments.done'))
  write(date(), file = file.path(opt$outputDir, 'buildFragments', 'xxx'))
}


core_createFauxSiteDoneFiles <- function(){
  if(! dir.exists(file.path(opt$outputDir, 'buildSites'))) dir.create(file.path(opt$outputDir, 'buildSites'), showWarnings = FALSE)
  write(date(), file = file.path(opt$outputDir,  'buildSites', 'sites.done'))
  write(date(), file = file.path(opt$outputDir,  'buildSites', 'xxx'))
}


# Check the system path for the presence of required system level software.

checkSoftware <- function(){
  s <- c('blat', 'cutadapt', 'hmmbuild', 'hmmsearch', 'blastn', 'makeblastdb', 'cd-hit-est', 'pear')
  
  o <- bind_rows(lapply(s, function(x){
         r <- tryCatch({
                system(paste0('which ', x), intern = TRUE)
              },
              error=function(cond) {
                return(NA)
              },
              warning=function(cond) {
                return(NA)
              })
           tibble(program = x, path = r)
        }))   
  
        o$program <- paste0('                       ', o$program)
        
        write(paste0('Paths of expected software packages:'), file = file.path(opt$outputDir, 'log'), append = TRUE)
        write(pandoc.table.return(o, style = "simple", split.tables = Inf, plain.ascii = TRUE, justify = 'left'), file.path(opt$outputDir, 'log'), append = TRUE)
  
       if(any(is.na(o$path))){
         write(c(paste('Error - one or more of the expected softwares are not in your search path.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
         updateMasterLog()
         closeAllConnections()
         quit(save = "no", status = 0, runLast = FALSE)
       } 
}


CD_HIT_clusters <- function(x, dir, params){
  f <- tmpFile()

  writeXStringSet(x, file.path(dir, paste0(f, '.fasta')))
  
  comm <- paste0("cd-hit-est -i ", file.path(dir, paste0(f, '.fasta')), 
                " -o ", file.path(dir, f), " -T ", opt$buildStdFragments_CPUs, 
                "  ", params, ' 1> /dev/null')
  
  m <- system(comm)
  
  if(m != 0) quitOnErorr(paste0('cd-hit-est failed in run dir: ', dir))
  
  r <- paste0(readLines(paste0(file.path(dir, f), '.clstr')), collapse = '')
  invisible(file.remove(list.files(dir, pattern = f, full.names = TRUE)))

  o <- unlist(strsplit(r, '>Cluster'))
  o[2:length(o)]
}


# Last Path Element -- return the last element of a file path delimited by slashes.
lpe <- function(x){
  o <- unlist(strsplit(x, '/'))
  o[length(o)]
}


# Convert a ShortRead object to a BioString DNA object and
# removing trailing lane information from reads ids.
shortRead2DNAstringSet <- function(x){
  r <- x@sread
  names(r) <- sub('\\s+.+$', '', as.character(x@id))
  r
}


# Helper function to wait for a file to appear on the file system.
waitForFile <- function(f, seconds = 1){
  repeat
  {
    Sys.sleep(seconds)
    if(file.exists(f)) break
  }
  return(TRUE)
}


demultiplex <- function(x){
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ShortRead))
 
  source(file.path(opt$softwareDir, 'lib', 'main.R'))
  
  chunk <- x$n[1]
  
  f <- file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', x$file[grepl('index1', x$file)])
  if(! file.exists(f)) waitForFile(f)
  index1Reads <- readFastq(f)
  invisible(file.remove(f))
  
  waitForMemory(stepDesc = 'Demultiplex', minMem = opt$system_minMemThreshold, maxWaitSecs = opt$system_minMemWaitTime, sleepSecs = opt$system_minMemSleepTime)
  
  f <- file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', x$file[grepl('anchor', x$file)])
  if(! file.exists(f)) waitForFile(f)
  anchorReads <- readFastq(f)
  invisible(file.remove(f))
  
  waitForMemory(stepDesc = 'Demultiplex', minMem = opt$system_minMemThreshold, maxWaitSecs = opt$system_minMemWaitTime, sleepSecs = opt$system_minMemSleepTime)
  
  f <- file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', x$file[grepl('adrift', x$file)])
  if(! file.exists(f)) waitForFile(f)
  adriftReads <- readFastq(f)
  invisible(file.remove(f))
  
  waitForMemory(stepDesc = 'Demultiplex', minMem = opt$system_minMemThreshold, maxWaitSecs = opt$system_minMemWaitTime, sleepSecs = opt$system_minMemSleepTime)
  
  anchorReads <- trimTailw(anchorReads, opt$demultiplex_qualtrim_events, opt$demultiplex_qualtrim_code, opt$demultiplex_qualtrim_halfWidth)
  anchorReads <- anchorReads[width(anchorReads) >= opt$demultiplex_qualtrim_minLength]
  if(length(anchorReads) == 0) return()
  
  adriftReads <- trimTailw(adriftReads, opt$demultiplex_qualtrim_events, opt$demultiplex_qualtrim_code, opt$demultiplex_qualtrim_halfWidth)
  adriftReads <- adriftReads[width(adriftReads) >= opt$demultiplex_qualtrim_minLength]
  if(length(adriftReads) == 0) return()
  
  index1Reads <- shortRead2DNAstringSet(index1Reads)
  anchorReads <- shortRead2DNAstringSet(anchorReads)
  adriftReads <- shortRead2DNAstringSet(adriftReads)
  
  reads <- syncReads(index1Reads, anchorReads, adriftReads)
  index1Reads <- reads[[1]];  anchorReads  <- reads[[2]]; adriftReads  <- reads[[3]]
  
  if(length(adriftReads) == 0) return()
  
  # Correct Golay encoded barcodes if requested.
  if(opt$demultiplex_correctGolayIndexReads){
    source(file.path(opt$softwareDir, 'lib', 'golay.R'))
    g <- correctGolay12(index1Reads)

    gg <- DNAStringSet(ifelse(g$uncorrectable, g$input, g$corrected))
    names(gg) <- names(index1Reads)

    percentChanged <- (sum(! as.character(index1Reads) == as.character(gg)) / length(index1Reads))*100

    index1Reads <- gg
    rm(gg)

    updateLog(paste0('Golay correction complete. ', sprintf("%.2f", percentChanged), '% of reads updated via Golay correction.'))
    
    reads <- syncReads(index1Reads, anchorReads, adriftReads)
  
    index1Reads <- reads[[1]];  anchorReads  <- reads[[2]]; adriftReads  <- reads[[3]]
  }
  
  
  # Loop through samples in sample data file to demultiplex and apply read specific filters.
  invisible(lapply(1:nrow(samples), function(r){
    r <- samples[r,]
      
    v0 <- rep(TRUE, length(anchorReads))
    if('anchorReadStartSeq' %in% names(r)){
      v0 <- vcountPattern(r$anchorReadStartSeq, subseq(anchorReads, 1, nchar(r$anchorReadStartSeq)), max.mismatch = opt$demultiplex_anchorReadStartSeq.maxMisMatch) == 1
    }
    
    # Create barcode demultiplexing vectors.
    v1 <- vcountPattern(r$index1Seq, index1Reads, max.mismatch = opt$demultiplex_index1ReadMaxMismatch) > 0
    
    log.report <- tibble(sample = r$uniqueSample, demultiplexedIndex1Reads = sum(v1))
    
    # Create break read linker barcode demultiplexing vector.
    v2 <- rep(TRUE, length(adriftReads))
    if(opt$demultiplex_useAdriftReadUniqueLinkers){
      testSeq <- substr(r$adriftReadLinkerSeq, r$adriftRead.linkerBarcode.start, r$adriftRead.linkerBarcode.end)
      v2 <- vcountPattern(testSeq, subseq(adriftReads, r$adriftRead.linkerBarcode.start, r$adriftRead.linkerBarcode.end), max.mismatch = opt$demultiplex_adriftReadLinkerBarcodeMaxMismatch) > 0
      log.report$demultiplexedLinkerReads <- sum(v2)
    } else {
      log.report$demultiplexedLinkerReads <- NA
    }
    
    # Test to see if any reads demultiplex to this row of the sample table and then subset reads to this sample.
    i <- Reduce(base::intersect, list(which(v0), which(v1), which(v2)))
    if(length(i) == 0){
      log.report$demultiplexedReads <- 0
    } else {
      reads <- syncReads(index1Reads[i], anchorReads[i], adriftReads[i])
      index1Reads <- reads[[1]];  anchorReads  <- reads[[2]]; adriftReads  <- reads[[3]]
      
      if(length(index1Reads) == 0){
        log.report$demultiplexedReads <- 0
      } else {
        writeFasta(subseq(adriftReads, r$adriftRead.linkerRandomID.start, r$adriftRead.linkerRandomID.end), 
                   file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0(r$uniqueSample, '.', chunk, '.randomAdriftReadIDs')), compress = opt$compressDataFiles)
        
        writeFasta(anchorReads, file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0(r$uniqueSample, '.', chunk, '.anchorReads')), compress = opt$compressDataFiles)
        writeFasta(adriftReads, file.path(opt$outputDir, opt$demultiplex_outputDir, 'tmp', paste0(r$uniqueSample, '.', chunk, '.adriftReads')), compress = opt$compressDataFiles)
        log.report$demultiplexedReads <- length(index1Reads)
      }
    }
    
    write.table(log.report, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = file.path(opt$outputDir, opt$demultiplex_outputDir, 'logs', paste0(r$uniqueSample, '.', chunk, '.logReport')))
  }))
}


blat <- function(y, ref, dir){
  f <- file.path(dir, tmpFile())
  write(paste0('>', y$id, '\n', y$seq), f)
  
  occFlag <- ''
  if(opt$alignReads_blatUseOocFile){
    occFlag <- paste0(' -ooc=', file.path(opt$outputDir, opt$alignReads_outputDir, paste0(opt$alignReads_genomeAlignment_blatTileSize, '.ooc'))) 
  }
  
  system(paste0('blat ', ref, ' ', f, ' ', paste0(f, '.psl'), 
                ' -tileSize=', opt$alignReads_genomeAlignment_blatTileSize, 
                ' -stepSize=', opt$alignReads_genomeAlignment_blatStepSize, 
                ' -repMatch=', opt$alignReads_genomeAlignment_blatRepMatch, 
                ' -out=psl -t=dna -q=dna -minScore=0 -minIdentity=0 -noHead -noTrimA',
                occFlag, ifelse(opt$alignReads_blat_fastMap, ' -fastMap', '')), ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  write('done', paste0(f, '.done'))
}


alignReadEndsToVector <- function(y){
  suppressPackageStartupMessages(library(Biostrings))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(lubridate))
  
  s <- DNAStringSet(y$testSeq)
  names(s) <- y$readID
  
  b <- blastReads(s, tmpDirPath = file.path(opt$outputDir, opt$prepReads_outputDir, 'tmp'), dbPath = file.path(opt$outputDir, opt$prepReads_outputDir, 'dbs', 'd'))
  
  if(nrow(b) > 0){
    b$alignmentLength <- b$qend - b$qstart + 1
    b <- dplyr::filter(b, pident >= opt$prepReads_vectorAlignmentTest_minPercentID, alignmentLength >= floor(opt$prepReads_vectorAlignmentTestLength * (opt$prepReads_vectorAlignmentTest_minPercentCoverage/100)), gapopen <= 1)
  }
  
  if(nrow(b) > 0){
    b <- left_join(b, data.table(qname = names(s), testSeq = as.character(s)), by = 'qname')
    b$start  <- ifelse(b$sstart > b$send, b$send, b$sstart)
    b$end    <- ifelse(b$sstart > b$send, b$sstart, b$send)
    b$strand <- ifelse(b$sstart > b$send, '-', '+')
    b <- group_by(b, qname) %>% dplyr::top_n(1, wt = bitscore) %>% dplyr::arrange(desc(strand)) %>% dplyr::slice(1) %>% ungroup()
  }
  
  b
}


reconstructSampleEnvironment <- function(dbConn, trial, subject, sample, replicate, refGenome, outputDir){
  o <- dbGetQuery(dbConn, paste0("select * from samples where trial = '", trial, "' and subject = '", subject, 
                                 "' and sample = '", sample, "' and replicate = '", replicate, "'"))
  if(nrow(o) == 1){
    writeBin(unserialize(o$environment[[1]]), file.path(outputDir, 'environment.tar.xz'))
  } else {
    message('Error - ', nrow(o), ' rows of data were returned when retrieving data for reconstructSampleEnvironment().')
  }
}


reconstructDBtable <- function(o, tmpDirPath){
  r <- tibble()
  
  if(nrow(o) > 0){
    r <- bind_rows(lapply(split(o, 1:nrow(o)), function(y){
           f <- tmpFile()
           writeBin(unserialize(y$data[[1]]), file.path(tmpDirPath, paste0(f, '.xz')))
           system(paste0('unxz ', file.path(tmpDirPath, paste0(f, '.xz'))))
           d <- readr::read_tsv(file.path(tmpDirPath, f), show_col_types = FALSE)
           invisible(file.remove(file.path(tmpDirPath, f)))
           d
         }))
  }
  
  r
}


pullDBsubjectFrags <- function(dbConn, trial, subject, tmpDirPath){
  o <- dbGetQuery(dbConn, paste0("select * from fragments where trial = '", trial, "' and subject = '", subject, "'"))
  reconstructDBtable(o, tmpDirPath)
}


pullDBsubjectSites <- function(dbConn, trial, subject, tmpDirPath){
  o <- dbGetQuery(dbConn, paste0("select * from sites where trial = '", trial, "' and subject = '", subject, "'"))
  reconstructDBtable(o, tmpDirPath)
}


syncReads <- function(...){
  arguments <- list(...)
  n <- Reduce(base::intersect, lapply(arguments, names))
  lapply(arguments, function(x) x[match(n, names(x))])
}


createDBconnection <- function(){
  suppressPackageStartupMessages(library(RMariaDB))
  
  tryCatch({
    dbConnect(RMariaDB::MariaDB(), group = opt$database_configGroup, default.file = opt$database_configFile)
  },
  error=function(cond) {
    quitOnErorr('Error - could not connect to the database.')
  })
}


uploadSitesToDB <- function(sites){
  suppressPackageStartupMessages(library(RMariaDB))
  
  updateLog('Writing sites to the database.')

  if(! dir.exists(file.path(opt$outputDir, 'tmp'))) dir.create(file.path(opt$outputDir, 'tmp'), showWarnings = FALSE)
  
  conn <- createDBconnection()
  
  invisible(lapply(split(sites, paste(sites$trial, sites$subject, sites$sample, sites$refGenome)), function(x){
    updateLog(paste0('Uploading sites for: ', sites$trial[1], '/', sites$subject[1], '/', sites$sample[1], '/',sites$refGenome[1]))
    
    dbExecute(conn, paste0("delete from sites where trial='", x$trial[1], "' and subject='", x$subject[1],
                           "' and sample='", x$sample[1], "' and refGenome='", x$refGenome[1], "'"))
    
    f <- tmpFile()
    readr::write_tsv(x, file.path(opt$outputDir, 'tmp', f))
    system(paste0('xz --best ', file.path(opt$outputDir, 'tmp', f)))
    
    fp <- file.path(opt$outputDir, 'tmp', paste0(f, '.xz'))
    
    updateLog(paste0('Reading tar ball as a byte stream (', file.size(fp), ' bytes).'))
    
    tab <- readBin(fp, "raw", n = as.integer(file.info(fp)["size"])+100)
    
    invisible(file.remove(list.files(file.path(opt$outputDir, 'tmp'), pattern = f, full.names = TRUE)))
    
    updateLog('Pushing byte object to database.')
    
    r <- dbExecute(conn,
                   "insert into sites values (?, ?, ?, ?, ?)",
                   params = list(x$trial[1], x$subject[1], x$sample[1], x$refGenome[1], list(serialize(tab, NULL))))
    if(r == 0){
      quitOnErorr(paset0('Error -- could not upload site data for ', x$trial[1], '~', x$subject[1], '~', x$sample[1], ' to the database.'))
    } else {
      updateLog(paste0('Uploaded sites data for ',  x$trial[1], '~', x$subject[1], '~', x$sample[1], ' to the database.'))
    }
  }))
  
  dbDisconnect(conn)
}


loadSamples <- function(){
  
  samples <- readr::read_tsv(opt$demultiplex_sampleDataFile, show_col_types = FALSE, comment = "#")
  
  if(nrow(samples) == 0){
    updateLog('Error - no lines of information was read from the sample configuration file.')
    updateMasterLog()
    closeAllConnections()
    quit(save = "no", status = 0, runLast = FALSE) 
  }
  
  if('refGenome' %in% names(samples)){
    if(! all(sapply(unique(file.path(opt$softwareDir, 'data', 'referenceGenomes', 'blat', paste0(samples$refGenome, '.2bit'))), file.exists))){
      updateLog("Error - one or more blat database files could not be found in AAVengeR's data/referenceGenomes directory.")
      updateMasterLog()
      closeAllConnections()
      quit(save = "no", status = 0, runLast = FALSE) 
    }
  } else {
    samples$refGenome <- NA
  }
  
  if('vectorFastaFile' %in% names(samples)){
    if(! all(sapply(unique(file.path(opt$softwareDir, 'data', 'vectors', samples$vectorFastaFile)), file.exists))){
      updateLog("Error - one or more vector FASTA files could not be found in AAVengeR's data/vectors directory.")
      updateMasterLog()
      closeAllConnections()
      quit(save = "no", status = 0, runLast = FALSE) 
    }
  } else {
    samples$vectorFastaFile <- NA
  }
  
  if('leaderSeqHMM' %in% names(samples)){
    samples$leaderSeqHMM <- file.path(opt$softwareDir, 'data', 'hmms', samples$leaderSeqHMM)
    
    if(! all(sapply(unique(samples$leaderSeqHMM), file.exists))){
      updateLog("Error - one or more leader sequence HMM files could not be found in AAVengeR's data/hmms directory.")
      updateMasterLog()
      closeAllConnections()
      quit(save = "no", status = 0, runLast = FALSE) 
    }
  }
  
  if(! 'leaderSeqHMM' %in% names(samples) & opt$mode == 'integrase'){
      updateLog("Error - a leaderSeqHMM column is required when running in integrase mode.")
      updateMasterLog()
      closeAllConnections()
      quit(save = "no", status = 0, runLast = FALSE) 
  }
  
  if(! 'anchorReadStartSeq' %in% names(samples) & opt$mode == 'AAV'){
    updateLog("Error - a anchorReadStartSeq column is required when running in AAV mode.")
    updateMasterLog()
    closeAllConnections()
    quit(save = "no", status = 0, runLast = FALSE) 
  }
  
  # Create missing columns and set to NA.
  if(! 'adriftRead.linkerBarcode.start' %in% names(samples)){
    samples$adriftRead.linkerBarcode.start <- NA
    samples$adriftRead.linkerBarcode.end   <- NA
  }
  
  if(! 'adriftRead.linkerRandomID.start' %in% names(samples)){
    samples$adriftRead.linkerRandomID.start <- NA
    samples$adriftRead.linkerRandomID.end   <- NA
  }
  
  samples <- bind_rows(lapply(split(samples, 1:nrow(samples)), function(x){
    Ns <- stringr::str_extract(x$adriftReadLinkerSeq, 'N+')
    
    # Detecting Ns in the adrift linker and having adriftRead.linkerBarcode.start 
    # set to NA triggers the determination of positions.
    
    if(! is.na(Ns) & is.na(x$adriftRead.linkerBarcode.start)){  
      o <- stringr::str_locate(x$adriftReadLinkerSeq, Ns)
      x$adriftRead.linkerBarcode.start  <- 1
      x$adriftRead.linkerBarcode.end    <- o[1,1] - 1
    }
    
    # Detecting Ns in the adrift linker and having adriftRead.linkerRandomID.start 
    # set to NA triggers the determination of positions.
    
    if(! is.na(Ns) & is.na(x$adriftRead.linkerRandomID.start)){ 
      o <- stringr::str_locate(x$adriftReadLinkerSeq, Ns)
      x$adriftRead.linkerRandomID.start  <- o[1,1]
      x$adriftRead.linkerRandomID.end    <- o[1,2]
    }
    x
  }))
  
  requiredColumns <- c("trial", "subject", "sample", "replicate",  
                       "adriftReadLinkerSeq", "index1Seq", "refGenome", "vectorFastaFile", "flags")
  
  if(! all(requiredColumns %in% names(samples))){
    missingCols <- paste0(requiredColumns[! requiredColumns %in% names(samples)], collapse = ', ')
    updateLog(paste0('Error - the following columns were missing from the sample data file: ', missingCols))
    closeAllConnections()
    quit(save = "no", status = 0, runLast = FALSE)
  }
  
  if(any(grepl('\\s|~|\\||\\.', paste0(samples$trial, samples$subject, samples$sample, samples$replicate)))){
    updateLog("Error -- spaces, tildas (~), pipes (|), and dots (.) are reserved characters and can not be used in the trial, subject, sample, or replicate sample configuration columns.")
    updateMasterLog()
    closeAllConnections()
    quit(save = "no", status = 0, runLast = FALSE)
  }
  
  if(opt$demultiplex_useAdriftReadUniqueLinkers){
    if(any(is.na(samples$adriftRead.linkerBarcode.start)) | any(is.na(samples$adriftRead.linkerBarcode.end))){
      updateLog('Error - adriftRead.linkerBarcode.start or adriftRead.linkerBarcode.end is set to NA when requesting demultiplex_useAdriftReadUniqueLinkers')
      updateMasterLog()
      closeAllConnections()
      quit(save = "no", status = 0, runLast = FALSE)
    }
  }
  
  if(any(is.na(samples$adriftRead.linkerRandomID.start)) | any(is.na(samples$adriftRead.linkerRandomID.end))){
      updateLog('Error - adriftRead.linkerRandomID.start or adriftRead.linkerRandomID.end is set to NA.')
      updateMasterLog()
      closeAllConnections()
      quit(save = "no", status = 0, runLast = FALSE)
  }
  
  samples$uniqueSample <- paste0(samples$trial, '~', samples$subject, '~', samples$sample, '~', samples$replicate)
  
  if(any(duplicated(samples$uniqueSample))){
    updateLog('Error -- There are one ore more duplications of trial~subject~sample~replicate ids in the sample data file.')
    updateMasterLog()
    closeAllConnections()
    quit(save = "no", status = 0, runLast = FALSE)
  }
  
  samples
}


golayCorrection <- function(x, tmpDirPath = NA){
  suppressPackageStartupMessages(library(ShortRead))
  suppressPackageStartupMessages(library(dplyr))
  
  f <- file.path(tmpDirPath, paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') )
  writeFasta(x, file = f)
  
  system(paste('python2', file.path(opt$softwareDir, 'bin', 'golayCorrection', 'processGolay.py'), f))
  invisible(file.remove(f))
  
  corrected <- readFasta(paste0(f, '.corrected'))
  invisible(file.remove(paste0(f, '.corrected')))
  
  a <- left_join(tibble(id = names(x), seq = as.character(x)),
                 tibble(id = as.character(corrected@id), seq2 = as.character(sread(corrected))), by = 'id')
  a <- rowwise(a) %>% mutate(editDist = as.integer(adist(seq, seq2))) %>% ungroup()
  
  i <- which(a$editDist <= 2)
  a[i,]$seq <- a[i,]$seq2
  
  if(! all(names(x) == a$id)){
    updateLog('Error - There was an ordering error during the Golay correction step.')
    updateMasterLog()
    closeAllConnections()
    quit(save = "no", status = 0, runLast = FALSE) 
  }
  
  r <- DNAStringSet(a$seq)
  names(r) <- a$id
  r
}


buildRearrangementModel <- function(b, seqsMinAlignmentLength){
  r <- vector()
  counter <- 0
  b <- b[(b$qend - b$qstart + 1) >= seqsMinAlignmentLength,]
  b$qstartBin <- cut(b$qstart, c(seq(0, max(b$qstart), 5), Inf), labels = FALSE)
  b <- arrange(b, qstartBin, desc(bitscore))

  while(nrow(b) != 0){
    counter <- counter + 1
    b <- arrange(b, qstart, evalue)
    x <- b[1,]
    r <- paste0(r, ';', x$qstart, '..', x$qend, '[', x$sstart2, x$strand, x$send2, ']')

    b <- b[b$qstart >= x$qend - 5 & b$qstartBin != x$qstartBin,] 
    b <- arrange(b, qstartBin, desc(bitscore))
    
    if(counter == 1000) break
  } 
  
  sub('^;', '', r)
}


blast2rearangements <- function(b, maxMissingTailNTs, minLocalAlignmentLength){
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(IRanges))
  suppressPackageStartupMessages(library(data.table))
  
  bind_rows(lapply(split(b, b$qname), function(b){
    
    # Sort BLAST results by query start position and evalue (low to high).
    b$strand <- ifelse(b$send < b$sstart, '-', '+')
    
    # Alignment to the negative strand will result in the subject end to come before the start.
    # Switch it back so that they are sequential.
    b$sstart2 <- ifelse(b$sstart > b$send, b$send, b$sstart)
    b$send2   <- ifelse(b$sstart > b$send, b$sstart, b$send)
    
    r <- tibble(qname = b$qname[1], rearrangement = buildRearrangementModel(b, minLocalAlignmentLength))
    
    k <- unlist(stringr::str_extract_all(r$rearrangement, '\\.\\.\\d+'))
    lastRangeEnd <- as.integer(sub('\\.\\.', '', stringr::str_extract(k[length(k)], '\\.\\.\\d+')))
    
    # Add missing internal segments.
    o <- unlist(strsplit(r$rearrangement, ';'))
    
    if(length(o) > 1){
      for(i in c(2:length(o))){
        n1 <- as.integer(sub('\\[', '', stringr::str_extract(o[i-1], '\\d+\\[')))
        n2 <- as.integer(stringr::str_extract(o[i], '^\\d+'))
        if(n2 - n1 >= 10) o[i-1] <- paste0(o[i-1], ';', n1+1, '..', n2-1, '[x]')
      }
      
      r$rearrangement <- paste0(o, collapse = ';')
    }
      
    # Add an [x] segment to the end of models if it falls short of the end of the analyzed sequence.
    if(b$qlen[1] - lastRangeEnd >= maxMissingTailNTs) r$rearrangement <- paste0(r$rearrangement, ';', (lastRangeEnd+1), '..', b$qlen[1], '[x]')

    r
  }))
}


captureHMMleaderSeq <- function(reads, hmm, tmpDirPath = NA){
  outputFile <- file.path(tmpDirPath, tmpFile())
  
  if(opt$prepReads_useDefaultHMMsettings){
    localOpt <- opt # Force local copy opt opt.
    
    f <- sub('\\.hmm$', '.settings', hmm)
    if(! file.exists(f)) stop(paste0('Error -- hmm settings file: ', f, ' does not exists'))
    
    # Update local copy.
    y <- yaml::read_yaml(f)  
    for(i in names(y)){
      localOpt[[i]] <- y[[i]]
    }
    
    # rename local copy to original which will now be local and updated.
    opt <- localOpt
  }
  
  endPos <- ifelse(width(reads) > opt$prepReads_HMMsearchReadEndPos, opt$prepReads_HMMsearchReadEndPos, width(reads))
  writeXStringSet(subseq(reads, opt$prepReads_HMMsearchReadStartPos, endPos), outputFile)
  
  ### comm <- paste0('nhmmer --max -T 0 --incT 0 --tblout ', outputFile, '.tbl ', hmm, ' ', outputFile, ' > ', outputFile, '.hmmsearch')
  comm <- paste0('nhmmer --F1 1 --F2 1 --F3 1 -T -5 --incT -5 --nobias --popen 0.15 --pextend 0.05 --tblout ', outputFile, '.tbl ', hmm, ' ', outputFile, ' > ', outputFile, '.hmmsearch')
  system(comm)
  
  o <- readr::read_table(paste0(outputFile, '.tbl'), col_names = FALSE, col_types = NULL, comment = "#", show_col_types = FALSE)
  
  if(nrow(o) == 0) return(tibble())
  
  names(o) <- c('targetName', 'targetAcc', 'queryName', 'queryAcc', 'hmmStart', 'hmmEnd', 'targetStart', 'targetEnd', 'envStart', 'envEnd', 'seqLength', 'strand', 'fullEval', 'fullScore', 'bias', 'desc')
  
  # Collapse duplicate hits.
  o <- group_by(o, targetName) %>% dplyr::slice_max(fullScore, n = 1, with_ties = FALSE) %>% ungroup()

  # Subset the data based on user scoring thresholds.
  o <- subset(o, targetStart <= opt$prepReads_HMMmaxStartPos & fullScore >= as.numeric(opt$prepReads_HMMminFullBitScore))
  
  if(nrow(o) == 0) return(tibble())
  
  h <- readLines(hmm)
  hmmLength <- as.integer(unlist(strsplit(h[grepl('^LENG', h)], '\\s+'))[2])
  hmmName <- unlist(strsplit(h[grepl('^NAME', h)], '\\s+'))[2]
  
  # If requested, limit HMM hits to those with alignments near the end of the HMM.
  if(opt$prepReads_HMMmatchEnd) o <- o[abs(hmmLength - o$hmmEnd) <= opt$prepReads_HMMmatchEndRadius,]
  
  if(nrow(o) == 0) return(tibble())
  
  # Limit reads to those with matching HMM hits.
  reads2 <- reads[names(reads) %in% o$targetName]
  
  rm(reads)
  suppressMessages(gc())
  
  # Arrange reads to match o data frame.
  reads2 <- reads2[match(o$targetName, names(reads2))]
  
  if(length(reads2) == 0) return(tibble())
  
  if(! grepl('none', opt$prepReads_HMMmatchTerminalSeq, ignore.case = TRUE)){
    # Retrieve the NTs from sequences that correspond to the end of the HMM alignment range (not Env range).
    terminalSeqs <- unname(as.character(subseq(reads2, o$targetEnd-(nchar(opt$prepReads_HMMmatchTerminalSeq)-1), o$targetEnd)))
    
    a <- reads2[terminalSeqs == opt$prepReads_HMMmatchTerminalSeq]    # matches expected terminal sequence.
    b <- reads2[! terminalSeqs == opt$prepReads_HMMmatchTerminalSeq]  # does not match terminal sequence.
    
    if(length(b) > 0 & opt$prepReads_HMMmatchEndRadius >= 1){
      k <- select(o, targetName, targetEnd) %>% dplyr::filter(targetName %in% names(b))
      k <- k[match(k$targetName, names(b)),]
      
      k$start <- k$targetEnd - opt$prepReads_HMMmatchEndRadius
      k$end   <- k$targetEnd + opt$prepReads_HMMmatchEndRadius
      
      k$readLength <- width(b)
      k <- k[k$end < k$readLength,]
      b <- b[names(b) %in% k$targetName]  # Limit b to k.
      
      k$seq <- as.character(subseq(b, k$start, k$end))
      k$matchEnd <- unlist(lapply( stringr::str_locate_all(k$seq, opt$prepReads_HMMmatchTerminalSeq), function(x){ ifelse(length(x) == 0, NA, max(x[,2])) }))
      k <- k[! is.na(k$matchEnd),]        # Remove entries without a match to the expected terminal seq.
      b <- b[names(b) %in% k$targetName]  # Limit b to k.
      
      if(nrow(k) > 0){
        k$newTargetEnd <- k$targetEnd + (k$matchEnd - opt$prepReads_HMMmatchEndRadius - 1)
        k$newSeq <- as.character(subseq(b, 1, k$newTargetEnd))
        i <- match(k$targetName, o$targetName)
        o[i,]$targetEnd <- k$newTargetEnd
      }             
    } else {
      b <- DNAStringSet()
    }
    
    reads2 <- c(a, b)
  }
  
  o <- o[o$targetName %in% names(reads2),]
  
  o <- o[match(names(reads2), o$targetName),]
  
  tibble(id = names(reads2),
         LTRname = hmmName,
         LTRseq = as.character(subseq(reads2, 1, o$targetEnd)))
}


blastReads <- function(reads, wordSize = 5, evalue = 100, tmpDirPath = NA, dbPath = NA){
  suppressPackageStartupMessages(library(Biostrings))
  suppressPackageStartupMessages(library(dplyr))
  
  f <- tmpFile()
  writeXStringSet(reads,  file.path(tmpDirPath, paste0(f, '.fasta')))
  
  comm <- paste0('blastn -dust no -soft_masking false -word_size ', wordSize, ' -evalue ', evalue,' -outfmt 6 -query ',
                 file.path(tmpDirPath, paste0(f, '.fasta')), ' -db ',
                 dbPath, ' -num_threads 1 -out ', file.path(tmpDirPath, paste0(f, '.blast')))
  
  system(comm, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  waitForFile(file.path(tmpDirPath, paste0(f, '.blast')))
  
  if(file.info(file.path(tmpDirPath, paste0(f, '.blast')))$size > 0){
    b <- read.table(file.path(tmpDirPath, paste0(f, '.blast')), sep = '\t', header = FALSE)
    names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
    
    o <- reads[names(reads) %in% b$qname]
    d <- tibble(qname = names(o), qlen = width(o))
    b <- left_join(b, d, by = 'qname')
    return(b)
  } else {
    return(tibble())
  }
}



nearestGene <- function(posids, genes, exons, CPUs = 20){
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(parallel))
  
  o <- base::strsplit(posids, '[\\+\\-]')
  d <- tibble(chromosome = unlist(lapply(o, '[', 1)),
              position = unlist(lapply(o, '[', 2)),
              strand = stringr::str_extract(posids, '[\\+\\-]'))
  d$n <- ntile(1:nrow(d), CPUs)
  d$position <- as.integer(d$position)
  
  cluster <- makeCluster(CPUs)
  clusterSetRNGStream(cluster, 1)
  clusterExport(cluster, c('genes', 'exons'), envir = environment())
  
  r <- bind_rows(parLapply(cluster, split(d, d$n), function(x){
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(GenomicRanges))
    
    x$exon <- NA
    x$nearestGene <- NA
    x$nearestGeneStrand <- NA
    x$nearestGeneDist <- NA
    x$inGene <- FALSE
    x$inExon <- FALSE
    x$beforeNearestGene <- NA
    
    r <- GenomicRanges::makeGRangesFromDataFrame(x, 
                                                 seqnames.field = 'chromosome',
                                                 start.field = 'position',
                                                 end.field = 'position')
    
    o <- suppressWarnings(GenomicRanges::findOverlaps(r, exons, select='all', ignore.strand=TRUE, type='any'))
    
    if(length(queryHits(o)) > 0){
      for(i in unique(queryHits(o))){
        x[i,]$exon <- paste0(unique(exons[subjectHits(o[queryHits(o) == i]),]$name2), collapse = ', ')
        x[i,]$nearestGeneStrand <- paste0(unique(as.character(strand(exons[subjectHits(o[queryHits(o) == i])]))), collapse = ',')
      }
      
      if(any(! is.na(x$exon))){
        x[! is.na(x$exon),]$nearestGene <- x[! is.na(x$exon),]$exon
        x[! is.na(x$exon),]$nearestGeneDist <- 0
      }
    }
    
    x1 <- tibble()
    x2 <- tibble()
    
    if(any(! is.na(x$exon))) x1 <- x[! is.na(x$exon),]  # Exon hits
    if(any(is.na(x$exon)))   x2 <- x[is.na(x$exon),]    # Non-exon hits
    
    if(nrow(x2) > 0){
      r <- GenomicRanges::makeGRangesFromDataFrame(x2, 
                                                   seqnames.field = 'chromosome',
                                                   start.field = 'position',
                                                   end.field = 'position')
    
      o <- suppressWarnings(GenomicRanges::distanceToNearest(r, genes, select='all', ignore.strand=TRUE))
      
      if(length(queryHits(o)) > 0){
        for(i in unique(queryHits(o))){
          x2[i,]$nearestGene <- paste0(unique(genes[subjectHits(o[queryHits(o) == i]),]$name2), collapse = ', ')
          x2[i,]$nearestGeneStrand <- paste0(unique(as.character(strand(genes[subjectHits(o[queryHits(o) == i])]))), collapse = ',')
          x2[i,]$nearestGeneDist <- unique(mcols(o[queryHits(o) == i])$distance) # Always single?
          x2[i,]$beforeNearestGene <- ifelse(x2[i,]$position < min(start(genes[subjectHits(o[queryHits(o) == i])])), TRUE, FALSE)
        }
      }
      
      x2 <- x2[! is.na(x2$nearestGene),]
      
      if(any(x2$nearestGeneDist == 0)) x2[which(x2$nearestGeneDist == 0),]$beforeNearestGene <- NA
    }
    
    if(nrow(x1) > 0  & nrow(x2) > 0)  x <- bind_rows(x1, x2)
    if(nrow(x1) == 0 & nrow(x2) > 0)  x <- x2
    if(nrow(x1) > 0  & nrow(x2) == 0) x <- x1
    
    x
  }))
  
  stopCluster(cluster)
  
  if(any(r$nearestGeneDist == 0)) r[which(r$nearestGeneDist == 0),]$inGene <- TRUE
  if(any(! is.na(r$exon))) r[which(! is.na(r$exon)),]$inExon <- TRUE
  if('n' %in% names(r)) r$n <- NULL
  if('exon' %in% names(r)) r$exon <- NULL
  r
}


parseBLAToutput <- function(f){
  if(! file.exists(f) | file.info(f)$size == 0) return(tibble::tibble())
  b <- readr::read_delim(f, delim = '\t', col_names = FALSE, col_types = 'iiiiiiiicciiiciiiiccc')
  
  x <- read.table(textConnection(system(paste(file.path(opt$softwareDir, 'bin', 'pslScore.pl') ,  f), intern = TRUE)), sep = '\t')
  names(x) <- c('tName', 'tStart', 'tEnd', 'hit', 'pslScore', 'percentIdentity')
  
  names(b) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand',
                'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')
  
  b$queryPercentID <- x$percentIdentity
  b$pslScore <- x$pslScore
  b$tStart <- x$tStart
  b$tEnd <- x$tEnd - 1 # Correct for zero-based half open coordinate system.
  
  dplyr::select(b, qName, matches, strand, qSize, qStart, qEnd, tName, tNumInsert, qNumInsert, tBaseInsert, qBaseInsert, tStart, tEnd, queryPercentID, pslScore)
}


blast2rearangements_worker <-  function(b){
  library(IRanges)
  library(dplyr)
  library(data.table)
  
  data.table::rbindlist(lapply(split(b, b$qname), function(b2){
    
    # Alignment to the negative strand will result in the subject end to come before the start.
    # Switch it back so that they are sequential.
    b2 <- arrange(b2, qstart, evalue)
    b2$sstart2 <- ifelse(b2$sstart > b2$send, b2$send, b2$sstart)
    b2$send2   <- ifelse(b2$sstart > b2$send, b2$sstart, b2$send)
    
    
    ir <- IRanges(start = b2$qstart, end = b2$qend)
    if(length(ir) == 0) return(data.table())
    
    r <- IRanges::reduce(ir, min.gapwidth = opt$prepReads_mapLeaderSeqsMaxGapBetweenAlignments)
    r <- r[start(r) <= opt$prepReads_buildReadMaps_minMapStartPostion]
    if(length(r) == 0) return(data.table())
    
    o <- data.table(qname = b2$qname[1], end = IRanges::end(r), width = width(r))
    o <- o[o$width == max(o$width),]
    o[1,]
  }))
}



buildRepLeaderSeqTable <- function(xx, collapseReplicates = FALSE){
  xx$posid2 <- sub('\\.\\d+$', '', xx$posid)
  
  leaderSeqReps <- data.table()

  if(collapseReplicates){
    xx <- group_by(xx, trial, subject, sample, posid2) %>% 
            mutate(g = cur_group_id(), replicate = '*') %>%
            ungroup()
  } else {
    xx <- group_by(xx, trial, subject, sample, posid2, replicate) %>% 
          mutate(g = cur_group_id()) %>%
          ungroup()
  }
  
  o <- split(xx, xx$g)
  
  for(x in o){
    x$w <- x$fragEnd - x$fragStart + 1
  
    if(n_distinct(x$posid) > 1){
      # First remnant
      updateLog(paste0('Fa ', x$g[1]))
      f <- x[grepl('\\.1$', x$posid),] 
    
      repSeq <- group_by(f, repLeaderSeq) %>%
                summarise(frags = n_distinct(w), reads = sum(reads)) %>% 
                ungroup() %>%
                arrange(desc(frags), desc(reads)) %>%
                dplyr::slice(1) %>%
                dplyr::pull(repLeaderSeq)
      
      leaderSeqReps <- distinct(data.table::rbindlist(list(leaderSeqReps, 
                                                           data.table(trial = x$trial[1], subject = x$subject[1], sample = x$sample[1], 
                                                                      posid = f$posid[1], replicate = f$replicate[1], repLeaderSeq = repSeq))))
      
      i <- as.integer(stringr::str_extract(unique(x$posid), '\\d+$'))
      i <- sort(i[i != 1])
    
      invisible(lapply(i, function(n){
        lowerPosID = paste0(x$posid2[1], '.', n-1)
        lowerRepSeq <- subset(leaderSeqReps, trial == x$trial[1] & subject == x$subject[1] & sample == x$sample[1] & posid == lowerPosID)$repLeaderSeq
        theseRepSeqs <- unique(subset(x, posid == paste0(x$posid2[1], '.', n))$repLeaderSeq)
      
        if(length(theseRepSeqs) > 1){
          m <- stringdist::stringdistmatrix(lowerRepSeq, theseRepSeqs)
        
          candidateReps <-  theseRepSeqs[which(m == max(m))]
        
          repSeq <- group_by(subset(x, posid == paste0(x$posid2[1], '.', n) & repLeaderSeq %in% candidateReps), repLeaderSeq) %>%
                    summarise(frags = n_distinct(w), reads = sum(reads)) %>% 
                    ungroup() %>%
                    arrange(desc(frags), desc(reads)) %>%
                    dplyr::slice(1) %>%
                    dplyr::pull(repLeaderSeq)
        } else {
          repSeq <- theseRepSeqs
        }

        leaderSeqReps <- distinct(data.table::rbindlist(list(leaderSeqReps, data.table(trial = x$trial[1], subject = x$subject[1], sample = x$sample[1], 
                                                                                        posid = paste0(x$posid2[1], '.', n), replicate = x$replicate[1], repLeaderSeq = repSeq))))
        
      }))
    } else {
      # Only one posid in this fragment data chunk.
      repSeq <- group_by(x, repLeaderSeq) %>%
                summarise(frags = n_distinct(w), reads = sum(reads)) %>% 
                ungroup() %>%
                arrange(desc(frags), desc(reads)) %>%
                dplyr::slice(1) %>%
                dplyr::pull(repLeaderSeq)
    
      leaderSeqReps <- distinct(data.table::rbindlist(list(leaderSeqReps, data.table(trial = x$trial[1], subject = x$subject[1], sample = x$sample[1], 
                                                                                     posid = x$posid[1], replicate = x$replicate[1], repLeaderSeq = repSeq))))
    }
  }

  as_tibble(leaderSeqReps)
}



ensureDataTypes <- function(opt){
  
 opt$core_CPUs <- as.integer(opt$core_CPUs)
 
 opt$demultiplex_CPUs <- as.integer(opt$demultiplex_CPUs)
 opt$demultiplex_qualtrim_halfWidth <- as.integer(opt$demultiplex_qualtrim_halfWidth)
 opt$demultiplex_qualtrim_events    <- as.integer(opt$demultiplex_qualtrim_events)
 opt$demultiplex_qualtrim_minLength <- as.integer(opt$demultiplex_qualtrim_minLength)
 opt$demultiplex_sequenceChunkSize  <- as.integer(opt$demultiplex_sequenceChunkSize)
 opt$demultiplex_index1ReadMaxMismatch <- as.integer(opt$demultiplex_index1ReadMaxMismatch)
 opt$demultiplex_adriftReadLinkerBarcodeMaxMismatch <- as.integer(opt$demultiplex_adriftReadLinkerBarcodeMaxMismatch)
 opt$demultiplex_anchorReadStartSeq.maxMisMatch <- as.integer(opt$demultiplex_anchorReadStartSeq.maxMisMatch)
 opt$demultiplex_requirePostUmiLinker_maxMismatch <- as.integer(opt$demultiplex_requirePostUmiLinker_maxMismatch)

 opt$prepReads_CPUs <- as.integer(opt$prepReads_CPUs)
 opt$prepReads_vectorAlignmentTestLength <- as.integer(opt$prepReads_vectorAlignmentTestLength)
 opt$prepReads_vectorAlignmentTest_minPercentCoverage <- as.numeric(opt$prepReads_vectorAlignmentTest_minPercentCoverage)
 opt$prepReads_vectorAlignmentTest_minPercentID       <- as.numeric(opt$prepReads_vectorAlignmentTest_minPercentID)
 opt$prepReads_minAnchorReadLength <- as.integer(opt$prepReads_minAnchorReadLength)
 opt$prepReads_minAdriftReadLength <- as.integer(opt$prepReads_minAdriftReadLength)
 opt$prepReads_mapLeaderSeqsMinAlignmentLength <- as.integer(opt$prepReads_mapLeaderSeqsMinAlignmentLength)
 opt$prepReads_mapLeaderSeqsMinPercentID       <- as.numeric(opt$prepReads_mapLeaderSeqsMinPercentID)
 opt$prepReads_mapLeaderSeqsMaxGapBetweenAlignments <- as.integer(opt$prepReads_mapLeaderSeqsMaxGapBetweenAlignments)
 opt$prepReads_buildReadMaps_minMapStartPostion     <- as.integer(opt$prepReads_buildReadMaps_minMapStartPostion)
 opt$prepReads_HMMsearchReadStartPos <- as.integer(opt$prepReads_HMMsearchReadStartPos)
 opt$prepReads_HMMsearchReadEndPos   <- as.integer(opt$prepReads_HMMsearchReadEndPos)
 opt$prepReads_HMMmaxStartPos        <- as.integer(opt$prepReads_HMMmaxStartPos)
 opt$prepReads_HMMminFullBitScore    <- as.numeric(opt$prepReads_HMMminFullBitScore)
 opt$prepReads_HMMmatchEndRadius     <- as.integer(opt$prepReads_HMMmatchEndRadius)
 opt$prepReads_cutAdaptErrorRate     <- as.numeric(opt$prepReads_cutAdaptErrorRate)

 opt$alignReads_CPUs <- as.integer(opt$alignReads_CPUs)
 opt$alignReads_genomeAlignment_anchorRead_maxStartPos <- as.integer(opt$alignReads_genomeAlignment_anchorRead_maxStartPos)
 opt$alignReads_genomeAlignment_adriftRead_maxStartPos <- as.integer(opt$alignReads_genomeAlignment_adriftRead_maxStartPos)
 opt$alignReads_genomeAlignment_minPercentID <- as.numeric(opt$alignReads_genomeAlignment_minPercentID)
 opt$alignReads_genomeAlignment_blatStepSize <- as.integer(opt$alignReads_genomeAlignment_blatStepSize)
 opt$alignReads_genomeAlignment_blatTileSize <- as.integer(opt$alignReads_genomeAlignment_blatTileSize)
 opt$alignReads_genomeAlignment_blatRepMatch <- as.integer(opt$alignReads_genomeAlignment_blatRepMatch) 
 opt$alignReads_genomeAlignment_blatMaxtNumInsert <- as.integer(opt$alignReads_genomeAlignment_blatMaxtNumInsert)
 opt$alignReads_genomeAlignment_blatMaxqNumInsert <- as.integer(opt$alignReads_genomeAlignment_blatMaxqNumInsert)
 opt$alignReads_genomeAlignment_blatMaxtBaseInsert <- as.integer(opt$alignReads_genomeAlignment_blatMaxtBaseInsert)
 opt$alignReads_genomeAlignment_blatMaxqBaseInsert <- as.integer(opt$alignReads_genomeAlignment_blatMaxqBaseInsert)
 opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned <- as.integer(opt$alignReads_genomeAlignment_anchorReadEnd_maxUnaligned)
 opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned <- as.integer(opt$alignReads_genomeAlignment_adriftReadEnd_maxUnaligned)
 
 opt$buildFragments_CPUs <- as.integer(opt$buildFragments_CPUs)
 opt$buildFragments_idGroup_size      <- as.integer(opt$buildFragments_idGroup_size)
 opt$buildFragments_maxReadAlignments <- as.integer(opt$buildFragments_maxReadAlignments)
 opt$buildFragments_minFragLength     <- as.integer(opt$buildFragments_minFragLength)
 opt$buildFragments_maxFragLength     <- as.integer(opt$buildFragments_maxFragLength)
 opt$buildFragments_randomLinkerID_minReadCountToSegreagate       <- as.integer(opt$buildFragments_randomLinkerID_minReadCountToSegreagate)
 opt$buildFragments_randomLinkerID_minSingleSampleMajorityPercent <- as.numeric(opt$buildFragments_randomLinkerID_minSingleSampleMajorityPercent)
 
 opt$buildStdFragments_CPUs <- as.integer(opt$buildStdFragments_CPUs)
 opt$buildStdFragments_randomIDdupReadMult <- as.numeric(opt$buildStdFragments_randomIDdupReadMult)
 opt$buildStdFragments_randomIDdupAbundMult <- as.numeric(opt$buildStdFragments_randomIDdupAbundMult)
 opt$buildStdFragments_intSiteStdWindowWidth <- as.integer(opt$buildStdFragments_intSiteStdWindowWidth)
 opt$buildStdFragments_breakPointStdWindowWidth <- as.integer(opt$buildStdFragments_breakPointStdWindowWidth)
 opt$buildStdFragments_minReadsPerFrag <- as.integer(opt$buildStdFragments_minReadsPerFrag)
 opt$buildStdFragments_maxFragLength <- as.integer(opt$buildStdFragments_maxFragLength)
 opt$buildStdFragments_fragEvalAdriftReadTestLen <- as.integer(opt$buildStdFragments_fragEvalAdriftReadTestLen)
 opt$buildStdFragments_UMIminPercentReads <- as.numeric(opt$buildStdFragments_UMIminPercentReads) 
 opt$buildStdFragments_fragEvalAnchorReadTestLen <- as.integer(opt$buildStdFragments_fragEvalAnchorReadTestLen)
 opt$buildStdFragments_fragEvalAnchorReadMinAbundDiff <- as.integer(opt$buildStdFragments_fragEvalAnchorReadMinAbundDiff)
 opt$buildStdFragments_fragEvalAnchorReadMinReadMult <- as.numeric(opt$buildStdFragments_fragEvalAnchorReadMinReadMult)
 
 opt$buildSites_CPUs <- as.integer(opt$buildSites_CPUs)
 opt$buildSites_dualDetectWidth <- as.integer(opt$buildSites_dualDetectWidth)
 opt$buildSites_integraseCorrectionDist <- as.integer(opt$buildSites_integraseCorrectionDist)
 
 opt
}