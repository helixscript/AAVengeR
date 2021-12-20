library(dplyr)
library(lubridate)

configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)

source(file.path(opt$softwareDir, 'lib.R'))

dir.create(file.path(opt$outputDir, opt$processHIVsites_outputDir))

sites <- readRDS(file.path(opt$outputDir, opt$processHIVsites_inputFile))
if(! 'flags' %in% names(sites)){
  write(c(paste(now(), 'Error -- the flags column is not in the sites data object.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}
  
if(! all(grepl('HIV_u5|HIV_u3', sites$flags))){
  write(c(paste(now(), 'Error -- all sites do not have either a HIV_u5 or HIV_u3 flag.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
  q(save = 'no', status = 1, runLast = FALSE) 
}

if(file.exists(file.path(opt$outputDir, opt$processHIVsites_outputDir, 'dual_detections.tsv'))) file.remove(file.path(opt$outputDir, opt$processHIVsites_outputDir, 'dual_detections.tsv'))

# Adjustments
# + strand alignment: +2
# - stramd alignment: -2
# flip strands for u3 experiments. Aligning to the negative strand in u3 means a positive insertion orientation.

sites.org <- sites
sites$n <- 1:nrow(sites)
a <- subset(sites, strand == '+')
b <- subset(sites, strand == '-')

a$position <- a$position + 2
b$position <- b$position - 2

sites <- bind_rows(a, b)
sites$posid <- paste0(sites$chromosome, sites$strand, sites$position)


#chr14+69415488

sites <- bind_rows(lapply(split(sites, sites$subject), function(x){
  merged_sites <- vector()
  new_sites <- vector()

  sites2 <- bind_rows(lapply(1:nrow(x), function(y){
    y <- x[y,]
    
    x <- subset(x, chromosome == y$chromosome)
    o <- subset(x, position >= (y$position - opt$processHIVsites_dualDetectDistance) & position <= (y$position + opt$processHIVsites_dualDetectDistance))
    
    o <- subset(o, ! posid %in% merged_sites)
    if(nrow(o) == 0) return(tibble())
    
    if(nrow(o) == 1) return(y) # No other sites nearby to compare with.
    
    if(length(unique(o$strand)) == 1) return(y)  # Want a mix of both + and - strand alignemnts.

    if(!(any(grepl('HIV_u5', o$flags)) & any(grepl('HIV_u3', o$flags)))) return(return(y)) # Want a mix of u5 and u3 flags.

    if(nrow(o) == 2){
          merged_sites <<- append(merged_sites, c(o[1,]$posid, o[2,]$posid))
          o <- arrange(o, position)
          y$position <- o[1,]$position + floor((o[2,]$position - o[1,]$position)/2)
          y$strand <- '*'
          y$reads <- sum(o$reads)
          y$posid <- paste0(y$chromosome, y$strand, y$position)
          new_sites <<- append(new_sites, y$posid)
          y$flags <- 'U3,U5 merged'
          
          if(opt$processHIVsites_dualDetectAubundanceMethod == 'max'){
            y$estAbund <- max(o$estAbund)
          } else if(opt$processHIVsites_dualDetectAubundanceMethod == 'average'){
            y$estAbund <- floor(mean(o$estAbund))
          } else {
            write(c(paste(now(), 'Error - unknown value for processHIVsites_dualDetectAubundanceMethod')), file = file.path(opt$outputDir, 'log'), append = TRUE)
            q(save = 'no', status = 1, runLast = FALSE) 
          }
          
          return(y)
    } else {
      
          o <- o[order(abs(mean(o$position) - o$position)),]  # order sites by distance from the center of the distance interval.
          tag <- stringr::str_extract(o$flags, 'HIV_u\\d')    # extract the HIV tags.
          o$s <- paste(o$strand, tag)                         # create new tag: strand + tag, eg.  + HIV_u5
          i <- which(! o$s[2:nrow(o)] == o$s[1])[1]+1         # determine the position of the next closest site with opposite strand and tag.
          o1 <- o[c(1, i),]
          
          merged_sites <<- append(merged_sites, c(o1[1,]$posid, o1[2,]$posid))
          
          y$position <- o1[1,]$position + floor((o1[2,]$position - o1[1,]$position)/2)
          y$strand <- '*'
          y$reads <- sum(o1[1,]$reads, o1[2,]$reads)
          y$posid <- paste0(y$chromosome, y$strand, y$position)
          
          new_sites <<- append(new_sites, y$posid)
          
          y$flags <- 'U3,U5 merged'
          
          if(opt$processHIVsites_dualDetectAubundanceMethod == 'max'){
            y$estAbund <- max(o1[1,]$estAbund, o1[2,]$estAbund)
          } else if(opt$processHIVsites_dualDetectAubundanceMethod == 'average'){
            y$estAbund <- floor(mean(o1[1,]$estAbund, o1[2,]$estAbund))
          } else {
            write(c(paste(now(), 'Error - unknown value for processHIVsites_dualDetectAubundanceMethod.')), file = file.path(opt$outputDir, 'log'), append = TRUE)
            q(save = 'no', status = 1, runLast = FALSE)
          }
          
          return(y)
        }
    }))
  
  if(length(new_sites) > 0){
    new_sites <- unlist(lapply(new_sites, rep, 2))
    a <- split( merged_sites, new_sites)
    b <- bind_rows(mapply(function(s, n){ data.frame(site1 = s[1], site2 = s[2], merged_site = n)}, a, names(a), SIMPLIFY = FALSE))
    b$subject <- x$subject[1]
    readr::write_tsv(b, file.path(opt$outputDir, opt$processHIVsites_outputDir, 'dual_detections.tsv'), append = TRUE)
  }
  
  sites2
}))


# Flip the strand for u3 detections to reflect orientaion rather than alignment strand.
a <- sites[grepl('HIV_u3', sites$flags, ignore.case = TRUE),]
b <- sites[! grepl('HIV_u3', sites$flags, ignore.case = TRUE),]

a$strand <- ifelse(a$strand == '+', '-', '+')
sites <- bind_rows(a, b)
sites$posid <- paste0(sites$chromosome, sites$strand, sites$position)
sites$n <- NULL

saveRDS(sites, file = file.path(opt$outputDir, opt$processHIVsites_outputDir, opt$processHIVsites_outputFile))

q(save = 'no', status = 0, runLast = FALSE) 
