library(lubridate)
library(dplyr)
library(openxlsx)
options(stringsAsFactors = FALSE)

# Read in configuration file.
configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)

# 3595
sites <- readRDS(file.path(opt$outputDir, opt$addSampleDetails_inputFile))
details <- readr::read_tsv(opt$addSampleDetails_detailsFile)

sites$k <- paste(sites$trial, sites$subject, sites$sample)
details$k <- paste(details$trial, details$subject, details$sample)

sites <- left_join(sites, select(details, -trial, -subject, -sample), by = 'k') %>% select(-k)

colNames <- names(select(details, -trial, -subject, -sample, -k))
colsNames2 <- c('nearestGene', 'nearestGeneStrand', 'nearestGeneDist', 'inGene', 'inExon',	'beforeNearestGene')

sites <- relocate(sites, colNames, .after = opt$addSampleDetails_addAfter)
sites <- relocate(sites, colsNames2, .after = nRepsObs)


a <- subset(sites, nearestGene != 'F8')
a <- rename(a, UMIs = fragments)
a <- rename(a, sonicLengths = fragmentWidths)
i <- which(is.na(a$PCRartifact1) & is.na(a$PCRartifact2) & a$ reads >= 3)

a1 <- a[i,]
a2 <- a[-i,]

wb <- createWorkbook()
addWorksheet(wb, "Passing filters")
addWorksheet(wb, "Failing filters")

writeData(wb, "Passing filters", arrange(a1, desc(sample), desc(sonicLengths)), startRow = 1, startCol = 1)
writeData(wb, "Failing filters", arrange(a2, desc(sample), desc(sonicLengths)), startRow = 1, startCol = 1)

saveWorkbook(wb, file = "sites.xlsx", overwrite = TRUE)



