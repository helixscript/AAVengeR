#!/usr/bin/Rscript
o <- readLines(commandArgs(trailingOnly=TRUE))
a <- stringr::str_extract(o, '[a-z]\\s+\\-\\s\\-')
a <- a[!is.na(a)]
a <- paste0(toupper(sub('\\s\\-\\s\\-', '', a)), collapse = '')
message(a)
