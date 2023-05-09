o <- readLines('HIV1_1-100_U3_RC.hmm')
a <- stringr::str_extract(o, '[a-z]\\s+\\-\\s\\-')
a <- a[!is.na(a)]
a <- paste0(toupper(sub('\\s\\-\\s\\-', '', a)), collapse = '')