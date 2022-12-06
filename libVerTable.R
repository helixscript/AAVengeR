path <- './'
f <- unique(unlist(lapply(system(paste0('find ', path, ' -name "*.R"'), intern = TRUE), function(x){
       r <- suppressWarnings(readLines(x))
       i <- grepl('^\\s*library\\(', r, ignore.case = TRUE)
       if(any(i)) return(unique(gsub('\\s*library\\((.*)\\).*', '\\1', r[i], perl = TRUE)))
      })))

ip <- data.frame(installed.packages()[,c(1,3)])
ip <- ip[ip$Package %in% f,]
rownames(ip) <- NULL
ip