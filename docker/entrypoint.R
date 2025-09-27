#!/usr/bin/Rscript
f <- Sys.getenv("CONFIG_PATH")
o <- yaml::read_yaml(f)
comm <- paste0('Rscript ', o$softwareDir, '/aavenger.R ', f)
message(comm)
system(comm)
q()
