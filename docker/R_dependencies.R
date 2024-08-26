
install.packages('remotes',repos='https://cloud.r-project.org/')

library(remotes)

install.packages('devtools',repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("lubridate", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("data.table", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("ggplot2", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("scales", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("tidyverse", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("RColorBrewer", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("gridExtra", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("png", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("RMySQL", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("readr", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages('BiocManager', repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages('RcppAlgos', repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages('stringdist', repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages('openxlsx', repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages('pander', repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages('dtplyr', repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages('argparse', repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages('RMariaDB', repos='https://cloud.r-project.org/', quiet=FALSE)

BiocManager::install(c("GenomicRanges", "ShortRead",'GenomicRanges','Biostrings' , 'rtracklayer'), update = FALSE)

install.packages("igraph", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("this.path", repos='https://cloud.r-project.org/', quiet=FALSE)

install.packages("optparse", repos='https://cloud.r-project.org/', quiet=FALSE)