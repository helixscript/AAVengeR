
.o <- readLines('~/chr22-1840213.reads')
.d <- readr::read_tsv(file.path(opt$outputDir, 'uniqueFasta', 'dupReadPairMap.tsv'))
.d_reads <- unlist(strsplit(.d$id_list, ','))
.o <- .o[! .o %in% .d_reads]
sprintf("%.2f%%", sum(names(anchorReads) %in% .o)/length(.o)*100)
sprintf("%.2f%%", sum(anchorReadAlignments$qName %in% .o)/length(.o)*100)
