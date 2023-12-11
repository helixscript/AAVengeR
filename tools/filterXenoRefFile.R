
refGene_filter_file <- 'hg38.refGene.txt.gz'
xenoRef_file <- 'canFam4.xenoRefGene.txt.gz'
output_file <- 'canFam4.xenoRefGene.hg38Filtered.txt.gz'

h <- readr::read_tsv(refGene_filter_file, col_names = FALSE)
names(h)  <- c('bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart',
               'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2',
               'cdsStartStat', 'cdsEndStat', 'exonFrames')

x <- readr::read_tsv(xenoRef_file, col_names = FALSE)
names(x)  <- c('bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart',
               'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2',
               'cdsStartStat', 'cdsEndStat', 'exonFrames')

x <- subset(x, name %in% h$name)

readr::write_tsv(x, output_file)