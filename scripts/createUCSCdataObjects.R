library(GenomicRanges)
library(tidyverse)
options(stringsAsFactors = FALSE)

humanGeneFilter <- function(d){
  humanGenes <- gsub('\\.\\d+$', '', gt23::hg38.refSeqGenesGRanges$name)
  d$name <- gsub('\\.\\d+$', '', d$name)
  subset(d, toupper(name) %in% toupper(humanGenes))
}

createRefSeqObjects <- function(genomeLabel, file, humanGeneFilter = FALSE){
  d <- read.table(file, sep = '\t', header = FALSE, quote = '')
  
  names(d) <- c('bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart',
                'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2',
                'cdsStartStat', 'cdsEndStat', 'exonFrames')
  
  if(humanGeneFilter) d <- humanGeneFilter(d)
  
  g <- makeGRangesFromDataFrame(d, seqnames.field = 'chrom', start.field = 'txStart',
                                end.field = 'txEnd', strand.field = 'strand',
                                keep.extra.columns = TRUE)
  
  d$exonStarts2 <- strsplit(d$exonStarts, ',')
  d$exonEnds2 <- strsplit(d$exonEnds, ',')
  
  e <- group_by(d, 1:nrow(d)) %>%
    unnest(cols = c(exonStarts2, exonEnds2)) %>%
    mutate(name = paste('exon', 1:n()),
           exonStarts = as.integer(exonStarts2),
           exonEnds = as.integer(exonEnds2)) %>%
    ungroup() %>%
    makeGRangesFromDataFrame(seqnames.field = 'chrom', start.field = 'exonStarts2',
                             end.field = 'exonEnds2', strand.field = 'strand',
                             keep.extra.columns = TRUE)
  
  saveRDS(g, file = paste0('../data/genomeAnnotations/', genomeLabel, '.TUs.rds'))
  saveRDS(e, file = paste0('../data/genomeAnnotations/', genomeLabel, '.exons.rds'))
}


createRefSeqObjects('hg38',    'geneReferences/June2019/hg38.refSeq.curated.txt.gz')
createRefSeqObjects('hg18',    'geneReferences/June2019/hg18.refGene.txt.gz')
createRefSeqObjects('mm9',     'geneReferences/June2019/mm9.refGene.txt.gz')
createRefSeqObjects('susScr3', 'geneReferences/June2019/susScr3.refGene.txt.gz')
createRefSeqObjects('macFas5', 'geneReferences/June2019/macFas5.refGene.txt.gz')
createRefSeqObjects('canFam3', 'geneReferences/June2019/canFam3.refGene.txt.gz')
createRefSeqObjects('canFam3.humanXeno', 'geneReferences/June2019/canFam3.xenoRefGene.txt.gz', humanGeneFilter = TRUE)
createRefSeqObjects('macFas5.humanXeno', 'geneReferences/June2019/macFas5.xenoRefGene.txt.gz', humanGeneFilter = TRUE)

createRefSeqObjects('sacCer3', 'geneReferences/Sept2021/sacCer3.refGene.curated.txt.gz')
