library(dplyr)
library(ShortRead)
library(GenomicRanges)
source('~/AAVengeR/stdPos.lib.R')

chromosome      <- 'chr1'
position        <- 13736069
leaderSeqLen    <- 62            # For HIV-1 HMMs -- u5: 69  u3: 62
linkerSeqLen    <- 48            # typically 48
#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

radius          <- 1000          # Radius around target position for anchor alignments
radius2         <- radius * 10   # Radius around target for subsequent adrift alignments
minMatches      <- 25
minReadLength   <- 30
minPercentSeqID <- 97
genome          <- '../../../AAVengeR/data/referenceGenomes/blat/hg38.2bit'
anchorReads     <- '../merged/data/Undetermined_S0_R2_001.fastq.gz' 
adriftReads     <- '../merged/data/Undetermined_S0_R1_001.fastq.gz' 
blat            <- '../../../software/blat'
faToTwoBit      <- '../../../software/faToTwoBit'


# Import genome sequence.
g <- rtracklayer::import.2bit(genome)
g <- g[names(g) == chromosome]

# Write out region around genome target.
t <- Biostrings::subseq(g, position - radius, position + radius)
Biostrings::writeXStringSet(t, 'tmpSeqFile')
system(paste0(faToTwoBit, ' tmpSeqFile tmpSeqFile.2bit'))
file.remove('tmpSeqFile')
rm(t)

# Convert fastq to fasta with shorter read ids.
r <- ShortRead::readFastq(anchorReads)
r <- ShortRead::trimTailw(r, 2, '0', 5)
r <- r[width(r) >= leaderSeqLen + 2]
r <- narrow(r, leaderSeqLen + 1, width(r))
r <- r[width(r) >= minReadLength]

# Write out anchor reads for alignment.
r@id <- BStringSet(sub('\\s.+$', '', as.character(r@id)))
writeFasta(r, 'anchorReads.fasta')

# Align all anchor reads to short reference.
system(paste0(blat, ' tmpSeqFile.2bit anchorReads.fasta anchorReads.psl ',
              ' -tileSize=11 -stepSize=9 -repMatch=3000 -out=psl -t=dna -q=dna -minScore=0 -minIdentity=0 -noHead -noTrimA'))

file.remove(c('tmpSeqFile.2bit', 'anchorReads.fasta'))

# Read in anchor alignments and filter.
b1 <- readr::read_delim('anchorReads.psl', delim = '\t', col_names = FALSE, col_types = readr::cols())

names(b1) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 
              'tNumInsert', 'tBaseInsert', 'anchor.strand','qName', 'qSize', 'qStart', 'qEnd', 'tName', 
              'tSize', 'anchor.tStart', 'anchor.tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')

b1$anchor.tEnd <- b1$anchor.tEnd - 1  # Correct for half-open coord system.
b1$queryPercentID       <- (b1$matches / b1$qSize)*100
b1$tAlignmentWidth      <- (b1$anchor.tEnd - b1$anchor.tStart) + 1
b1$queryWidth           <- (b1$qEnd - b1$qStart) + 1
b1$alignmentPercentID   <- (b1$matches / b1$tAlignmentWidth)*100

b1 <- dplyr::filter(b1, qStart <= 1, matches >= minMatches, alignmentPercentID >= minPercentSeqID,
                    tNumInsert <= 1, qNumInsert <= 1, tBaseInsert <= 2, qBaseInsert <= 2) %>%
      select(qName, anchor.strand, anchor.tStart, anchor.tEnd)

b1$anchor.tStart <- b1$anchor.tStart + (position - radius)
b1$anchor.tEnd <-   b1$anchor.tEnd + (position - radius)

file.remove('anchorReads.psl') 
 
# Convert adrift fastq to fasta with shorted ids.
r <- ShortRead::readFastq(adriftReads)
r <- ShortRead::trimTailw(r, 2, '0', 5)
r <- r[width(r) >= linkerSeqLen + 2]
r <- narrow(r, linkerSeqLen + 1, width(r))
r <- r[width(r) >= minReadLength]

# Write out anchor reads for alignment.
r@id <- BStringSet(sub('\\s.+$', '', as.character(r@id)))
writeFasta(r, 'adriftReads.fasta')

# Write out region around genome target.
t <- Biostrings::subseq(g, position - radius2, position + radius2)
Biostrings::writeXStringSet(t, 'tmpSeqFile')
system(paste0(faToTwoBit, ' tmpSeqFile tmpSeqFile.2bit'))
file.remove('tmpSeqFile')
rm(t)

# Align all anchor reads to short reference.
system(paste0(blat, ' tmpSeqFile.2bit adriftReads.fasta adriftReads.psl ',
              ' -tileSize=11 -stepSize=9 -repMatch=3000 -out=psl -t=dna -q=dna -minScore=0 -minIdentity=0 -noHead -noTrimA'))

file.remove(c('tmpSeqFile.2bit', 'adriftReads.fasta'))

# Read in anchor alignments and filter.
b2 <- readr::read_delim('adriftReads.psl', delim = '\t', col_names = FALSE, col_types = readr::cols())

names(b2) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 
               'tNumInsert', 'tBaseInsert', 'adrift.strand','qName', 'qSize', 'qStart', 'qEnd', 'tName', 
               'tSize', 'adrift.tStart', 'adrift.tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')

b2$adrift.tEnd <- b2$adrift.tEnd - 1  # Correct for half-open coord system.
b2$queryPercentID       <- (b2$matches / b2$qSize)*100
b2$tAlignmentWidth      <- (b2$adrift.tEnd - b2$adrift.tStart) + 1
b2$queryWidth           <- (b2$qEnd - b2$qStart) + 1
b2$alignmentPercentID   <- (b2$matches / b2$tAlignmentWidth)*100

b2 <- dplyr::filter(b2, qStart <= 1, matches >= minMatches, alignmentPercentID >= minPercentSeqID,
                    tNumInsert <= 1, qNumInsert <= 1, tBaseInsert <= 2, qBaseInsert <= 2) %>%
      select(qName, adrift.strand, adrift.tStart, adrift.tEnd)

b2$adrift.tStart <- b2$adrift.tStart + (position - radius2)
b2$adrift.tEnd   <- b2$adrift.tEnd + (position - radius2)

file.remove('adriftReads.psl') 

f <- tidyr::drop_na(left_join(b1, b2, by = 'qName'))

f <- mutate(f[! f$anchor.strand == f$adrift.strand,], 
            start  = ifelse(anchor.strand == '+', anchor.tStart + 1, adrift.tStart + 1),
            end    = ifelse(anchor.strand == '+', adrift.tEnd + 1,   anchor.tEnd + 1),
            strand     = ifelse(anchor.strand == '+', '+', '-'),
            chromosome = chromosome,  
            fragTest  = ifelse(anchor.strand == '+', anchor.tStart < adrift.tEnd, adrift.tStart < anchor.tEnd)) %>%
     dplyr::filter(fragTest == TRUE) %>%
     dplyr::select(chromosome, strand, start, end) %>%
     group_by(chromosome, strand, start, end) %>% summarise(reads = n()) %>% ungroup() %>% arrange(end)

k <-  standardize_sites(makeGRangesFromDataFrame(f, keep.extra.columns = TRUE), counts.col = 'reads', sata.gap = 5)

# update positions.
f2 <- group_by(data.frame(k), seqnames, start, end, strand) %>%
      summarise(reads = sum(reads)) %>%
      ungroup()

if(as.character(f2[1,]$strand) == '+'){
  f2 <- arrange(f2, end)
} else {
  f2 <- arrange(f2, start)
}

k2 <- makeGRangesFromDataFrame(f2, keep.extra.columns = TRUE)


r1 <- refine_breakpoints(k2, counts.col = 'reads', sata.gap = 1)
if(as.character(f2[1,]$strand) == '+'){
  message('Refinement with data.gap = 1 num breaks: ', n_distinct(end(r1)))
} else {
  message('Refinement with data.gap = 1 num breaks: ', n_distinct(start(r1)))
}

r2 <- refine_breakpoints(k2, counts.col = 'reads', sata.gap = 2)
if(as.character(f2[1,]$strand) == '+'){
  message('Refinement with data.gap = 2 num breaks: ', n_distinct(end(r2)))
} else {
  message('Refinement with data.gap = 2 num breaks: ', n_distinct(start(r2)))
}

r3 <- refine_breakpoints(k2, counts.col = 'reads', sata.gap = 3)
if(as.character(f2[1,]$strand) == '+'){
  message('Refinement with data.gap = 3 num breaks: ', n_distinct(end(r3)))
} else {
  message('Refinement with data.gap = 3 num breaks: ', n_distinct(start(r3)))
}

f2


