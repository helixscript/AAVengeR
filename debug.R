
o <- readLines('~/chr22-1840213.reads')


a <- readFastq('/home/everett/projects/AAVengeR/data/testData/R2.test.fastq.gz')
a.ids <- sub('\\s.+', '', as.character(a@id))
table(o %in% a.ids)  # all true

# where is M03249:365:000000000-C3CH4:1:2110:8655:4587?

r0 <- Reduce('append', lapply(list.files('out/demultiplex/', pattern = 'anchor', full.names = TRUE), readFasta))


r1 <- Reduce('append', lapply(list.files('out/uniqueFasta/', pattern = 'anchor', full.names = TRUE), readFasta))
table(o %in% as.character(r1@id))

missing <- o[! o %in% as.character(r1@id)]
length(missing) # 4091 missing


d <- readr::read_tsv('out/uniqueFasta/dupReadPairMap.tsv')
d_id_list <- unlist(strsplit(unlist(d$id_list), ','))
table(missing %in% d_id_list)
# FALSE  TRUE 
#  488  4413

still_missing <- missing[! missing %in% d_id_list]





o <- subset(o, ! o %in% d$id)