library(dplyr)
library(GenomicRanges)
library(ggplot2)

comms <- list(
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed1_0error  --mode AAV --nSites 1000 --seed 1  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed2_0error  --mode AAV --nSites 1000 --seed 2  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed3_0error  --mode AAV --nSites 1000 --seed 3  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed4_0error  --mode AAV --nSites 1000 --seed 4  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed5_0error  --mode AAV --nSites 1000 --seed 5  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed6_0error  --mode AAV --nSites 1000 --seed 6  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed7_0error  --mode AAV --nSites 1000 --seed 7  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed8_0error  --mode AAV --nSites 1000 --seed 8  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed9_0error  --mode AAV --nSites 1000 --seed 9  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed10_0error --mode AAV --nSites 1000 --seed 10 --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0 --threads 20',
              
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed1_1error  --mode AAV --nSites 1000 --seed 1  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.01 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed2_1error  --mode AAV --nSites 1000 --seed 2  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.01 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed3_1error  --mode AAV --nSites 1000 --seed 3  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.01 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed4_1error  --mode AAV --nSites 1000 --seed 4  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.01 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed5_1error  --mode AAV --nSites 1000 --seed 5  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.01 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed6_1error  --mode AAV --nSites 1000 --seed 6  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.01 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed7_1error  --mode AAV --nSites 1000 --seed 7  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.01 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed8_1error  --mode AAV --nSites 1000 --seed 8  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.01 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed9_1error  --mode AAV --nSites 1000 --seed 9  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.01 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed10_1error --mode AAV --nSites 1000 --seed 10 --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.01 --threads 20',
              
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed1_3error  --mode AAV --nSites 1000 --seed 1  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.03 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed2_3error  --mode AAV --nSites 1000 --seed 2  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.03 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed3_3error  --mode AAV --nSites 1000 --seed 3  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.03 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed4_3error  --mode AAV --nSites 1000 --seed 4  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.03 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed5_3error  --mode AAV --nSites 1000 --seed 5  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.03 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed6_3error  --mode AAV --nSites 1000 --seed 6  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.03 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed7_3error  --mode AAV --nSites 1000 --seed 7  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.03 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed8_3error  --mode AAV --nSites 1000 --seed 8  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.03 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed9_3error  --mode AAV --nSites 1000 --seed 9  --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.03 --threads 20',
              './buildAndTestSynData.py --outputDir /data/AAVengeR/tests/synData/hs1_AAV_1000sites_seed10_3error --mode AAV --nSites 1000 --seed 10 --refGenomePath /data/AAVengeR/data/referenceGenomes/blat/hs1.2bit --refGenomeID hs1 --R1_length 200 --R2_length 300 --percentGenomicError 0.03 --threads 20')

invisible(sapply(comms, function(x){ message(x); system(x) }))

d <- list.dirs(recursive = FALSE)
d <- d[grepl('hs1', d)]
o <- tibble(dir = d,
            seed = as.integer(sub('seed', '', stringr::str_extract(dir, 'seed\\d+'))),
            error = as.integer(sub('error', '', stringr::str_extract(dir, '\\d+error'))))

r <- readr::read_tsv('../../data/genomeAnnotations/hs1.repeatTable.gz', show_col_types = FALSE)
r$strand <- '*'
r <- makeGRangesFromDataFrame(r, seqnames.field = 'query_seq', start.field = 'query_start', end.field = 'query_end', keep.extra.columns = TRUE)
                              



d <- bind_rows(lapply(split(o, 1:nrow(o)), function(x){
       message(x)
       if(! file.exists(file.path(x$dir, 'eval', 'table2.tsv'))) return(tibble())
       tab1 <- readr::read_tsv(file.path(x$dir, 'eval', 'table1.tsv'), show_col_types = FALSE)
       tab2 <- readr::read_tsv(file.path(x$dir, 'eval', 'table2.tsv'), show_col_types = FALSE)
  
       x$percentUniqueRecovery <- as.numeric(sub('%', '', tab1$percentUniqueRecovery))
       x$percentTotalRecovery  <- as.numeric(sub('%', '', tab1$percentTotalRecovery ))
  
       truth <- readr::read_tsv(file.path(x$dir, 'truth.tsv'), show_col_types = FALSE)
  
       g <- makeGRangesFromDataFrame(tibble(seqnames = unlist(lapply(str_split(tab2$posid, '[\\+\\-]'), '[[', 1)),
                                            strand = '*',
                                            start = unlist(lapply(str_split(tab2$posid, '[\\+\\-]'), '[[', 2)),
                                            end = start))
   
       o <- findOverlaps(g, r, ignore.strand = TRUE)
       tab2$repeatClass <- 'NA'

       invisible(sapply(queryHits(o), function(qi){
         tab2[qi,]$repeatClass <<- names(sort(table( r[subjectHits(o[queryHits(o) == qi])]$repeat_class), decreasing = TRUE))[1]  
       }))
  
       if(file.exists(file.path(x$dir, 'sitesNotFound.tsv'))) file.remove(file.path(x$dir, 'sitesNotFound.tsv'))
       readr::write_tsv(tab2[tab2$found == FALSE,], file.path(x$dir, 'sitesNotFound.tsv'))
  
       x
     }))

d$percentUniqueRecovery <- d$percentUniqueRecovery / 100
d$percentTotalRecovery <- d$percentTotalRecovery / 100

d2 <- bind_rows(lapply(split(d, 1:nrow(d)), function(x){
        tibble(seed = x$seed, 
               cond = c(paste0('Unique recovery\n', x$error, '% error'), paste0('Total recovery\n', x$error, '% error')),
               value = c(x$percentUniqueRecovery, x$percentTotalRecovery))
       }))

d2$cond <- factor(d2$cond, levels = c("Total recovery\n0% error", "Unique recovery\n0% error",
                                      "Total recovery\n1% error", "Unique recovery\n1% error",
                                      "Total recovery\n3% error", "Unique recovery\n3% error"))


set.seed(2)
p1 <- ggplot(d2, aes(cond, value)) +
      theme_bw() +
      scale_y_continuous(limits = c(0.85, 1), labels = scales::percent, expand = c(0, 0)) +
      geom_boxplot(outliers = FALSE)+
      geom_jitter(width = 0.25, height = 0, shape = 21, color = 'black', fill = 'white', size = 3) +
      labs(x = '', y = 'Percent recovery') +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
            axis.title.y = element_text(margin = margin(r = 15)),
            plot.margin = unit(c(1, 1, 1, 1), "cm"))

ggsave('p1.pdf', p1, dpi = 300)


nf <- bind_rows(lapply(list.files(pattern = 'sitesNotFound.tsv', recursive = TRUE), readr::read_tsv, show_col_types = FALSE))
nf[is.na(nf$repeatClass),]$repeatClass <- 'None'

otherClasses <- c('DNA/hAT-Charlie', 'DNA/hAT-Tip100', 'DNA/TcMar-Mariner', 'DNA/TcMar-Tc2', 'LTR', 
                  'LTR/ERVL', 'LTR/Gypsy', 'Satellite/acro', 'srpRNA', 'Unknown', 'LINE/CR1')

nf$repeatClass <- ifelse(nf$repeatClass %in% otherClasses, 'Other', nf$repeatClass)
nf[nf$repeatClass == 'Satellite',]$repeatClass <- 'Satellite/centr'

d <- data.frame(table(nf$repeatClass)) %>% arrange(desc(Freq))
d$p <- d$Freq / sum(d$Freq)
d$Var1 <- factor(d$Var1, levels = d$Var1)
d$x <-  'x'
d$x <- factor(d$x)

fillColors <- colorRampPalette(brewer.pal(12, "Paired"))(nrow(d))

p3 <- ggplot(d, aes(x, p, fill = Var1)) +
  theme_bw() +
  scale_fill_manual(name = 'Repeat class', values = fillColors) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent, expand = c(0, 0)) +
  geom_col(color = 'black', width = 1.2, height = 1.2) + 
  labs(x = '', y = 'Percent unrecovered sites') +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))


ggsave('p3.pdf', p3, dpi = 300)
