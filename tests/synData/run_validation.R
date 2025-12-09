library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(RColorBrewer)


CPUs          <- 38
log           <- 'run_validation.log'
mode          <- 'AAV'
outputDir     <- '/data/AAVengeR/tests/synData/output'
refGenome     <- 'hs1'
aavengerDir   <- '/data/AAVengeR'
nSites        <- 500

seeds         <- 1:5
percentErrors <- c(0, 0.01, 0.02, 0.03)

# seeds <- 1
# percentErrors <- 0

R1_length     <- 200
R2_length     <- 300


for (seed in seeds){
  for(percentError in percentErrors){
    thisOutputDir <- paste0(outputDir, '/', mode, '_', refGenome, '_', seed, '_', percentError, '_', nSites)
  
    if(! file.exists(paste0(thisOutputDir, '/eval/table1.tsv'))){
      comm <- paste0('./buildAndTestSynData.py --outputDir ', thisOutputDir, ' --mode ', mode, ' --nSite ', nSites, 
                    ' --seed ', seed, ' --refGenomeID ',  refGenome, ' --refGenomePath ', paste0(aavengerDir, '/data/referenceGenomes/blat/', refGenome, '.2bit'), 
                    ' --R1_length ', R1_length, ' --R2_length ', R2_length, ' --percentGenomicError ',  percentError, ' --threads ', CPUs)
     write(paste0(date(), '\t', comm, '\n'), file = log, append = TRUE)
     system(comm)
    } else {
      message('Skipping ', thisOutputDir)
    }
  }
}

stop()
q()



d <- list.dirs('output', recursive = FALSE)
d <- d[grepl('hs1', d)]

o <- bind_rows(lapply(split(d, 1:length(d)), function(x){
       o <- unlist(strsplit(x, '_'))
       tibble(dir = x,
              seed = as.integer(o[3]),
              error = as.numeric(o[4]))
      }))

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


d$error <- sprintf("%.1f%%", (d$error*100))
d$error <- factor(d$error, levels = c("0.0%", "1.0%", "2.0%", "3.0%"))

d2 <- bind_rows(lapply(split(d, 1:nrow(d)), function(x){
        bind_rows(tibble(error = x$error, group = paste0('seed ', x$seed, ' unique'), group2 = paste0(x$error, '_unique'), value = x$percentUniqueRecovery),
                  tibble(error = x$error, group = paste0('seed ', x$seed, ' total'),  group2 = paste0(x$error, '_total'), value = x$percentTotalRecovery))
      }))


set.seed(1)
pd <- position_dodge(width = 0.5)

p1 <- ggplot(d2[grepl('total', d2$group),], aes(error, value, fill = group)) +
  theme_bw() +
  scale_fill_manual(name = 'Simulations', values = brewer.pal(n=5, "Set1")) + 
  geom_boxplot(
    aes(group = interaction(error, group2)),
    position = pd,
    alpha = 0.2,
    outlier.shape = NA,   # don’t double-plot outliers since you show points
    width = 0.6) + 
  geom_point(position = pd, shape = 21, color = 'black', size = 3) + 
  scale_y_continuous(limits = c(0.80, 1), labels = scales::percent, expand = c(0, 0)) +
  labs(x = 'Simulated error', y = 'Percent recovery') +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(r = 15)),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))

ggsave('p1.pdf', p1, dpi = 300)



p2 <- ggplot(d2[grepl('unique', d2$group),], aes(error, value, fill = group)) +
  theme_bw() +
  scale_fill_manual(name = 'Simulations', values = brewer.pal(n=5, "Set1")) + 
  geom_boxplot(
    aes(group = interaction(error, group2)),
    position = pd,
    alpha = 0.2,
    outlier.shape = NA,   # don’t double-plot outliers since you show points
    width = 0.6) + 
  geom_point(position = pd, shape = 21, color = 'black', size = 3) + 
  scale_y_continuous(limits = c(0.80, 1), labels = scales::percent, expand = c(0, 0)) +
  labs(x = 'Simulated error', y = 'Percent recovery') +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(r = 15)),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))

ggsave('p2.pdf', p2, dpi = 300)








d2 <- bind_rows(lapply(split(d, 1:nrow(d)), function(x){
        bind_rows(tibble(error = x$error, group = paste0('seed ', x$seed, ' unique'), group2 = paste0(x$error, '_unique'), value = x$percentUniqueRecovery),
                  tibble(error = x$error, group = paste0('seed ', x$seed, ' total'),  group2 = paste0(x$error, '_total'), value = x$percentTotalRecovery))
      }))

d2$group <- factor(d2$group, levels = c("seed 1 unique", "seed 1 total", "seed 2 unique", "seed 2 total", "seed 3 unique", "seed 3 total", "seed 4 unique", "seed 4 total", "seed 5 unique", "seed 5 total" ))
d2$group2 <- factor(d2$group2, levels = c("0.0%_unique", "0.0%_total", "1.0%_unique", "1.0%_total", "2.0%_unique", "2.0%_total", "3.0%_unique", "3.0%_total"))

set.seed(1)
pd <- position_dodge(width = 0.5)

p1 <- ggplot(d2, aes(error, value, fill = group)) +
      theme_bw() +
      scale_fill_manual(name = 'Syn data seeds', values = brewer.pal(n=10, "Paired")) + 
      geom_boxplot(
          aes(group = interaction(error, group2)),
          position = pd,
          alpha = 0.2,
          outlier.shape = NA,   # don’t double-plot outliers since you show points
          width = 0.6) + 
        geom_point(position = pd, shape = 21, color = 'black', size = 3) + 
        scale_y_continuous(limits = c(0.85, 1), labels = scales::percent, expand = c(0, 0)) +
        labs(x = 'Simulated error', y = 'Percent recovery') +
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
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 28),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))


ggsave('p3.pdf', p3, dpi = 300)
