#!/usr/bin/Rscript

# AAVengeR/buildAAVremnantPlots.R
# John K. Everett, Ph.D.
# 
# This script create ITR remnant plots for AAV experiments.
# Plots are exported as png, pdf, and rds files. 

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

# Read in the configuration file and perform basic sanity checks.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) stop('Expected at least one command line argument')
source(file.path(yaml::read_yaml(args[1])$softwareDir, 'lib.R'))
opt <- startModule(args)

dir.create(file.path(opt$outputDir, opt$buildAAVremnantPlots_outputDir), showWarnings = FALSE)

sites <- readRDS(file.path(opt$outputDir, opt$buildAAVremnantPlots_inputFile))

x <- lapply(split(sites, sites$sample), function(x){
       message('sample: ', x$sample[1])
  
       breaks <- c(seq(0, opt$buildAAVremnantPlots_ITRlength, opt$buildAAVremnantPlots_NTbinSize), Inf)
       
       d <- tibble(remnantLen = as.integer(sub('\\.\\.', '', stringr::str_extract(x$repLeaderSeqMap, '\\.\\.\\d+'))) + opt$buildAAVremnantPlots_ITRseqStart - 1,
                   bin = cut(remnantLen, breaks = breaks, right = TRUE),
                   r = stringr::str_count(x$repLeaderSeqMap, ';'),
                   r2 = ifelse(r >= 5, '≥ 5', r)) %>% 
                   group_by(bin, r2) %>% 
                      summarise(n = n()) %>% 
                   ungroup() %>% 
                   mutate(r2 = factor(r2, levels = rev(c('0', '1', '2', '3', '4', '≥ 5'))))
        
       d$label <- as.integer(stringr::str_extract(as.character(d$bin), '\\d+')) + ceiling(opt$buildAAVremnantPlots_NTbinSize / 2)
       plotLabels  <- as.integer(stringr::str_extract(as.character(levels(d$bin)), '\\d+')) + ceiling(opt$buildAAVremnantPlots_NTbinSize / 2)
        
       d$label <- factor(d$label, levels = plotLabels)
       
       tips <-  cut(c(opt$buildAAVremnantPlots_ITRdumbellTip1, opt$buildAAVremnantPlots_ITRdumbellTip2), breaks = breaks, right = TRUE)
       tipLabels <- as.integer(stringr::str_extract(as.character(tips), '\\d+')) + ceiling(opt$buildAAVremnantPlots_NTbinSize / 2)
       
       seqStartPos <- as.integer(stringr::str_extract(as.character(cut(opt$buildAAVremnantPlots_ITRseqStart, breaks = breaks, right = TRUE)), '\\d+')) + ceiling(opt$buildAAVremnantPlots_NTbinSize / 2)
       
       plotScaleLabels <- as.character(levels(d$label))
       plotScaleLabels <- plotScaleLabels[seq(1, length(plotScaleLabels), by = 2)]
       
       p <- ggplot(d[! is.na(d$label),], aes(label, n, fill = r2)) + 
             theme_bw() +
             geom_vline(xintercept = match(tipLabels, plotLabels), linetype = 'dashed', linewidth = 1.5) +
             geom_col(color = 'black', linewidth = 0.50, alpha = 0.85, show.legend = TRUE) +
             scale_fill_manual(name = 'Recombinations', 
                               values = rev(c('green4', 'green1', 'gold2', 'orange', 'orangered1', 'red4')), 
                               drop = FALSE) +
             scale_x_discrete(drop = FALSE, breaks = plotScaleLabels) +
             scale_y_continuous(limits = c(0, NA), expand = c(0, 0), labels = round) +
             labs(x = 'ITR position', y = 'Integrations') +
             guides(fill=guide_legend(nrow = 1, byrow = TRUE, reverse = TRUE)) +
             theme(text = element_text(size=16), plot.title = element_text(size = 14),
                  legend.title = element_text(size=12),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
                  legend.position="bottom", plot.margin=grid::unit(c(0.25, 0.25, 0.25, 0.25), "in")) 
      
   ggsave(file.path(opt$outputDir, opt$buildAAVremnantPlots_outputDir, paste0(x$trial[1], '~', x$subject[1], '~', x$sample[1], '.pdf')), p, device=cairo_pdf, width = opt$buildAAVremnantPlots_plotOutputWidthInches, units = 'in')
   ggsave(file.path(opt$outputDir, opt$buildAAVremnantPlots_outputDir, paste0(x$trial[1], '~', x$subject[1], '~', x$sample[1], '.png')), p, dpi = 300, width = opt$buildAAVremnantPlots_plotOutputWidthInches, units = 'in')
   saveRDS(p, file = file.path(opt$outputDir, opt$buildAAVremnantPlots_outputDir, paste0(x$trial[1], '~', x$subject[1], '~', x$sample[1], '.rds')))
   p    
})

saveRDS(x, file = file.path(opt$outputDir, opt$buildAAVremnantPlots_outputDir, 'plots.rds'))

q(save = 'no', status = 0, runLast = FALSE) 
