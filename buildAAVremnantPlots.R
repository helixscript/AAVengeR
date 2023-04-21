library(dplyr)
library(ggplot2)

# Read in configuration file.
configFile <- commandArgs(trailingOnly=TRUE)
if(! file.exists(configFile)) stop('Error - configuration file does not exists.')
opt <- yaml::read_yaml(configFile)

dir.create(file.path(opt$outputDir, opt$buildAAVremnantPlots_outputDir))

sites <- readRDS(file.path(opt$outputDir, opt$buildAAVremnantPlots_inputFile))

x <- lapply(split(sites, sites$sample), function(x){
        message('sample: ', x$sample[1])

        range <- seq(0, opt$buildAAVremnantPlots_ITRlength, opt$buildAAVremnantPlots_NTbinSize)
  
        d <- tibble(remnantLen = as.integer(sub('\\.\\.', '', stringr::str_extract(x$repLeaderSeqMap, '\\.\\.\\d+'))) + opt$buildAAVremnantPlots_ITRseqStart,
                    bin = cut(remnantLen, breaks = c(-Inf, range, Inf), labels = FALSE) - (opt$buildAAVremnantPlots_NTbinSize/2),
                    r = stringr::str_count(x$repLeaderSeqMap, ';'),
                    r2 = ifelse(r >= 5, '≥ 5', r)) %>% group_by(bin, r2) %>% summarise(n = n()) %>% ungroup() %>% mutate(r2 = factor(r2, levels = rev(c('0', '1', '2', '3', '4', '≥ 5'))))
  
       range2 <- (range * opt$buildAAVremnantPlots_NTbinSize) - opt$buildAAVremnantPlots_NTbinSize

       p <- ggplot(d, aes(bin, n, fill = r2)) + 
             theme_bw() +
             geom_col() +
             scale_fill_manual(name = 'Recombinations', 
                               values = rev(c('green4', 'green2', 'gold2', 'orange', 'orangered1', 'red4')), 
                               drop = FALSE) +
             scale_x_continuous(breaks = range,
                               labels = range2, 
                               limits = c(opt$buildAAVremnantPlots_NTbinSize, 
                                          cut(opt$buildAAVremnantPlots_ITRlength, breaks = c(-Inf, range, Inf), labels = FALSE) - (opt$buildAAVremnantPlots_NTbinSize/2))) +
             scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) + 
             geom_vline(xintercept = cut(opt$buildAAVremnantPlots_ITRdumbellTip1, breaks = c(-Inf, range, Inf), labels = FALSE), linetype = 'dashed') +
             geom_vline(xintercept = cut(opt$buildAAVremnantPlots_ITRdumbellTip2, breaks = c(-Inf, range, Inf), labels = FALSE), linetype = 'dashed') +
             ggtitle(paste0(x$subject[1], ' | ', x$sample[1], ' | ', formatC(n_distinct(x$posid), format="d", big.mark=","), ' sites')) + 
             labs(x = 'ITR position', y = 'Integrations') +
             guides(fill=guide_legend(nrow = 1, byrow = TRUE, reverse = TRUE)) +
             theme(text = element_text(size=16), plot.title = element_text(size = 14),
                  legend.title = element_text(size=12),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
                  legend.position="bottom", plot.margin=grid::unit(c(0.25, 0.25, 0.25, 0.25), "in")) +
         geom_point(data = tibble(x = cut(opt$buildAAVremnantPlots_ITRseqStart, breaks = c(-Inf, range, Inf), labels = FALSE) - (opt$buildAAVremnantPlots_NTbinSize/2)), aes(x, 0), 
                    size = 7, shape="\u27A1", inherit.aes = FALSE) +
         coord_cartesian(clip = "off")
         
   ### ggsave(file.path(opt$outputDir, opt$buildAAVremnantPlots_outputDir, paste0(x$trial[1], '~', x$subject[1], '~', x$sample[1], '.pdf')), p, width = opt$buildAAVremnantPlots_plotOutputWidthInches, units = 'in')
   ggsave(file.path(opt$outputDir, opt$buildAAVremnantPlots_outputDir, paste0(x$trial[1], '~', x$subject[1], '~', x$sample[1], '.png')), p, dpi = 300, width = opt$buildAAVremnantPlots_plotOutputWidthInches, units = 'in')
   saveRDS(p, file = file.path(opt$outputDir, opt$buildAAVremnantPlots_outputDir, paste0(x$trial[1], '~', x$subject[1], '~', x$sample[1], '.rds')))
   p    
})

saveRDS(x, file = file.path(opt$outputDir, opt$buildAAVremnantPlots_outputDir, 'plots.rds'))

q(save = 'no', status = 0, runLast = FALSE) 

