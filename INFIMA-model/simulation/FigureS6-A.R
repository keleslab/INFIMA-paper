library(data.table)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)
library(dplyr)

rm(list = ls())
load('~/Documents/AlanAttie4/EM/simulation/auc.RData')
auc$id <- 1:nrow(auc)
auc2 <- melt(auc, id.vars = c('id', 'type'), 
             measure.vars = c('auroc', 'aupr','tau','power'))

auc2$type <- recode(auc2$type, 'real-data' = 'MI',
                    'non-informative' = 'NI',
                    'super-informative' = 'HI')
auc2$type <- factor(auc2$type, levels = c('NI','MI','HI'))

setwd('~/Documents/INFIMA-manuscript/Figure6/')

#### p1: AUROC
p1 <- ggplot(auc2[variable == 'auroc'], aes(x = type, y = value, fill = type)) + 
  geom_boxplot() + theme_classic() + theme(legend.position = 'none', axis.text.x = element_blank()) + 
  labs(y = 'AUROC', x = '', fill = 'Simulation setting') + scale_fill_aaas()
my_comparisons <- list( c("MI", "HI"),
                        c("NI", "HI"))
p1 <- p1 + stat_compare_means(comparisons = my_comparisons)
p1

#### p2: AUPR
p2 <- ggplot(auc2[variable == 'aupr'], aes(x = type, y = value, fill = type)) + 
  geom_boxplot() + theme_classic() + theme(legend.position = 'none', axis.text.x = element_blank()) + 
  labs(y = 'AUPR', x = '', fill = 'Simulation setting') + scale_fill_aaas()
my_comparisons <- list( c("MI", "HI"),
                        c("NI", "HI"))
p2 <- p2 + stat_compare_means(comparisons = my_comparisons)
p2


#### p3: power at FDR 0.05
p3 <- ggplot(auc2[variable == 'power'], aes(x = type, y = value, fill = type)) + 
  geom_boxplot() + theme_classic() + theme(legend.position = 'none', axis.text.x = element_blank()) + 
  labs(y = 'Power', x = '', fill = 'Simulation setting') + scale_fill_aaas()
my_comparisons <- list( c("MI", "HI"),
                        c("NI", "MI"))
p3 <- p3 + stat_compare_means(comparisons = my_comparisons)
p3

library(gridExtra)
library(grid)
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

pdf('FigureS6-A.pdf', width = 8, height = 4,onefile = F)
grid_arrange_shared_legend(p1, p2, p3, nrow = 1)
dev.off()