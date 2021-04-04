library(data.table)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)
library(dplyr)



load('~/Documents/AlanAttie4/EM/simulation/eva.RData')
eva$id <- 1:nrow(eva)
eva2 <- melt(eva, id.vars = c('id','p.g','type'),
             measure.vars = c('worst','best','ctgene','ctmark','ran'))
eva2$group <- paste0(eva2$type,'-',eva2$variable)
eva2[variable == 'worst']$variable <- 'Least Likely'
eva2[variable == 'best']$variable <- 'Most Likely'
eva2[variable == 'ctgene']$variable <- 'Closest to Gene'
eva2[variable == 'ctmark']$variable <- 'Closest to Marker'
eva2[variable == 'ran']$variable <- 'Random'


eva2$type <- recode(eva2$type, 'real-data' = 'MI',
       'non-informative' = 'NI',
       'super-informative' = 'HI')
eva2$type <- factor(eva2$type, levels = c('NI','MI','HI'))
p4 <- ggplot(eva2, aes(x = variable, y = value/p.g, fill = type)) + geom_boxplot() +
  theme_classic() + theme(legend.position = 'top') +
  labs(y = 'Proportion of candidates required to cover \n the causal local-ATAC-MV',
       x = 'Strategy', fill = 'Simulation setting') + scale_fill_aaas()
#table(eva2$group)

my_comparisons <- list(
  c('non-informative-best', 'real-data-best'),
  c('real-data-best', 'super-informative-best')
)
test_pvals <- sapply(my_comparisons, function(x){
  sample_1 <- eva2[group == x[1]]$value/eva2[group == x[1]]$p.g
  sample_2 <- eva2[group == x[2]]$value/eva2[group == x[2]]$p.g
  wilcox.test(sample_1, sample_2)$p.value
})
# > test_pvals
# [1]  1.219003e-85 2.818550e-121
p4 <- p4 + geom_signif(annotations = c('p < 2.22e-16','p < 2.22e-16'),
                       y_position = c(1.05, 1.19), xmin=c(1.75, 2.0), xmax=c(2.0, 2.25))
p4 <- p4 + scale_y_continuous(breaks = seq(0, 1, 0.25))

pdf('~/Documents/INFIMA-manuscript/Figure6/Figure6A.pdf', width = 8, height = 4,onefile = F)
p4
dev.off()
