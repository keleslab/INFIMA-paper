#### further analysis of the simulated results ####



#######################################################################
#### Step 1: get ROC/PR curve and area under curve for V.g
#### detection rate of the causal gene
#### power analysis: TPR is the statistical power
rm(list = ls())
library(data.table)

make.data.table <- function(outdir, name){
  load(paste0(outdir, 'AUC.RData'))
  load(paste0(outdir, 'Evaluation.RData'))
  dt <- data.table(auroc = auroc,
                   aupr = aupr,
                   tau = tau,
                   power = unlist(sapply(eva, function(x) x$power)),
                   type = name)
  return(dt)
}


outdir <- '/p/keles/Collab_2014/volumeK/model/simulation/real-data/'
dt1 <- make.data.table(outdir, 'real-data')

outdir <- '/p/keles/Collab_2014/volumeK/model/simulation/non-informative/'
dt2 <- make.data.table(outdir, 'non-informative')

outdir <- '/p/keles/Collab_2014/volumeK/model/simulation/super-informative/'
dt3 <- make.data.table(outdir, 'super-informative')

auc <- rbind(dt1, dt2, dt3)
save(auc, file = '/p/keles/Collab_2014/volumeK/model/simulation/auc.RData')



#######################################################################
#### Step 2: for each simulation compare 5 baseline strategies
#### the minimum set required to contain the causal SNP

rm(list = ls())
library(data.table)
make.eva.dt <- function(outdir, name, n.sim = 3){
  load(paste0(outdir, 'Evaluation.RData'))
  dt <- data.table(p.g = unlist(sapply(eva[1:n.sim], function(x) x$p.g)),
                   worst = unlist(sapply(eva[1:n.sim], function(x) x$worst)),
                   best = unlist(sapply(eva[1:n.sim], function(x) x$best)),
                   ctgene = unlist(sapply(eva[1:n.sim], function(x) x$ctgene)),
                   ctmark = unlist(sapply(eva[1:n.sim], function(x) x$ctmark)),
                   ran = unlist(sapply(eva[1:n.sim], function(x) x$ran)),
                   type = name)
  return(dt)
}


outdir <- '/p/keles/Collab_2014/volumeK/model/simulation/real-data/'
dt1 <- make.eva.dt(outdir, 'real-data', 3)

outdir <- '/p/keles/Collab_2014/volumeK/model/simulation/non-informative/'
dt2 <- make.eva.dt(outdir, 'non-informative', 3)

outdir <- '/p/keles/Collab_2014/volumeK/model/simulation/super-informative/'
dt3 <- make.eva.dt(outdir, 'super-informative', 3)

eva <- rbind(dt1, dt2, dt3)
save(eva, file = '/p/keles/Collab_2014/volumeK/model/simulation/eva.RData')






#### ploting auc
library(data.table)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)
load('~/Dropbox/EM/simulation/auc.RData')
auc$id <- 1:nrow(auc)
auc2 <- melt(auc, id.vars = c('id', 'type'), 
             measure.vars = c('auroc', 'aupr','tau','power'))

#### p1: AUROC
p1 <- ggplot(auc2[variable == 'auroc'], aes(x = type, y = value, fill = type)) + 
  geom_boxplot() + theme_classic() + theme(legend.position = 'none', axis.text.x = element_blank()) + 
  labs(y = 'AUROC', x = '', fill = 'Group') + scale_fill_aaas()
my_comparisons <- list( c("real-data", "super-informative"),
                        c("non-informative", "super-informative"))
p1 <- p1 + stat_compare_means(comparisons = my_comparisons)
p1

#### p2: AUPR
p2 <- ggplot(auc2[variable == 'aupr'], aes(x = type, y = value, fill = type)) + 
  geom_boxplot() + theme_classic() + theme(legend.position = 'none', axis.text.x = element_blank()) + 
  labs(y = 'AUPR', x = '', fill = 'Group') + scale_fill_aaas()
my_comparisons <- list( c("real-data", "super-informative"),
                        c("non-informative", "super-informative"))
p2 <- p2 + stat_compare_means(comparisons = my_comparisons)
p2


#### p3: power at FDR 0.05
p3 <- ggplot(auc2[variable == 'power'], aes(x = type, y = value, fill = type)) + 
  geom_boxplot() + theme_classic() + theme(legend.position = 'none', axis.text.x = element_blank()) + 
  labs(y = 'Power', x = '', fill = 'Group') + scale_fill_aaas()
my_comparisons <- list( c("real-data", "super-informative"),
                        c("non-informative", "real-data"))
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

# pdf('~/Dropbox/EM/simulation/auc.pdf', width = 8, height = 4,onefile = F)
# grid_arrange_shared_legend(p1, p2, p3, nrow = 1)
# dev.off()

p1 + p2 + p3


#### plotting eva
load('~/Dropbox/EM/simulation/eva.RData')
eva$id <- 1:nrow(eva)
eva2 <- melt(eva, id.vars = c('id','p.g','type'),
             measure.vars = c('worst','best','ctgene','ctmark','ran'))
eva2$group <- paste0(eva2$type,'-',eva2$variable)
eva2[variable == 'worst']$variable <- 'Least Likely'
eva2[variable == 'best']$variable <- 'Most Likely'
eva2[variable == 'ctgene']$variable <- 'Closest to Gene'
eva2[variable == 'ctmark']$variable <- 'Closest to Marker'
eva2[variable == 'ran']$variable <- 'Random'


p4 <- ggplot(eva2, aes(x = variable, y = value/p.g, fill = type)) + geom_boxplot() +
  theme_classic() + theme(legend.position = 'top') + 
  labs(y = 'Proportion', x = 'Strategy', fill = 'Group') + scale_fill_aaas()
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
                 y_position = c(1.05, 1.20), xmin=c(1.75, 2.0), xmax=c(2.0, 2.25))
p4 <- p4 + scale_y_continuous(breaks = seq(0, 1, 0.25))

pdf('~/Dropbox/EM/simulation/Figure.pdf', width = 7, height = 7)
(p1 | p2 | p3)/p4
dev.off()
pdf('~/Documents/INFIMA-manuscript/Figure5/Figure5.pdf', width = 7, height = 7)
(p1 | p2 | p3)/p4
dev.off()




# pdf('~/Dropbox/EM/simulation/eva-1.pdf', width = 8, height = 6)
# ggplot(eva2, aes(x = variable, y = value, color = variable)) + geom_boxplot() +
#   facet_grid(~type) + labs(y = 'min set cover causal SNP')
# dev.off()

# pdf('~/Dropbox/EM/simulation/eva-2.pdf', width = 8, height = 4)
# ggplot(eva2, aes(x = variable, y = value/p.g, fill = type)) + geom_boxplot() +
#   theme_classic() + theme(legend.position = 'none') + 
#   labs(y = 'Proportion of Candidates to Cover Causal SNP', x = 'Strategy', fill = 'Group') + scale_fill_aaas()
# dev.off()
