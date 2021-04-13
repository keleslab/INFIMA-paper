## ----------------------------------------------------------
##         analyze the the simulations results
## ----------------------------------------------------------


#######################################################################
#### Step 1: get ROC/PR curve and area under curve for V.g
#### detection rate of the causal gene
#### power analysis: TPR is the statistical power
#######################################################################

library(ROCit)
library(caTools) # auc <- trapz(recall, precision)
library(parallel)

rm(list = ls())

# consider V.g >= tau as causal
direct.pprob <- function(p, fdr = 0.05){
  err <- 1 - p # the error rate
  err <- sort(err) # ascending order
  tau <- 1 - err[max(which( (cumsum(err)/1:length(err) <= fdr) == 1 ))]
  return(tau)
}

# compute AUC and tau for each simulation
compute.auc <- function(i, outdir){
  print(i)
  filename <- paste0(outdir,'/sim/sim-',i,'.RData')
  load(filename)
  tmp <- measureit(score = params.sim$V.g,
                   class = params.real$V.g,
                   measure = c('TPR','FPR','PPV'))
  tmp$PPV[1] <- 1
  auroc <- trapz(tmp$FPR, tmp$TPR)
  aupr <- trapz(tmp$TPR, tmp$PPV)
  tau <- direct.pprob(params.sim$V.g)
  # print(paste0('AUROC = ', format(auroc[i], digits = 6),
  #              ' , AUPR = ', format(aupr[i], digits = 6),
  #              ' , tau = ', format(tau[i], digits = 6)))
  return(list(auroc = auroc, aupr = aupr, tau = tau))
}

# save the AUC and direct pprob data
save.auc <- function(n.sim, outdir){
  res <- mclapply(1:n.sim, function(x){
    compute.auc(x, outdir)
  }, mc.cores = 20)
  auroc <- sapply(res, function(x) x$auroc)
  aupr <- sapply(res, function(x) x$aupr)
  tau <- sapply(res, function(x) x$tau)
  save(auroc, aupr, tau, file = paste0(outdir, 'AUC.RData'))
}

n.sim <- 1000
outdir <- 'INFIMA-model/simulation/NI'
save.auc(n.sim, outdir)
outdir <- 'INFIMA-model/simulation/MI'
save.auc(n.sim, outdir)
outdir <- 'INFIMA-model/simulation/HI'
save.auc(n.sim, outdir)







#######################################################################
#### Step 2: for each simulation compare 5 baseline strategies
#### the minimum set required to contain the causal SNP
#######################################################################

rm(list = ls())
library(parallel)
# load the R objects that summarize the distances to marker/gene
load('dist_summary.RData')
# dist.gene, dist.mark


compute.min.set <- function(i, .tau, outdir){
  print(i)
  filename <- paste0(outdir,'/sim/sim-',i,'.RData')
  load(filename)
  tot.actual.pos <- sum(params.real$V.g)
  indicator <- params.sim$V.g >= .tau
  params.sim$V.g <- params.sim$V.g[indicator]
  params.sim$Z.g <- params.sim$Z.g[indicator]
  params.real$V.g <- params.real$V.g[indicator]
  params.real$Z.g <- params.real$Z.g[indicator]
  d.gene <- dist.gene[indicator]
  d.mark <- dist.mark[indicator]
  
  # ~95% of the V.g after direct pprob approach are true positives, TPR
  # a.k.a. statistical power
  power <- sum(params.real$V.g)/tot.actual.pos
  
  # retain the cases where real V.g = 1
  indicator <- params.real$V.g == 1
  params.sim$V.g <- params.sim$V.g[indicator]
  params.sim$Z.g <- params.sim$Z.g[indicator]
  params.real$V.g <- params.real$V.g[indicator]
  params.real$Z.g <- params.real$Z.g[indicator]
  d.gene <- d.gene[indicator]
  d.mark <- d.mark[indicator]
  
  
  # compute the min set to contain causal SNP for each strategy
  G <- length(d.gene)
  p.g <- rep(0, G) # the total set
  best <- rep(0, G)
  worst <- rep(0, G)
  ctgene <- rep(0, G)
  ctmark <- rep(0, G)
  ran <- rep(0, G)
  
  for(g in 1:G){
    p.g[g] <- length(params.real$Z.g[[g]])
    ind <- which(params.real$Z.g[[g]] == 1)
    # 1: worst strategy
    worst[g] <- rank(params.sim$Z.g[[g]], ties.method = 'average')[ind]
    # 2: best strategy
    best[g] <- p.g[g] - worst[g]
    # 3: closest to gene
    ctgene[g] <- rank(d.gene[[g]], ties.method = 'average')[ind]
    # 4: closest to marker
    ctmark[g] <- rank(d.mark[[g]], ties.method = 'average')[ind]
    # 5: random (use expectation, p.g / 2)
    ran[g] <- sample(1:p.g[g],1)
  }

  strategy <- list(power = power,
                   p.g = p.g, 
                   worst = worst, best = best,
                   ctgene = ctgene, ctmark = ctmark, ran = ran)
  return(strategy)
}



# save the power and evaluation data
save.eva <- function(n.sim, outdir){
  load(paste0(outdir, 'AUC.RData')) # use tau to filter causal genes
  eva <- mclapply(1:n.sim, function(x){
    compute.min.set(x, tau[x], outdir)
  }, mc.cores = 20)
  # power <- unlist(sapply(res, function(x) x$power))
  # p.g <- unlist(sapply(res, function(x) x$p.g))
  # worst <- unlist(sapply(res, function(x) x$worst))
  # best <- unlist(sapply(res, function(x) x$best))
  # ctgene <- unlist(sapply(res, function(x) x$ctgene))
  # ctmark <- unlist(sapply(res, function(x) x$ctmark))
  # ran <- unlist(sapply(res, function(x) x$ran))
  
  save(eva, file = paste0(outdir, 'Evaluation.RData'))
}

n.sim <- 1000
set.seed(2020)
outdir <- 'INFIMA-model/simulation/NI'
save.eva(n.sim, outdir)

outdir <- 'INFIMA-model/simulation/MI'
save.eva(n.sim, outdir)

outdir <- 'INFIMA-model/simulation/HI'
save.eva(n.sim, outdir)


#######################################################################
#### Step 3: Make ROC/PR curve plots for one simulation
#######################################################################


get.roc.pr <- function(outdir, name, i){
  print(i)
  filename <- paste0(outdir,'/sim/sim-',i,'.RData')
  load(filename)
  tmp <- measureit(score = params.sim$V.g,
                   class = params.real$V.g,
                   measure = c('TPR','FPR','PPV'))
  tmp$PPV[1] <- 1
  auroc <- trapz(tmp$FPR, tmp$TPR)
  aupr <- trapz(tmp$TPR, tmp$PPV)
  dt <- data.table(TPR = tmp$TPR, FPR = tmp$FPR, PPV = tmp$PPV,
                   type = name)
  return(list(dt = dt, auroc = auroc, aupr = aupr))
}

outdir <- 'INFIMA-model/simulation/MI'
res1 <- get.roc.pr(outdir, 'MI', 1)

outdir <- 'INFIMA-model/simulation/NI'
res2 <- get.roc.pr(outdir, 'NI', 1)

outdir <- 'INFIMA-model/simulation/HI'
res3 <- get.roc.pr(outdir, 'HI', 1)

roc.pr <- rbind(res1$dt, res2$dt, res3$dt)

save(roc.pr, file = 'INFIMA-model/simulation/roc.pr.RData')



#######################################################################
##### Step 4: Make power curve plots for one simulation #####
#######################################################################

# compute FDR and power for each simulation
fdr.power <- function(outdir, name, i){
  print(i)
  filename <- paste0(outdir,'/sim/sim-',i,'.RData')
  load(filename)
  tmp <- data.table(score = params.sim$V.g,
                    class = params.real$V.g,
                    err = 1 - params.sim$V.g)
  tmp <- tmp[order(err)]
  tot.actual.pos <- sum(params.real$V.g)
  tmp[, fdr := cumsum(err)/1:length(err)]
  tmp[, power := cumsum(class)/tot.actual.pos]
  tmp$type <- name
  return(tmp)
}

outdir <- 'INFIMA-model/simulation/MI'
dt1 <- fdr.power(outdir, 'MI', 1)

outdir <- 'INFIMA-model/simulation/NI'
dt2 <- fdr.power(outdir, 'NI', 1)

outdir <- 'INFIMA-model/simulation/HI'
dt3 <- fdr.power(outdir, 'HI', 1)

fdr.power <- rbind(dt1, dt2, dt3)

save(fdr.power, file = 'INFIMA-model/simulation/fdr.power.RData')




library(data.table)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)
library(dplyr)

# library(data.table)
# library(ggplot2)
# pdf('fdr.power.pdf', width = 8, height = 4)
# ggplot(fdr.power, aes(x = fdr, y = power, color = type)) + geom_line()
# dev.off()

# pdf('ROC.pdf', width = 8, height = 4)
# ggplot(roc.pr, aes(x = FPR, y = TPR, color = type)) + geom_line()
# dev.off()
# pdf('PR.pdf', width = 8, height = 4)
# ggplot(roc.pr, aes(x = TPR, y = PPV, color = type)) + geom_line()
# dev.off()


load('eva.RData')
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

pdf('eva.pdf', width = 8, height = 4,onefile = F)
p4
dev.off()



load('auc.RData')
auc$id <- 1:nrow(auc)
auc2 <- melt(auc, id.vars = c('id', 'type'), 
             measure.vars = c('auroc', 'aupr','tau','power'))

auc2$type <- recode(auc2$type, 'real-data' = 'MI',
                    'non-informative' = 'NI',
                    'super-informative' = 'HI')
auc2$type <- factor(auc2$type, levels = c('NI','MI','HI'))


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

pdf('eva-supp.pdf', width = 8, height = 4,onefile = F)
grid_arrange_shared_legend(p1, p2, p3, nrow = 1)
dev.off()