
#######################################################################
##### Make ROC/PR curve plots for one simulation #####

library(ROCit)
library(caTools) # auc <- trapz(recall, precision)
library(parallel)
library(data.table)

rm(list = ls())

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

outdir <- '/p/keles/Collab_2014/volumeK/model/simulation/real-data/'
res1 <- get.roc.pr(outdir, 'real-data', 1)

outdir <- '/p/keles/Collab_2014/volumeK/model/simulation/non-informative/'
res2 <- get.roc.pr(outdir, 'non-informative', 1)

outdir <- '/p/keles/Collab_2014/volumeK/model/simulation/super-informative/'
res3 <- get.roc.pr(outdir, 'super-informative', 1)

roc.pr <- rbind(res1$dt, res2$dt, res3$dt)

save(roc.pr, file = '/p/keles/Collab_2014/volumeK/model/simulation/roc.pr.RData')

#### plotting ROC/PR curve
library(data.table)
library(ggplot2)
pdf('~/Dropbox/EM/simulation/ROC.pdf', width = 8, height = 4)
ggplot(roc.pr, aes(x = FPR, y = TPR, color = type)) + geom_line()
dev.off()
pdf('~/Dropbox/EM/simulation/PR.pdf', width = 8, height = 4)
ggplot(roc.pr, aes(x = TPR, y = PPV, color = type)) + geom_line()
dev.off()


#######################################################################
##### Make power curve plots for one simulation #####

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

outdir <- '/p/keles/Collab_2014/volumeK/model/simulation/real-data/'
dt1 <- fdr.power(outdir, 'real-data', 1)

outdir <- '/p/keles/Collab_2014/volumeK/model/simulation/non-informative/'
dt2 <- fdr.power(outdir, 'non-informative', 1)

outdir <- '/p/keles/Collab_2014/volumeK/model/simulation/super-informative/'
dt3 <- fdr.power(outdir, 'super-informative', 1)

fdr.power <- rbind(dt1, dt2, dt3)

save(fdr.power, file = '/p/keles/Collab_2014/volumeK/model/simulation/fdr.power.RData')

#### plotting FDR-power curve
library(data.table)
library(ggplot2)
pdf('~/Dropbox/EM/simulation/fdr.power.pdf', width = 8, height = 4)
ggplot(fdr.power, aes(x = fdr, y = power, color = type)) + geom_line()
dev.off()
