## ----------------------------------------------------------
##         analyze the the simulations results
## ----------------------------------------------------------


#######################################################################
#### Step 1: get ROC/PR curve and area under curve for V.g
#### detection rate of the causal gene
#### power analysis: TPR is the statistical power
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

rm(list = ls())
library(parallel)
# load the R objects that summarize the distances to marker/gene
load('/p/keles/Collab_2014/volumeK/model/run032719/dist_summary.RData')
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

