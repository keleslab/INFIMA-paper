## ----------------------------------------------------------
##         simulate datasets under 3 different schemes
## ----------------------------------------------------------


#### simulation of one dataset ####
library(data.table)
library(GenomicAlignments)
rm(list = ls())

#### get the input data and simulate 
setwd('INFIMA-model/RData')
load('EM_input_main.RData')
load('EM_results.RData') 

# get the pseudocounts
# real-data (MI)
load('pseudo_counts.RData')
load('pseudo_counts_components.RData')
save(PI, dist, F.g, cor.A.R, cor.A.E, 
     file = 'INFIMA-model/MI/pseudo_counts.RData')

# non-informative (NI)
# neutralilze the effects of pseudocounts
G <- length(PI)
for(g in 1:G){
  if(p.g[g] > 0){
    for(k in 1:p.g[g]){
      # F.g[[g]][k] <- 0
      # distance prior
      # PI[[g]][k] <- dist[[g]][k] + F.g[[g]][k] + abs(cor.A.R[[g]][k]) + abs(cor.A.E[[g]][k])
      PI[[g]][k] <- 1
    }
  }
}
save(PI, dist, F.g, cor.A.R, cor.A.E, 
     file = 'INFIMA-model/NI/pseudo_counts.RData')

# super-informative (HI)
# randomly choose one SNP contains footprint
# set the footprint score to be 10, which is large enough to drive the signal
G <- length(PI)
set.seed(2020)
for(g in 1:G){
  if(p.g[g] > 0){
    flag <- sample(1:p.g[g],1)
    for(k in 1:p.g[g]){
      F.g[[g]][k] <- 0
      if(k == flag){
        F.g[[g]][k] <- 10
      }
      # distance prior
      PI[[g]][k] <- dist[[g]][k] + F.g[[g]][k] + abs(cor.A.R[[g]][k]) + abs(cor.A.E[[g]][k])
    }
  }
}
save(PI, dist, F.g, cor.A.R, cor.A.E, 
     file = 'INFIMA-model/HI/pseudo_counts.RData')

