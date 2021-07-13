## ----------------------------------------------------------
##           compute nulls of FPD changes from atSNP
## ----------------------------------------------------------


library(data.table)
library(parallel)

n.cores <- detectCores()

setwd('data-analysis/ATAC-seq/mouse/RData')
load('atSNP_output.RData')
motif.scores <- motif_score$motif.scores
is.sig <- motif.scores$pval_rank <= 0.05
atsnp.insig <- motif.scores[!is.sig]
save(atsnp.insig, file = 'atsnp.insig.RData')

snp.info <-
  fread('data-analysis/ATAC-seq/mouse/footprint/atSNP_input.txt')
load('motif_lib_total.RData')
motif.names <- names(motif_lib_total)
motif.index.insig <- rep(0, nrow(atsnp.insig))
for (i in 1:length(motif.names)) {
  print(i)
  motif.index.insig[atsnp.insig$motif == motif.names[i]] <- i
}

## tmp <- mclapply(1:length(motif.names),
##                 function(x) {print(x); motif.index.insig[atsnp.insig$motif == motif.names[x]] <<- x}, mc.cores = n.cores)

snpids <- snp.info$snpid
indexs <- mclapply(1:length(snpids),
                   function(x) {
                     print(x)
                     return(which(atsnp.insig$snpid == snpids[x]))
                   }, mc.cores = n.cores)
snp.index.insig <- rep(0, nrow(atsnp.insig))
for (i in 1:length(snpids)) {
  print(i)
  snp.index.insig[indexs[[i]]] <- i
}



#### get NULL distribution ####

load('atSNP_output.RData')
load('fpd.129.RData')
load('fpd.AJ.RData')
load('fpd.B6.RData')
load('fpd.Cast.RData')
load('fpd.NOD.RData')
load('fpd.NZO.RData')
load('fpd.PWK.RData')
load('fpd.WSB.RData')
motif.scores <- motif_score$motif.scores
is.insig <- motif.scores$pval_rank > 0.05

fpd.dt <- data.table(
  fpd.129 = fpd.129,
  fpd.AJ = fpd.AJ,
  fpd.B6 = fpd.B6,
  fpd.Cast = fpd.Cast,
  fpd.NOD = fpd.NOD,
  fpd.NZO = fpd.NZO,
  fpd.PWK = fpd.PWK,
  fpd.WSB = fpd.WSB
)
fpd.dt <- fpd.dt[is.insig]
# > dim(fpd.dt)
# [1] 57430844        8
# save(fpd.dt, file = 'fpd.dt.RData')

load('ATAC-QTL.RData')

genotype <- snpData[snp.index.insig, 8:15]
genotype <- as.matrix(genotype)
rownames(genotype) <- NULL
fpd.mat <- as.matrix(fpd.dt)
colnames(fpd.mat) <- colnames(genotype)

fpd.ref <- mclapply(1:nrow(genotype), function(i) {
  print(i)
  val <- fpd.mat[i, ]
  group <- genotype[i, ]
  return(mean(val[group == 0]))
}, mc.cores = 20)
fpd.ref <- unlist(fpd.ref)
# save(fpd.ref, file = 'fpd.ref.RData')

# > summary(fpd.ref)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# -15.07143  -0.09995   0.26460   0.16630   0.55411   1.00000
# > length(fpd.ref)
# [1] 4502748

fpd.snp <- mclapply(1:nrow(genotype), function(i) {
  print(i)
  val <- fpd.mat[i, ]
  group <- genotype[i, ]
  return(mean(val[group == 2]))
}, mc.cores = 20)
fpd.snp <- unlist(fpd.snp)

# save(fpd.snp, file = 'fpd.snp.RData')

atsnp.insig$fpd.ref <- fpd.ref
atsnp.insig$fpd.snp <- fpd.snp
atsnp.insig$fpd.change <- fpd.snp - fpd.ref

### compute the null distributions of fpd.change for each motif
### then get the p-value for each significant case.

# length(unique(atsnp.sig$motif))
# 1316 all the motifs are here.

OneMotifComputeEcdf <- function(motif.no) {
  print(motif.no)
  fpd <- atsnp.insig[motif == motif.names[motif.no]]$fpd.change
  ecdf(fpd)
}

null.dist <-
  mclapply(1:length(motif.names), OneMotifComputeEcdf, mc.cores = 10)
FpdNullDist <- null.dist
save(FpdNullDist, file = 'nullDistFpdChange.RData')
