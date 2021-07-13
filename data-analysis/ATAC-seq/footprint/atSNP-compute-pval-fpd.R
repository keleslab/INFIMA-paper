## ----------------------------------------------------------
##           compute p-values of FPD changes from atSNP
## ----------------------------------------------------------


## combine the footprint depth FPD ###
## compute pval_fpd for each row of atsnp.sig ###


library(data.table)
library(parallel)

setwd('data-analysis/ATAC-seq/mouse/RData')
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
is.sig <- motif.scores$pval_rank <= 0.05
# > sum(is.sig)
# [1] 4502748
atsnp.sig <- motif.scores[is.sig]
save(atsnp.sig, file = 'atsnp.sig.RData')

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
fpd.dt <- fpd.dt[is.sig]
# > dim(fpd.dt)
# [1] 4502748       8
# save(fpd.dt, file = './fpd.dt.RData')

load('ATAC-QTL.RData')
load('snp.index.RData')
genotype <- snpData[snp.index, 8:15]
genotype <- as.matrix(genotype)
rownames(genotype) <- NULL
fpd.mat <- as.matrix(fpd.dt)
colnames(fpd.mat) <- colnames(genotype)

library(parallel)

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

atsnp.sig$fpd.ref <- fpd.ref
atsnp.sig$fpd.snp <- fpd.snp
atsnp.sig$fpd.change <- fpd.snp - fpd.ref

### compute the null distributions of fpd.change for each motif
### then get the p-value for each significant case.

# length(unique(atsnp.sig$motif))
# 1316 all the motifs are here.

load('atsnp.sig.RData')
load('fpd.ref.RData')
load('fpd.snp.RData')

atsnp.sig$fpd.change <- fpd.snp - fpd.ref

load('motif_lib_total.RData')
motif.names <- names(motif_lib_total)
load('nullDistFpdChange.RData')

### for each row of atsnp.sig, compute the p-value for each fpd.change
OneFpdChangeComputePval <- function(row.no) {
  print(row.no)
  Fn <- FpdNullDist[[which(atsnp.sig$motif[row.no] == motif.names)]]
  fpd.change <- atsnp.sig$fpd.change[row.no]
  pval_fpd <- Fn(-abs(fpd.change)) + 1 - Fn(abs(fpd.change))
  pval_fpd
}

pval_fpd <-
  unlist(mclapply(1:nrow(atsnp.sig), OneFpdChangeComputePval, mc.cores = 10))
atsnp.sig$pval_fpd <- pval_fpd
# save(pval_fpd, file = 'pval_fpd.RData')
# > summary(pval_fpd)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.2541  0.5041  0.5030  0.7528  1.0000

atsnp.sig$fpd.ref <- fpd.ref
atsnp.sig$fpd.snp <- fpd.snp
save(atsnp.sig, file = 'atsnp.sig.RData')





# atsnp.sig[pval_fpd <= 0.05, .N] # 218280
# library(ggplot2)
# png('./histogram_pval_fpd.png', height = 600, width = 600)
# pval_fpd <- atsnp.sig$pval_fpd
# qplot(pval_fpd, geom = 'histogram')
# dev.off()


### check the gain.of.func subset

gain.of.func2 <-
  atsnp.sig[pval_ref > 0.05 & pval_snp <= 0.05 & pval_fpd <= 0.05]
nrow(gain.of.func2) # 55488
# gain.of.func + creation of a footprint
gain.of.func2[fpd.change < 0, .N]
# [1] 30820
gain.of.func2[fpd.change == 0, .N]
# [1] 0
gain.of.func2[fpd.change > 0, .N]
# [1] 24668 the number of consistent cases



### check the loss.of.func subset

loss.of.func2 <-
  atsnp.sig[pval_ref <= 0.05 & pval_snp > 0.05 & pval_fpd <= 0.05]
nrow(loss.of.func2) # 55774
loss.of.func2[fpd.change < 0, .N]
# [1] 35056 the number of consistent cases
loss.of.func2[fpd.change == 0, .N]
# [1] 0
loss.of.func2[fpd.change > 0, .N]
# [1] 20718



# motif enhance & footprint enhance
is.enhance <-
  atsnp.sig$pval_ref > 0.05 &
  atsnp.sig$pval_snp <= 0.05 & atsnp.sig$fpd.change > 0
sum(is.enhance) # 640038
# motif disrupt & footprint disrupt
is.disrupt <-
  atsnp.sig$pval_ref <= 0.05 &
  atsnp.sig$pval_snp > 0.05 & atsnp.sig$fpd.change < 0
sum(is.disrupt) # 571769
is.consistent <- is.enhance | is.disrupt
sum(is.consistent) # 1211807

save(is.consistent, is.enhance, is.disrupt, file = 'is.consistent.RData')




atsnp <- atsnp.sig[is.consistent] # 1,211,807 cases

dt <- atsnp
snps <- unique(dt$snpid) # 47,057 out of 47,062 are involved

# use FDR adjustments for each SNP separately
# adjust for pval_fpd not for pval_rank

one_snp_fdr <- function(i) {
  print(i)
  current_snp <- snps[i]
  tmp <- dt[snpid == current_snp]
  pval_fpd_adj <- p.adjust(tmp$pval_fpd, method = 'BH')
  min_pval_rank <- min(tmp$pval_rank)
  min_pval_fpd_adj <- min(pval_fpd_adj)
  pval_rank_motifs <- tmp[pval_rank < min_pval_rank + 1e-2]$motif
  pval_fpd_motifs <-
    tmp[pval_fpd_adj < min_pval_fpd_adj + 1e-2]$motif
  
  return(
    list(
      pval_rank = min_pval_rank,
      pval_rank_motifs = pval_rank_motifs,
      pval_fpd = min_pval_fpd_adj,
      pval_fpd_motifs = pval_fpd_motifs
    )
  )
}

dt.list <- mclapply(1:length(snps), one_snp_fdr, mc.cores = 20)


## FDR adjustment across all SNPs
# we have already filtered pval_rank < 0.05
# it is meaningless to adjust for pval_rank again

pval_fpds <- unlist(lapply(dt.list, function(x)
  x$pval_fpd))



indicator <- pval_fpds < 0.05 # only 1,451 cases survives

dt.list <- dt.list[indicator]
snps <- snps[indicator]

# for each SNP, we check the overlapping between motif sets predicted by pval_rank
# as well as the motif sets predicted by pval_fpd

is.common <- rep(0, length(dt.list))
for (i in 1:length(dt.list)) {
  pval_rank_motifs <- as.character(dt.list[[i]]$pval_rank_motifs)
  pval_fpd_motifs <- as.character(dt.list[[i]]$pval_fpd_motifs)
  common <- intersect(pval_rank_motifs, pval_fpd_motifs)
  dt.list[[i]]$pval_rank_motifs <- common
  dt.list[[i]]$pval_fpd_motifs <- common
  is.common[i] <- length(common)
}

# > sum(is.common == 0)
# [1] 101
dt.list <- dt.list[is.common > 0] # 1,350 SNPs left
snps <- snps[is.common > 0] # 1,350 SNPs left

# then create the data.table for all the final results

res <- NULL
for (i in 1:length(dt.list)) {
  res <-
    rbind(res, dt[snpid == snps[i]][motif %in% dt.list[[i]]$pval_fpd_motifs])
}
res <- as.data.table(res) # 8,029 cases in total
atsnp.final <- res

# > length(unique(atsnp.final$snpid))
# [1] 1350
# > length(unique(atsnp.final$motif))
# [1] 1196

save(atsnp.final, file = 'atsnp.final.RData')



## pval_rank vs. pval_fpd
library(ggplot2)
rm(list = ls())
load('atsnp.final.RData')
load('atsnp.sig.RData')

atsnp.sig <- atsnp.sig[fpd.change != 0]
dt <- atsnp.sig[, .(motif, snpid, pval_rank, pval_fpd)]
dt$delta_fpd <- ifelse(atsnp.sig$fpd.change > 0, 1,-1)
dt$delta_motif <-
  ifelse(atsnp.sig$pval_ref > atsnp.sig$pval_snp, 1,-1)
dt$pval_rank <- abs(-log10(dt$pval_rank)) * dt$delta_motif
dt$pval_fpd <- abs(-log10(dt$pval_fpd)) * dt$delta_fpd

dt$gain_of_func <-
  atsnp.sig$pval_ref > 0.05 & atsnp.sig$pval_snp <= 0.05
dt$loss_of_func <-
  atsnp.sig$pval_ref <= 0.05 & atsnp.sig$pval_snp > 0.05

dt <- dt[gain_of_func | loss_of_func]

table(paste0(dt$delta_motif, dt$delta_fpd))
# correct!
# 11    -11   -1-1    1-1
# 640038 577262 571769 525327

query.rows <- paste(atsnp.final$motif, atsnp.final$snpid)
subject.rows <- paste(dt$motif, dt$snpid)
ind <- match(subject.rows, query.rows)

dt$category <- 'other'
dt$category[dt$delta_motif > 0 &
              dt$delta_fpd > 0 & !is.na(ind)] <- 'enhance'
dt$category[dt$delta_motif < 0 &
              dt$delta_fpd < 0 & !is.na(ind)] <- 'disrupt'
table(dt$category)
# disrupt enhance   other
# 5467    2562 2306367



#### test ####
intercept <- -log10(0.05)
g <-
  ggplot(dt[abs(pval_fpd) > -log10(0.05)],
         aes(
           x = pval_fpd,
           y = pval_rank,
           color = category,
           alpha = category,
           size = category
         )) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c(
    'enhance' = 'blue',
    'disrupt' = 'red',
    'other' = 'grey'
  )) +
  geom_vline(xintercept = c(-intercept, intercept),
             linetype = 'dashed') +
  geom_hline(yintercept = c(-intercept, intercept),
             linetype = 'dashed')
pdf('pvalrank_vs_pvalfpd.pdf',
    width = 5,
    height = 5)
g + labs(x = expression(paste('|-log'[10], '(pval_fpd)|')),
         y = expression(paste('|-log'[10], '(pval_rank)|'))) +
  theme(legend.position = 'none') +
  annotate(
    geom = "text",
    x = 4,
    y = 7,
    label = "Enhance",
    color = "blue"
  ) +
  annotate(
    geom = "text",
    x = -4,
    y = -7,
    label = "Disrupt",
    color = "red"
  ) +
  scale_alpha_manual(values = c(
    'enhance' = 0.8,
    'disrupt' = 0.8,
    'other' = 0.4
  )) +
  scale_size_manual(values = c(
    'enhance' = 1,
    'disrupt' = 1,
    'other' = 0.8
  ))
dev.off()
