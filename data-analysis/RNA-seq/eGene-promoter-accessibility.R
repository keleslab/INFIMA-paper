## ----------------------------------------------------------
##         Compute ATAC-seq signals at eGene promoters
## ----------------------------------------------------------


### get the cut size density near the TSS for all the 13568 genes ###
library(data.table)
library(GenomicAlignments)
library(parallel)

### Part 1: get the gene locations
# gene location file:
load('data-analysis/RNA-seq/RData/countmat_13568_protein_coding_UQUA.RData') # countmat.n, dt

gene.loc <- GRanges(
  seqnames = dt$chr,
  ranges = IRanges(start = dt$start, end = dt$end),
  strand = dt$strand
)
promoter.loc <-
  promoters(gene.loc, upstream = 2000, downstream = 500)

### Part 2: get the cutsize info for all the promoters wrt all 8 strains.
file.list <-
  list.files(path = 'data-analysis/ATAC-seq/mouse/RData/', pattern = '\\_ATAC_cutsize.RData$')

cutsizes <- vector('list', length(file.list))
libsizes <- rep(0, length(file.list))
coveragesizes <- rep(0, length(file.list))
for (i in 1:length(file.list)) {
  fname <- file.list[i]
  strain <- strsplit(fname, '_')[[1]][1]
  load(fname)
  cvgsum <- cvglist[['-']] + cvglist[['+']] ## simple Rle list
  cutsizes[[i]] <- mclapply(1:length(promoter.loc),
                            function(x) {
                              print(x)
                              mean(cvgsum[[as.character(seqnames(promoter.loc[x]))]][start(promoter.loc[x]):end(promoter.loc[x])])
                            },
                            mc.cores = 20)
  libsizes[i] <- libSize
  coveragesizes[i] <- coverageSize
  cutsizes[[i]] <- unlist(cutsizes[[i]])
}

# data table for integrating ATAC-seq and RNA-seq
pro.atac <- NULL
for (i in 1:8) {
  pro.atac <- cbind(pro.atac, cutsizes[[i]])
}
pro.atac <- as.data.table(pro.atac)

for (i in 1:8) {
  fname <- file.list[i]
  strain <- strsplit(fname, '_')[[1]][1]
  colnames(pro.atac)[i] <- strain
}

save(pro.atac, libsizes, coveragesizes, file = 'data-analysis/RNA-seq/RData/promoter_ATAC_13568_genes.RData')

# normalization
# libSize <- length(bamIn) # the total number of cuts
# coverageSize <- sum(as.numeric(width(reduce(bamIn, ignore.strand=TRUE))))
# libFactor <- libSize / coverageSize # this should be larger than one, the average cut number for cut size.
libfactors.1 <- libsizes / coveragesizes
#[1] 1.679502 1.966592 2.177750 2.012178 2.076068 1.732596 1.623145 1.856826

# normalize by the sequencing depth
load('AlignmentResults.RData')
# designInfo.ref
libfactors.2 <- rep(0, 8)
for (i in 1:8) {
  libfactors.2[i] <-
    (designInfo.ref$depth[2 * i - 1] + designInfo.ref$depth[2 * i]) / 1e8
}
libfactors.2
# [1] 1.558235 1.430264 1.517541 1.920028 1.955475 1.590481 1.028524 1.428704
libfactors <- libfactors.1 * libfactors.2



pro.atac.n <- pro.atac
pro.atac.n <- as.data.frame(pro.atac.n)
for (i in 1:8) {
  pro.atac.n[, i] <- log(pro.atac.n[, i] / libfactors.2[i])
  # pro.atac.n[,i] <- log(pro.atac.n[,i])/(libfactors.1[i] + libfactors.2[i])
  # standardize each column
  # mu <- mean(pro.atac.n[,i])
  # sd <- sd(pro.atac.n[,i])
  # pro.atac.n[,i] <- (pro.atac.n[,i] - mu)/sd
}
summary(pro.atac.n)
pro.atac.n <- as.data.table(pro.atac.n)


setwd('RNA-seq/DESeq2/')
# all 8 strains
strains <- c('129', 'AJ', 'B6', 'Cast', 'NOD', 'NZO', 'PWK', 'WSB')
# color scheme for all 8 strains
colors <- c(
  rgb(240, 128, 128, maxColorValue = 255, alpha = 255),
  rgb(218, 165, 32, maxColorValue = 255, alpha = 255),
  rgb(128, 128, 128, maxColorValue = 255, alpha = 255),
  rgb(0, 160, 0, maxColorValue = 255, alpha = 255),
  rgb(16, 16, 240, maxColorValue = 255, alpha = 255),
  rgb(0, 160, 240, maxColorValue = 255, alpha = 255),
  rgb(240, 0, 0, maxColorValue = 255, alpha = 255),
  rgb(144, 0, 224, maxColorValue = 255, alpha = 255)
)

s1 <- 'B6'
for (s2 in c('129', 'AJ', 'Cast', 'NOD', 'NZO', 'PWK', 'WSB')) {
  dir <- paste0('RNA-seq/DESeq2/', s1, '-', s2, '/')
  load(paste0(dir, 'DESeqResultsNew.RData'))
  
  s1.col <- colors[which(strains == s1)]
  s2.col <- colors[which(strains == s2)]
  
  # adjust p-value using BH method
  dt <- data.table(fc = res$log2FoldChange, p = res$padj)
  dt$col <- 0
  dt$col[dt$fc < -2 & dt$p < 0.05] <- 1 # s1 specific genes
  dt$col[dt$fc > 2 & dt$p < 0.05] <- 2 # s2 specific genes
  dt$col <- as.factor(dt$col)
  
  # s1.specific genes on s1, s2 ATAC-seq data
  s1.eGenes <- pro.atac.n[dt$col == 1]
  s1.eGenes.mat <- as.matrix(s1.eGenes)
  s1.eGenes <-
    as.data.table(s1.eGenes.mat[, which(colnames(pro.atac.n) %in% c(s1, s2))])
  s1.dt <- melt(s1.eGenes, measure.vars = c(s1, s2))
  s1.dt$eGenes <- s1
  # s2.specifc genes on s1, s2 ATAC-seq data
  s2.eGenes <- pro.atac.n[dt$col == 2]
  s2.eGenes.mat <- as.matrix(s2.eGenes)
  s2.eGenes <-
    as.data.table(s2.eGenes.mat[, which(colnames(pro.atac.n) %in% c(s1, s2))])
  s2.dt <- melt(s2.eGenes, measure.vars = c(s1, s2))
  s2.dt$eGenes <- s2
  
  dt.s1s2 <- rbind(s1.dt, s2.dt)
  
  dt.s1s2$eGenes <- factor(dt.s1s2$eGenes, levels = c(s1, s2))
  dt.s1s2$variable <- factor(dt.s1s2$variable, levels = c(s1, s2))
  
  
  # # wilcox.test for the left panel: eGenes for s1
  # x <- dt.s1s2[eGenes == s1 & variable == s1]$value
  # y <- dt.s1s2[eGenes == s1 & variable == s2]$value
  # p.val.1 <- -log10(wilcox.test(x, y, alternative = "two.sided")$p.value)
  #
  # # wilcox.test for the left panel: eGenes for s2
  # x <- dt.s1s2[eGenes == s2 & variable == s1]$value
  # y <- dt.s1s2[eGenes == s2 & variable == s2]$value
  # p.val.2 <- -log10(wilcox.test(x, y, alternative = "two.sided")$p.value)
  
  pdf(paste0(dir, 'promoter_ATAC.pdf'),
      width = 5,
      height = 5)
  print(
    ggplot(dt.s1s2, aes(
      x = eGenes, y = value, fill = variable
    )) + geom_boxplot() + labs(title = paste0(s1, ' v.s. ', s2)) + scale_fill_manual(values = c(s1.col, s2.col))
  )
  dev.off()
}
