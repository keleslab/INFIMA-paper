## ----------------------------------------------------------
##                  run DESeq pairwisely
## ----------------------------------------------------------

### Pairwise DESeq analysis as well as computation of the percentage of concordant pairs
# run DESeq2 between two strains
# all 8 strains

library(pheatmap, quietly = T)
library(ggplot2, quietly = T)
library(DESeq2, quietly = T)
library(DEFormats, quietly = T)

setwd('data-analysis/RNA-seq/RData/')
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


for (i in 1:8) {
  for (j in 1:8) {
    if (j > i) {
      s1 <- strains[i]
      s2 <- strains[j]
      print(paste0(s1, '-', s2))
      outdir <- paste0('data-analysis/RNA-seq/DESeq2/', s1, '-', s2)
      dir.create(outdir)
      
      load('countmat_13568_protein_coding_UQUA.RData')
      
      strain <-
        unlist(lapply(as.character(colnames(countmat)), function(x) {
          strsplit(x, '-')[[1]][1]
        }))
      gender <-
        unlist(lapply(as.character(colnames(countmat)), function(x) {
          strsplit(x, '-')[[1]][2]
        }))
      
      keep <- which(strain %in% c(s1, s2))
      countmat.s1.s2 <- countmat[, keep]
      strain <- strain[keep]
      gender <- gender[keep]
      
      storage.mode(countmat.s1.s2) <-
        'integer' # countmat has to be integer
      dds <-
        DESeqDataSetFromMatrix(countmat.s1.s2,
                               DataFrame(strain, gender),
                               ~ gender + strain)
      
      dds$strain <- factor(dds$strain, levels = unique(strain))
      dds$gender <- factor(dds$gender, levels = unique(gender))
      
      dds$strain <- relevel(dds$strain, ref = s1)
      dds$gender <- relevel(dds$gender, ref = "M")
      
      dds <- DESeq(dds, betaPrior = FALSE)
      res <- results(dds)
      res$padj[is.na(res$padj)] <- 1
      save(res,
           file = paste0(outdir, '/DESeqResults.RData'))
      
    }
  }
}




#### Analyze the output from PairwiseDESeq
library(data.table)
library(pheatmap)
library(corrplot)
library(ggplot2)


prop.concord <- matrix(0, nrow = 8, ncol = 8)

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



for (i in 1:8) {
  for (j in 1:8) {
    if (j > i) {
      s1 <- strains[i]
      s2 <- strains[j]
      print(paste0(s1, '-', s2))
      outdir <- paste0('data-analysis/RNA-seq/DESeq2/', s1, '-', s2)
      load(paste0(outdir, '/DESeqResults.RData'))
      
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
      
      n.s1 <- nrow(s1.eGenes)
      n.concord.s1 <- sum(s1.eGenes[, 1] >= s1.eGenes[, 2])
      
      # s2.specifc genes on s1, s2 ATAC-seq data
      s2.eGenes <- pro.atac.n[dt$col == 2]
      s2.eGenes.mat <- as.matrix(s2.eGenes)
      s2.eGenes <-
        as.data.table(s2.eGenes.mat[, which(colnames(pro.atac.n) %in% c(s1, s2))])
      
      n.s2 <- nrow(s2.eGenes)
      n.concord.s2 <- sum(s2.eGenes[, 2] >= s2.eGenes[, 1])
      
      # save(n.s1, n.concord.s1,
      #      n.s2, n.concord.s2,
      #      file = paste0(outdir, '/ConcordanceSummary.RData'))
      prop.concord[i, j] <-
        (n.concord.s1 + n.concord.s2) / (n.s1 + n.s2)
      prop.concord[j, i] <- prop.concord[i, j]
      
    }
  }
}

save(prop.concord, file = 'prop.concord.RData')

summary(as.vector(prop.concord[upper.tri(prop.concord)]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3563  0.5903  0.6634  0.6337  0.7022  0.7765
# median 2/3 of the eGenes are concordant


strains <- c('129', 'AJ', 'B6', 'Cast', 'NOD', 'NZO', 'PWK', 'WSB')
rownames(prop.concord) <- colnames(prop.concord) <- strains
prop.concord[!lower.tri(prop.concord)] <- NA
pheatmap(prop.concord,
         cluster_rows = F,
         cluster_cols = F)

prop.concord[!lower.tri(prop.concord)] <- NA
melted_cormat <- melt(prop.concord, na.rm = TRUE)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = (min(melted_cormat$value) + max(melted_cormat$value)) /
      2,
    limit = c(min(melted_cormat$value), max(melted_cormat$value)),
    space = "Lab",
    name = "Pearson\nCorrelation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    size = 12,
    hjust = 1
  ))



col1 <-
  colorRampPalette(
    c(
      "#7F0000",
      "red",
      "#FF7F00",
      "yellow",
      "white",
      "cyan",
      "#007FFF",
      "blue",
      "#00007F"
    )
  )
col2 <-
  colorRampPalette(
    c(
      "#67001F",
      "#B2182B",
      "#D6604D",
      "#F4A582",
      "#FDDBC7",
      "#FFFFFF",
      "#D1E5F0",
      "#92C5DE",
      "#4393C3",
      "#2166AC",
      "#053061"
    )
  )
col3 <- colorRampPalette(c("red", "white", "blue"))
col4 <-
  colorRampPalette(
    c(
      "#7F0000",
      "red",
      "#FF7F00",
      "yellow",
      "#7FFF7F",
      "cyan",
      "#007FFF",
      "blue",
      "#00007F"
    )
  )

pdf('corrplot.pdf', onefile = F)
corrplot(
  prop.concord,
  addrect = 2,
  is.corr = F,
  type = 'lower',
  diag = F,
  tl.col = 'black',
  tl.cex = 1.3,
  col = col4(20)
)
dev.off()