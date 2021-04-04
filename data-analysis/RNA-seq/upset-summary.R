## ----------------------------------------------------------
##         Compute ATAC-seq signals at eGene promoters
## ----------------------------------------------------------

# compute the number of top local-ATAC-MVs associated with genes

library(data.table)
library(UpSetR)
library(parallel)

load('data-analysis/RNA-seq/RData/cis.RData')
load('data-analysis/ATAC-seq/mouse/RData/ATAC-QTL.RData')
cis.genes <- unique(cis$gene) # 6,418 genes in total
cis.snps <-
  unlist(mclapply(cis.genes, function(x)
    cis[gene == x][order(FDR)][1]$SNP, mc.cores = 20))

genotypes <- NULL
for (i in 1:length(cis.snps)) {
  print(i)
  genotypes <-
    rbind(genotypes, snpData[snpData$Name == cis.snps[i], 8:15])
}

genotypes <- as.data.table(genotypes)
genotypes <- as.matrix(genotypes)
mode(genotypes) <- "integer"
genotypes[genotypes == 2] <- 1
genotypes <- as.data.table(genotypes)

# all 8 strains
strains <- c('129', 'AJ', 'B6', 'CAST', 'NOD', 'NZO', 'PWK', 'WSB')
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

# only include the sets with more than 50 genes
metadata <- data.table(strains = strains)
metadata$names <- as.character(metadata$strains)

colnames(genotypes) <- gsub('Cast', 'CAST', colnames(genotypes))

pdf(
  'upset-plot.pdf',
  onefile = F,
  width = 8,
  height = 5.5
)
upset(
  genotypes,
  sets = c('CAST', 'PWK', 'WSB', 'NOD', 'AJ', 'NZO', '129', 'B6'),
  #sets = c('129', 'AJ', 'B6', 'CAST', 'NOD', 'NZO', 'PWK', 'WSB'),
  main.bar.color = 'black',
  keep.order = T,
  order.by = 'freq',
  nintersects = 23,
  text.scale = 1.5,
  matrix.dot.alpha = 0.1,
  sets.bar.color = colors[c(4, 7, 8, 5, 2, 6, 1, 3)],
  mainbar.y.label = "Number of top local-ATAC-QTL associations\n with the allele patterns",
  sets.x.label = "Set size",
  set.metadata = list(data = metadata,
                      plots = list(
                        list(
                          type = "matrix_rows",
                          column = "names",
                          colors = c(
                            `129` = colors[1],
                            AJ = colors[2],
                            B6 = colors[3],
                            CAST = colors[4],
                            NOD = colors[5],
                            NZO = colors[6],
                            PWK = colors[7],
                            WSB = colors[8]
                          ),
                          alpha = 1
                        )
                      ))
)
dev.off()
