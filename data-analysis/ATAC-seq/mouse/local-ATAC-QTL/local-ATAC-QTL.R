## ----------------------------------------------------------
##            Identification of local-ATAC-QTL
## ----------------------------------------------------------

library(VariantAnnotation)
library(snpStats)
library(parallel)
library(lmPerm)
library(Biostrings)
library(XVector)
library(GenomicAlignments)

## VCF files for the mouse strains:
## ftp://ftp-mouse.sanger.ac.uk/current_snps/README
## wget ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
## wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz.tbi

variant_dr <- 'data-analysis/ATAC-seq/mouse/local-ATAC-QTL/'
fl <-
  file.path(variant_dr, 'mgp.v5.merged.snps_all.dbSNP142.vcf.gz')
hdr <- scanVcfHeader(fl)

## Extract SNPs in ATAC-seq peaks
# load('data-analysis/ATAC-seq/mouse/RData/IDR_master_peak_list_trimmed_with_counts.RData')
load(
  'data-analysis/ATAC-seq/mouse/RData/IDR_master_peak_list_DESeq_with_counts.RData'
)

tab <- TabixFile(fl)
vcf <- readVcf(tab, 'mm10', param = peaks)
res <- genotypeToSnpMatrix(vcf)

# retain only the founder strains
snpMat.numeric <- t(as(res$genotype, 'numeric'))
snpMat.numeric <- snpMat.numeric[, c(2, 5, 12, 16, 26, 28, 30, 35)]
# > colnames(snpMat.numeric)
# [1] '129S1_SvImJ' 'A_J'         'C57BL_6NJ'   'CAST_EiJ'    'NOD_ShiLtJ'
# [6] 'NZO_HlLtJ'   'PWK_PhJ'     'WSB_EiJ'
# Note that C57BL_6NJ is not C57BL_6J (B6)
snpMat.numeric[, 3] <- 0 # always homozygous

## Here are the SNP coordinates and more
snpRanges <- rowRanges(vcf)

## only retain SNPs with good quality
snpMat.numeric <- snpMat.numeric[snpRanges$FILTER == 'PASS',]
snpRanges <- snpRanges[snpRanges$FILTER == 'PASS']

# sum of alternative alleles
MAC <- rowSums(snpMat.numeric)
# highest quality
row.indicator <- snpRanges$QUAL == 999
# remove NA's
row.indicator <-
  as.logical(row.indicator & (!is.na(MAC)) & (MAC > 0) & (MAC < 16))

snpMat.numeric <- snpMat.numeric[row.indicator,]
MAC <- MAC[row.indicator]
snpRanges <- snpRanges[row.indicator,]


## start computing the SNPs within each ATAC-seq peak
## in retrospect, can be simplified using findOverlaps()
chr.list <- gsub('chr', '', unique(data$Chr))
idr.list.with.snp.info <- vector("list", length(chr.list))

foo <- function(i,
                peaks,
                snpRanges,
                snpMat.numeric,
                chr.list) {
  current.chr <- chr.list[i]
  current.peaks <- peaks[seqnames(peaks) == current.chr]
  snp.subsetting.vector <-
    as.logical(seqnames(snpRanges) == current.chr)
  #current.snpMat <- snpMat[snp.subsetting.vector,]
  current.snpMat.numeric <- snpMat.numeric[snp.subsetting.vector,]
  current.snpRanges <- snpRanges[snp.subsetting.vector]
  #current.map <- map[snp.subsetting.vector,]
  
  # preparation of data for a specific chromosome finished
  prev <- 0
  current <- 0
  current.peak.list <- vector("list", length(current.peaks))
  for (j in 1:length(current.peaks)) {
    while (current < length(current.snpRanges)) {
      current <- current + 1
      # check whether current snp is within the range of the peak
      if (start(current.snpRanges)[current] > start(current.peaks)[j] + width(current.peaks)[j]) {
        current <- current - 1
        break
      }
    }
    if (prev < current) {
      prev <- prev + 1
      current.peak.list[[j]] <-
        list('snpMat.numeric' = current.snpMat.numeric[prev:current,], 'snpRanges' =
               current.snpRanges[prev:current])
    }
    prev <- current
  }
  
  # idr.list.with.snp.info[[i]] <- current.peak.list
  current.peak.list
}

idr.list.with.snp.info <-
  mclapply(1:length(chr.list),
           foo,
           peaks,
           snpRanges,
           snpMat.numeric,
           chr.list,
           mc.cores = 20)

idr.list <- unlist(idr.list.with.snp.info, recursive = FALSE)

save(countMatScaled, idr.list, file = 'data-analysis/ATAC-seq/mouse/RData/IDR_SNP_DESeq.RData')

load('data-analysis/ATAC-seq/mouse/RData/IDR_SNP_DESeq.RData')
# countMatScaled: the normalized count matrix of 38,749 differential ATAC-seq master peaks across 16 samples (38749 * 16).
# idr.list: a list with length 38749 containing SNPs within the ATAC-seq master peaks.

snps.per.peak <-
  unlist(lapply(idr.list, function(x) {
    length(x$snpRanges)
  }))
# > summary(snps.per.peak)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00    7.00   12.00   16.27   21.00  854.00
# if no SNP then IDRlist should be NULL
nosnp <- unlist(lapply(idr.list, is.null))
# > sum(nosnp)
# [1] 1883


## Permutation test for testing the significance of the model:
## normalized ATAC-seq signal ~ genotype + (gender)
onePeak.perm.wG <- function(pNo) {
  y <- countMatScaled[pNo, ]
  gender <- as.factor(rep(c('F', 'M'), 8))
  out <- NULL
  if (!nosnp[pNo]) {
    # snps.per.peak[pNo] > 0
    if (snps.per.peak[pNo] > 1) {
      x <- idr.list[[pNo]]$snpMat.numeric
      for (sNo in 1:nrow(x)) {
        geno.x <- rep(x[sNo, ], each = 2)
        m1 <- lmp(y ~ as.factor(geno.x) + gender,  perm = 'Exact')
        out[sNo] <- summary(m1)$coefficients[2, 3]
      }
    }
    else{
      # snps.per.peak[pNo] == 1
      x <- idr.list[[pNo]]$snpMat.numeric
      geno.x <- rep(x, each = 2)
      m1 <- lmp(y ~ as.factor(geno.x) + gender,  perm = 'Exact')
      out[1] <- summary(m1)$coefficients[2, 3]
    }
  }
  return(out)
}

## Compute the results for all peaks
numCores <- detectCores()
# makeCluster(numCores-1)
lmOut.wG <-
  mclapply(1:nrow(countMatScaled), onePeak.perm.wG, mc.cores = numCores)
# lmOut.wG.2 <- mclapply(1:50, onePeak.perm.wG, mc.cores = numCores)
# > sum(unlist(lapply(lmOut.wG, is.null)))
# [1] 1883

min.pval.wG <-
  lapply(lmOut.wG, function(x) {
    ifelse(!is.null(x), min(x), NA)
  })
min.pval.wG <- unlist(min.pval.wG)

## Apply multiple testing corrections
# min.pval.wG <- min.pval.wG[!is.na(min.pval.wG)]
# length(min.pval.wG)
a.min.pval.wG <- p.adjust(min.pval.wG, method = 'BH')
# a.min.pval.wG.2 <- p.adjust(min.pval.wG[!is.na(min.pval.wG)], method = 'BH')
# > sum(a.min.pval.wG[!is.na(a.min.pval.wG)] < 0.05)
# [1] 22200

a.min.pval.wG[is.na(a.min.pval.wG)] <- 1
sum(a.min.pval.wG < 0.05)
# 22200

# the SNP effect indicator
snp.effect.indicator <- a.min.pval.wG < 0.05

min.pval.wG <- min.pval.wG + 0.001
idr.list <- idr.list[snp.effect.indicator]
a.min.pval.wG <- a.min.pval.wG[snp.effect.indicator]
min.pval.wG <- min.pval.wG[snp.effect.indicator]
countMatScaled <- countMatScaled[snp.effect.indicator,]
countMat <- countMat[snp.effect.indicator,]
lmOut.wG <- lmOut.wG[snp.effect.indicator]
min.pval.wG <- min.pval.wG[snp.effect.indicator]
snps.per.peak <- snps.per.peak[snp.effect.indicator]
peaks <- peaks[snp.effect.indicator]
data <- data[snp.effect.indicator,]
rownames(data) <-
  rownames(countMat) <- rownames(countMatScaled) <- 1:22200

## We extract the most significant SNPs for each peak
min.pval.wG <- min.pval.wG + 0.001 # allow for a small offset

## in retrospect, can be optimized using findOverlaps()
for (i in 1:22200) {
  if (snps.per.peak[i] > 1) {
    idr.list[[i]]$snpMat.numeric <-
      idr.list[[i]]$snpMat.numeric[lmOut.wG[[i]] <= min.pval.wG[i],]
    idr.list[[i]]$snpRanges <-
      idr.list[[i]]$snpRanges[lmOut.wG[[i]] <= min.pval.wG[i]]
    lmOut.wG[[i]] <- lmOut.wG[[i]][lmOut.wG[[i]] <= min.pval.wG[i]]
  }
}

#idr.list.sig.snp <- mclapply(1:22200, foo, mc.cores = 20)
snps.per.peak <-
  unlist(lapply(idr.list, function(x) {
    length(x$snpRanges)
  }))
sum(snps.per.peak)
# [1] 47062

whichpeak <- rep(1:22200, snps.per.peak)
# > length(whichpeak)
# [1] 47062

# snpRanges <- mclapply(idr.list, function(x) {return(x$snpRanges)}, mc.cores = 20)
# snpRanges <- do.call('c', snpRanges)
# snpMatrix <- mclapply(idr.list, function(x) {return(x$snpMat.numeric)}, mc.cores = 20)
# snpMatrix <- do.call('rbind', snpMatrix)
# snpPvalue <- unlist(lmOut.wG)

## ----------------------------------------------------------
## 47,062 local-ATAC-QTLs in 22,200 peaks
## With data processing, the results are conveniently saved as
## 'data-analysis/ATAC-seq/mouse/RData/ATAC-QTL.RData'
## ----------------------------------------------------------
