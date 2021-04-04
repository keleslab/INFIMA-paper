## ----------------------------------------------------------
##   Compute footprint depth (FPD) for motif binding sites
## ----------------------------------------------------------

library(data.table)
library(GenomicRanges)
library(parallel)

setwd('data-analysis/ATAC-seq/mouse/RData/')

## We summarized the standard FIMO outputs for B6 ATAC-seq samples and computed the aggregated Footprint Depth (FPD) for each motif (1,316 in total).
## FIMO outputs are huge in space, we processed them and saved them as RData files.
## for each motif
rdatadir <- 'fimoRes/'
setwd(rdatadir)
## FIMO outputs {motif_index}.fimoRes.RData, compute the footprint depth for each grange object first
outdir <- 'fimoResFpd/'


load('B6_ATAC_cutsize.RData')
cvglist.B6 <- cvglist
load('jasparfix.RData')


## input GRange, motif.index
## return the same GRange with new column FPD
ComputeFPDsub <- function(profile.info) {
  # la lb bind rb ra
  la <- profile.info[1]
  lb <- profile.info[2]
  bind <- profile.info[3]
  rb <- profile.info[4]
  ra <- profile.info[5]
  fpda <- 0
  fpdb <- 0
  if (la + ra > 0) {
    fpda <- 1 - 2 * bind / (la + ra)
  }
  if (lb + rb > 0) {
    fpdb <- 1 - 2 * bind / (lb + rb)
  }
  return(max(fpda, fpdb))
}

ComputeFPD <- function(grange, motif.index, cvglist, window = 50) {
  fpd <- NULL
  if (!is.null(grange) && length(grange) > 0) {
    cvgfs <- cvglist[['+']]
    cvgrs <- cvglist[['-']]
    
    half.window <- as.integer(window / 2)
    # 1 - 25, 26 - 50, 51 - (50 + width), (51 + width) - (75 + width), (76 + width) - (100 + width)
    breaks <- c(
      0,
      half.window,
      half.window * 2,
      half.window * 2 + width[motif.index],
      half.window * 3 + width[motif.index],
      half.window * 4 + width[motif.index]
    )
    
    
    start <- start(grange) - window
    end <- end(grange) + window
    seqname <- as.character(seqnames(grange))
    fpd <- rep(0, length(grange))
    
    for (i in 1:length(grange)) {
      profile <-
        cvgfs[[seqname[i]]][start[i]:end[i]] + cvgrs[[seqname[i]]][start[i]:end[i]]
      profile.info <- rep(0, 5)
      for (j in 1:5) {
        profile.info[j] <- mean(profile[(breaks[j] + 1):breaks[j + 1]])
      }
      fpd[i] <- ComputeFPDsub(profile.info)
    }
  }
  
  return(fpd)
}

OneMotifComputeAggregateProfile <- function(i) {
  infile <- paste0('./', i, '.fimoRes.RData')
  outfile <- paste0(outdir, i, '.fimoResFpd.RData')
  load(infile)
  print(i)
  if (!file.exists(outfile)) {
    if (!is.null(gr)) {
      gr$fpd.B6 <- ComputeFPD(gr, i, cvglist.B6)
      gr$fpd.Cast <- ComputeFPD(gr, i, cvglist.Cast)
    }
    save(gr, file = outfile)
  }
}

## Compute FPD for aggregated profiles to get an idea of overall TF binding
## within ATAC-seq peaks
results <-
  mclapply(1:1316, OneMotifComputeAggregateProfile, mc.cores = 20)


## Full atSNP footprint analysis #####

args <- commandArgs(TRUE)
n.cores <- detectCores()

load('footprints_from_atSNP.RData')

## Take 129 strain as an example, we perform the same computations for all
## 8 founder strains: 129, AJ, B6, Cast, NOD, NZO, PWK, WSB

load('129_ATAC_cutsize.RData')
cvglist.129 <- cvglist
start <- start.ref
end <- end.ref

fpd.129 <-
  mclapply(1:length(start), function(x)
    ComputeFPD(x, cvglist.129), mc.cores = n.cores)
fpd.129 <- unlist(fpd.129)
save(fpd.129, file = 'fpd.129.RData')
