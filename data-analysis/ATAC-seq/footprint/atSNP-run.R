## ----------------------------------------------------------
##           Large-scale computing with atSNP
## ----------------------------------------------------------

library(atSNP)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BiocParallel)

## We suggest running atSNP analysis on a HPC
## and speed up the program with batches of input SNPs (e.g. by chromosomes)
## It takes a few days to run
args <- commandArgs(TRUE)

## Please refer to the vignettes of atSNP R Package for details
input.file <- as.character(args[1])
n.cores <- as.integer(args[2])
output.file <- as.character(args[3])

## the JASPAR motif library
load('data-analysis/ATAC-seq/mouse/RData/motif_lib_total.RData')

snpInfo <-
  LoadSNPData(
    input.file,
    genome.lib = "BSgenome.Mmusculus.UCSC.mm10",
    half.window.size = 30,
    default.par = FALSE,
    mutation = FALSE
  )


motif_library <- motif_lib_total

motif_score <-
  ComputeMotifScore(motif_library, snpInfo, ncores = n.cores)
motif_score$motif.scores <-
  ComputePValues(
    motif.lib = motif_library,
    snp.info = snpInfo,
    motif.scores = motif_score$motif.scores,
    ncores = n.cores
  )

save(motif_score, file = output.file)
## the orignal output from atSNP can be huge
## atSNP_output.RData 2.6 GB