## ----------------------------------------------------------
##           get footprints from atSNP
## ----------------------------------------------------------

## Please note that the putative TFBS with maximum likelihood for reference allele
## and alternative allele may not be exactly the same from atSNP

library(data.table)
library(parallel)
library(GenomicAlignments)
library(BiocParallel)

args <- commandArgs(TRUE)
n.cores <- as.integer(args[1])

setwd('data-analysis/ATAC-seq/mouse/RData/')
load('atSNP_output.RData') # the combined output results from atSNP
snp.info <-
  fread('data-analysis/ATAC-seq/mouse/footprint/atSNP_input.txt')
motif.scores <- motif_score$motif.scores

# system.time(
#   index <- unlist(bplapply(1:nrow(motif.scores),
#                 function(x) {if(x%%10000 == 0) {print(x)}; grep(motif.scores$snpid[x], snp.info$snpid)}, BPPARAM = SnowParam(workers = n.cores, type = "SOCK")))
# )

index <- unlist(mclapply(1:nrow(motif.scores),
                         function(x) {
                           if (x %% 1000 == 0) {
                             print(x)
                           }
                           grep(motif.scores$snpid[x], snp.info$snpid)
                         }, mc.cores = n.cores))

# TOO SLOW!!!! about 16 hours


seqname <- as.character(snp.info[index]$chr)
snp <- snp.info[index]$snp
width <- motif.scores$motif_len

window <- 50
start.ref <- snp + motif.scores$ref_start - 31 - window
end.ref <- snp + motif.scores$ref_end - 31 + window
start.snp <- snp + motif.scores$snp_start - 31 - window
end.snp <- snp + motif.scores$snp_end - 31 + window

save(start.ref, end.ref, start.snp, end.snp, seqname, width, file = 'footprints_from_atSNP.RData') # 1.2 GB
