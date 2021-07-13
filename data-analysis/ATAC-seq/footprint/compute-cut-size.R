## ----------------------------------------------------------
##   Compute Tn5 insertions for all ATAC-seq samples
## ----------------------------------------------------------


## 5' end only - Tn5 insertions
## [strain-name]_ATAC_cutsize.RData
## objects
## cvglist: the coverage object for both strands
##  cvglist$'+' or cvglist$'-'
##  each strand is a integer RLEList contains 21 chromosomes
##  cvglist$'+'$chr1


## Prepare cut size files for BAM files ####
library(data.table)
library(GenomicAlignments)
library(BSgenome.Mmusculus.UCSC.mm10)
library(parallel)

genome <- BSgenome.Mmusculus.UCSC.mm10
seqlev <- seqlevels(genome)[1:21]
seqlen <- seqlengths(genome)[1:21]
bamdir <- 'bam/'
bamfileList <- list.files(path = bamdir, pattern = '\\.bam$')
bamfileList <- paste0(bamdir, '/', bamfileList)
strains <- c('129', 'AJ', 'B6', 'CAST', 'NOD', 'NZO', 'PWK', 'WSB')

## read in bam file with input seqlev specified by users
which <- as(Seqinfo(seqlev, seqlen), "GRanges")

## this parameter restrict the genome regions of the aligned reads we want.
param <- ScanBamParam(which = which)

for (i in 1:8) {
  index <- bamfileList[i]
  bamIn <-
    mapply(function(.b, .i)
      readGAlignments(.b, .i, param = param),
      bamfileList[i],
      index,
      SIMPLIFY = FALSE)
  
  bamIn <- lapply(bamIn, as, Class = "GRanges")
  
  if (!is(bamIn, "GRangesList"))
    bamIn <- GRangesList(bamIn)
  bamIn <- unlist(bamIn)
  seqlevelsStyle(bamIn) <- seqlevelsStyle(genome)
  names(bamIn) <- NULL
  
  ## keep 5' end as cutting sites
  bamIn <- promoters(bamIn, upstream = 0, downstream = 1)
  ## also count the 3'end as cutting sites
  # bamIn <- c(promoters(bamIn, upstream=0, downstream=1),
  #            flank(bamIn, width = 1, start = F))
  
  libSize <- length(bamIn) # the total number of cuts
  coverageSize <-
    sum(as.numeric(width(reduce(
      bamIn, ignore.strand = TRUE
    ))))
  libFactor <- libSize / coverageSize
  # this should be larger than one, the average cut number for cut size.
  
  ## split into positive strand and negative strand
  bamIn <- split(bamIn, strand(bamIn))
  ## get coverage
  cvglist <- sapply(bamIn, coverage)
  cvglist <- cvglist[c("+", "-")]
  cvglist <- lapply(cvglist, function(.ele)
    .ele[names(.ele) %in% seqlev])
  
  
  save(
    cvglist,
    libSize,
    coverageSize,
    file = paste0(
      'data-analysis/ATAC-seq/mouse/RData/',
      strains[i],
      '_ATAC_cutsize.RData'
    )
  )
}
