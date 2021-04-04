## ----------------------------------------------------------
##                  run MatrixEQTL
## ----------------------------------------------------------


## Matrix eQTL by Andrey A. Shabalin
## http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
## Be sure to use an up to date version of R and Matrix eQTL.

library(MatrixEQTL)

## Location of the input data files.
base.dir = 'data-analysis/RNA-seq/MatrixEQTL'
setwd(base.dir)

## Settings
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR
# modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste(base.dir, "SNP.txt", sep = "")

snps_location_file_name = paste(base.dir, "snpsloc.txt", sep = "")


# Gene expression file name
expression_file_name = paste(base.dir, "GE.txt", sep = "")

gene_location_file_name = paste(base.dir, "geneloc.txt", sep = "")


# Covariates file name
# Set to character() for no covariates
# covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="");
covariates_file_name = character()

# Output file name
output_file_name_cis =  paste(base.dir, "cisEQTL.txt", sep = "")

output_file_name_tra = paste(base.dir, "transEQTL.txt", sep = "")


# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-5

pvOutputThreshold_tra = 1e-5


# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric()

# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6


## Load genotype data

snps = SlicedData$new()

snps$fileDelimiter = "\t"
# the TAB character
snps$fileOmitCharacters = "NA"
# denote missing values;
snps$fileSkipRows = 1
# one row of column labels
snps$fileSkipColumns = 1
# one column of row labels
snps$fileSliceSize = 2000
# read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)


## Load gene expression data

gene = SlicedData$new()

gene$fileDelimiter = "\t"
# the TAB character
gene$fileOmitCharacters = "NA"
# denote missing values;
gene$fileSkipRows = 1
# one row of column labels
gene$fileSkipColumns = 1
# one column of row labels
gene$fileSliceSize = 2000
# read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)


## Load covariates

cvrt = SlicedData$new()

cvrt$fileDelimiter = "\t"
# the TAB character
cvrt$fileOmitCharacters = "NA"
# denote missing values;
cvrt$fileSkipRows = 1
# one row of column labels
cvrt$fileSkipColumns = 1
# one column of row labels
if (length(covariates_file_name) > 0) {
  cvrt$LoadFile(covariates_file_name)
  
}

## Run the analysis
snpspos = read.table(snps_location_file_name,
                     header = TRUE,
                     stringsAsFactors = FALSE)

genepos = read.table(gene_location_file_name,
                     header = TRUE,
                     stringsAsFactors = FALSE)


me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)


# unlink(output_file_name_tra);
# unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')

cat('Detected local eQTLs:', '\n')

cis <- me$cis$eqtls
#show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n')

trans <- me$trans$eqtls
#show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values

setwd(base.dir)

pdf('./qqplot.pdf')
plot(me)
dev.off()

save(me, cis, trans, file = 'results.RData')


# plot histogram of p-values
output_file_name_cis =  paste(base.dir, "/output/cisEQTL2.txt", sep = "")

output_file_name_tra = paste(base.dir, "/output/transEQTL2.txt", sep = "")

me2 = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = 100,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)


# pdf('./histogram.pdf')
# plot(me2)
# dev.off()



# # eQTL and DO-eQTL genes summary

# load('data-analysis/RNA-seq/RData/DO_eQTL.RData')
# cis.genes <- unique(as.character(do.eqtl[cis == TRUE]$symbol))
# trans.genes <- unique(as.character(do.eqtl[cis == FALSE]$symbol))

# cisEQTL <- fread(file = 'cisEQTL.txt')
# cisEQTL$Symbol <- unlist(lapply(1:nrow(cisEQTL), function(x) {strsplit(cisEQTL$gene[x],'_')[[1]][2]}))
# cisEQTL$gene <- unlist(lapply(1:nrow(cisEQTL), function(x) {strsplit(cisEQTL$gene[x],'_')[[1]][1]}))
# cistotal <- length(unique(cisEQTL$Symbol))
# cisEQTL <- cisEQTL[Symbol %in% cis.genes]
# meta.cis <- list(total = cistotal, overlap = length(unique(cisEQTL$Symbol)))

# 6412 -> 5097
# nrow = 82188

# transEQTL <- fread(file = 'transEQTL.txt')
# transEQTL$Symbol <- unlist(lapply(1:nrow(transEQTL), function(x) {strsplit(transEQTL$gene[x],'_')[[1]][2]}))
# transEQTL$gene <- unlist(lapply(1:nrow(transEQTL), function(x) {strsplit(transEQTL$gene[x],'_')[[1]][1]}))
# transtotal <- length(unique(transEQTL$Symbol))
# transEQTL <- transEQTL[Symbol %in% trans.genes]
# meta.trans <- list(total = transtotal, overlap = length(unique(transEQTL$Symbol)))

# 12020 -> 8708
# nrow = 56261052
# save(cisEQTL, transEQTL, meta.cis, meta.trans, file = 'eQTL_DO_genes.RData')

# length(unique(union(transEQTL$Symbol, cisEQTL$Symbol)))
# 10160