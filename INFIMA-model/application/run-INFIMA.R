## ----------------------------------------------------------
##         fit the INFIMA model on the real data
## ----------------------------------------------------------


library(INFIMA)

setwd('INFIMA-model/application/')
load('inputs.RData')

## dt1: DO mouse eQTL results
## dt2: founder RNA-seq count matrix
## dt3: SNP information (with footprint annotation) and the local ATAC-seq signal
## For more information about the input data format, please refer to the vignette or ?raw_input_data
raw_data <- raw_input_data(dt1, dt2, dt3)
model_data <- model_input_data(raw_data) # 40 seconds

## Dirichlet prior for the INFIMA model
prior <- compute_prior(raw_data, model_data, verbose = F)

## Model fitting
infima <- model_fitting(model_data, prior, verbose = F)

summary(infima)
# INFIMA converges after 15 iterations!
#   Null generative model parameters:
#   a0 = (0.64788, 0.30744, 0.04468)
# b0 = (0.54889, 0.40517, 0.04594)
# Causal generative model parameters:
#   a1 = (8.611e-01, 1.389e-01, 2.047e-28)
# b1 = (8.012e-01, 1.988e-01, 6.953e-12)
# gamma = 0.6196 

## After fitting the INFIMA model, we focus on the DO genes containing a causal SNP
infima_results <- snp_link_gene(infima, raw_data, model_data, prior, fdr = 0.1, cum.pprob = 0.85, cred.set = 0.5, verbose = F)

## Obtain the fine-mapping results
## Link mouse SNPs to target genes
results <- as.data.frame(infima_results)

## An example input data query
input_query <- query_input_data(infima_results, snp_id = 'rs51076312', ensembl = 'ENSMUSG00000037995', qtl_marker = '1_172713578')

save(raw_data, model_data, prior, infima, infima_results, results,
     file = 'results.RData') # 


#### adaptive liftover results linking GWAS to effector genes

library(data.table)
library(GenomicAlignments)
library(parallel)
library(INFIMA)

load('examples.RData')
# the final 499 cases validated by INFIMA
# dt, atac.mm10, res, gwas.peak, gwas.snps, gwas.df, gwas.mm10.peak.infima

# GWAS - mouse ATAC-seq peaks (get the causal SNPs) 
# get SNP rsid, marker id, ensembl id

# remove the infima results where 

setwd('INFIMA-model/application/RData')
load('inputs.RData')
load('results.RData') 
# raw_data, model_data, prior, infima, infima_results, results (as.data.table(infima_results))

#### Step 1: filter DO genes ####
x <- infima_results
infima <- x$infima
res <- x$results
res <- res[sapply(res, function(x)
  ! is.null(x))] # 5720 DO genes after filtering

genes.all <- unique(unlist(sapply(strsplit(dt$genes, '-'), function(x) x)))
res <- res[sapply(res, function(x) toupper(x$input_data$do.eqtl$symbol) %in% genes.all)]
# only 212 DO genes after filtering

#### Step 2: filter SNPs
peaks <- GRanges(seqnames = dt$`mm-ATAC-chr`,
                 IRanges(start = dt$`mm-ATAC-start`,
                         end = dt$`mm-ATAC-end`),
                 strand = '*',
                 name = paste0(dt$`mm-ATAC-chr`, ':', dt$`mm-ATAC-start`, '-', dt$`mm-ATAC-end`))
peaks <- reduce(peaks) # only 186 peaks
snps <- GRanges(seqnames = dt3$chr,
                IRanges(start = dt3$snp_pos,
                        end = dt3$snp_pos),
                strand = '*',
                name = dt3$snp_id)
ii <- findOverlaps(snps, peaks)
snps <- snps[queryHits(ii)] # only 304 SNPs
snps.all <- snps$name


###################################################

infima_res <- list(infima = infima, results = res)
class(infima_res) <- 'infima_results'

library(INFIMA)
df <- as.data.frame(infima_res)
df <- df[snp_id %in% snps$name] # 547 input querys
# > length(unique(df$symbol))
# [1] 212
library(doParallel)
n_cores <- 30
registerDoParallel(cores = n_cores)

input_queries <- foreach(i = 1:nrow(df)) %dopar% {
query_input_data(infima_res,                              
                 snp_id = df$snp_id[i],
                 ensembl = df$ensembl[i],
                 qtl_marker = df$qtl_marker[i])
}

# reduce the size
input_queries_small <- vector('list', length(input_queries))
for(i in 1:length(input_queries_small)){
  input_queries_small[[i]] <- list(Y = as.numeric(input_queries[[i]]$Y),
                                   Y.t = as.numeric(input_queries[[i]]$Y.t),
                                   A = as.numeric(input_queries[[i]]$A),
                                   A.t = as.numeric(input_queries[[i]]$A.t),
                                   B.avg = as.numeric(input_queries[[i]]$B.avg),
                                   B.t = as.numeric(input_queries[[i]]$B.t),
                                   D = as.numeric(input_queries[[i]]$D),
                                   E.t = as.numeric(input_queries[[i]]$E.t),
                                   genotype = as.numeric(input_queries[[i]]$snpData[, `129-genotype`:`WSB-genotype`])
                                   )
}

save(dt, df, input_queries_small, file = 'INFIMA-shinyapp-data.RData')
