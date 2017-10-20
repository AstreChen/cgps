library(gage)
library(Biobase)
library(limma)
library(EnrichmentBrowser)
library(edgeR)

######
args <- commandArgs(TRUE)
if (length(args) != 4){
  stop("Rscript test.R expfile phefile datatype=[ma/rseq] output_directory \n")
}
expf <- args[1]
phef <- args[2]
dtype <- args[3]  # read count
#spe <- args[4]
outdir <- args[4]
spe <- 'hsa'

source('./individual_methods.R')

expdata <- read.exp(expf)  # return matrix
pdata <- read.phe(phef) # return AnnotateDF
write.csv(expdata,paste0(outdir,'./test.csv'))