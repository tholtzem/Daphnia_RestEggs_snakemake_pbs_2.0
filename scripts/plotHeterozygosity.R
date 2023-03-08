getwd()
setwd("/home/uibk/c7701178/scratch/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0/")

args <- commandArgs(trailingOnly=TRUE)


a <-scan(args[1])
a[2]/sum(a)
