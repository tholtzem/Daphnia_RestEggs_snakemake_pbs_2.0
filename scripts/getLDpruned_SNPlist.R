#rm(list=ls())
detach(package:stats)
library(tidyverse)
#getwd()
#setwd("/mnt/data/snakemake_mach2/Daphnia_RestEggs_snakemake_pbs/")

basedir="/home/uibk/c7701178/scratch/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0/"

args <- commandArgs(trailingOnly=TRUE)

#pruned_position <- read_lines(paste0(basedir, "ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.unlinked.id"))
pruned_position <- read_lines(args[1])
pruned_position2 <- sub('.*:', '', pruned_position) %>%
  as.integer()

#pruned_snp_list <- read_tsv(paste0(basedir, "angsd/angsd_LC_GL2_cutoff_maf0018_nInd55.mafs.gz")) %>%
pruned_snp_list <- read_tsv(args[2])
pruned_snp_list2 <- dplyr::select(pruned_snp_list, 1:4)
pruned_snp_list3 <- pruned_snp_list2 %>% filter(pruned_snp_list2$position %in% pruned_position2)

#write_tsv(pruned_snp_list, "list/LDpruned_snps.list", col_names = F)
write_tsv(pruned_snp_list3, args[3], col_names = F)
