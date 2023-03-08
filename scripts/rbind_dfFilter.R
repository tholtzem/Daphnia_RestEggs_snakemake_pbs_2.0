
library(data.table)
library(dplyr)
library(purrr)
library(stringr)

args <- commandArgs(trailingOnly=T)

setwd("/scratch/uibk/c7701178/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0/depth/stats/")

# list AND read files in directory with specific pattern
df <- list.files(path = "/scratch/uibk/c7701178/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0/depth/stats/", pattern = "*_depthFilter.list") %>% map_df(~fread(.))
df

# list file names
pops <- list.files(path = "/scratch/uibk/c7701178/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0/depth/stats/", pattern = "*_depthFilter.list")
pops <- str_replace(pops, "_depthFilter.list", "")
pops


df$pops <- pops
df

# write new df to file
write.table(df, file = args[1], sep = "\t", quote = F, row.names=F)
#write.table(df, "depthFilter.list", sep = '\t', quote = F, row.names=F)
