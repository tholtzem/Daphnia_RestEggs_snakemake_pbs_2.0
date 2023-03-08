#library(tidyverse) #load the tidyverse package for formatting and plotting
library(viridis)
library(plyr)
library(readr)
library(ggplot2)

#In this script we will plot heterozygosity computed by angsd using the -doSaf with the reference genome providing the ancestral state.
#This was achieved by dividing the number of heterozygous sites (the second entry in the SFS)
#by the total number of sites (first entry + second entry in the SFS)
#to obtain the proportion of heterozygous sites for each genome.


setwd("/mnt/data/snakemake_mach2/Daphnia_RestEggs_snakemake_pbs_2.0/")

#args <- commandArgs(trailingOnly=TRUE)

# population assingments
sample_data <- read.csv("list/samples184_metadata.tsv", header = T, sep='\t')
sample_data$sample_id
sample_data$species

# species
species <- as.vector(sample_data$species)
species

# IDs
ID <- as.vector(sample_data$sample_id)
ID


basedir <- "saf" # Make sure to edit this to match your $BASEDIR
bam_list <- list.files(path = basedir, pattern = "*est.ml", full.names = TRUE)

bam_list


# ldply: Split list, apply function, and return results in a data frame.
data_tsv = ldply(bam_list, scan)

# total number of sites
totalNR_sites <- (data_tsv$V1 + data_tsv$V2)
totalNR_sites

# number of homozygote sites
NR_HOMsites <- data_tsv[1]

# number of heterozygote sites
NR_HETsites <- data_tsv[2]
NR_HETsites

# individual heterozygosity
IND_heterozygosity <- (NR_HETsites/totalNR_sites)

output <- cbind(ID, species, NR_HOMsites, NR_HETsites, IND_heterozygosity)
output

names(output) = c("ID", "species", "NR_HOMsites", "NR_HETsites", "IND_heterozygosity")
options(scipen = 100)
write.table(output,"pop_stats/heterozygosity_angsd.txt", sep ="\t", quote = F)


# plot bar charts of heterozygosity estimates

barplot(t(as.matrix(output$IND_heterozygosity)), names=output$ID, border=NA, las=2, cex.names = 0.5, cex.axis = 0.8, xlab = "Individuals", ylab = "Heterozygosity")


png("pop_stats/heterozygosity_angsd.png", width = 4, height = 4, units = 'in', res = 300)
# barplot
p<-ggplot(data=output, aes(x=ID, y=IND_heterozygosity, color=species, fill=species, label=ID))+
  geom_bar(stat="identity", width=0.5)+
  theme(axis.text = element_text(size = 6))+
  scale_color_manual(values=c("#009E73", "#000000", "#D55E00", "#56B4E9"))
p2 <- p + scale_fill_manual(values=c("#009E73", "#000000", "#D55E00", "#56B4E9"))+
  xlab("Indiviuals")+
  ylab("Heterozygosity")+
  coord_flip()
p2 + theme(panel.background = element_blank(), axis.line = element_line(colour = "grey"))+
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
dev.off()





