library(ggplot2)
library(grid)
library(gridExtra)
library(ggrepel)
library(viridis)

getwd()
setwd("/home/uibk/c7701178/scratch/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0")

args <- commandArgs(trailingOnly=TRUE)

#Load the covariance matrix
#cov <- as.matrix(read.table("pcangsd/PCAngsd_LDpruned_angsd_GL2_minInd92_maf0.05_minDepth184_maxDepth8656.cov", header = F))
cov <- as.matrix(read.table(args[1], header = F))

#We will also add a column with population assingments
#sample_data <- read.csv("list/samples184_metadata.tsv", header = T, sep='\t')
sample_data <- read.csv(args[2], header = T, sep='\t')
sample_data$sample_id


# IDs
pop <- as.vector(sample_data$sample_id)
pop
mme.pca <- eigen(cov) #perform the pca using the eigen function
eigenvectors = mme.pca$vectors #extract eigenvectors
pca.vectors = as.data.frame(cbind(pop, eigenvectors)) #combine with our population assignments
df = type.convert(pca.vectors)
head(df)


# groups(species)

groups <- as.vector(sample_data$species)
groups
mme.pca <- eigen(cov) #perform the pca using the eigen function
eigenvectors = mme.pca$vectors #extract eigenvectors

pca.vectors2 = as.data.frame(cbind(groups, eigenvectors)) #combine with our population assignments
df2 = type.convert(pca.vectors2)
head(df2)

pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
#percentage <- round((mme.pca$values/pca.eigenval.sum)*100, 2)
#percentage <- paste( colnames(df), "(", paste( as.character(percentage), "%", ")", sep="") )
varPC1 <- round((mme.pca$values[1]/pca.eigenval.sum)*100, 2) #Variance explained by PC1
varPC2 <- round((mme.pca$values[2]/pca.eigenval.sum)*100, 2) #Variance explained by PC2
varPC3 <- round((mme.pca$values[3]/pca.eigenval.sum)*100, 2) #Variance explained by PC3
varPC4 <- round((mme.pca$values[4]/pca.eigenval.sum)*100, 2) #Variance explained by PC4

PC1 <- paste( "PC1", "(", paste( as.character(varPC1), "%", ")", sep=" ") )
PC2 <- paste( "PC2", "(", paste( as.character(varPC2), "%", ")", sep=" ") )
PC3 <- paste( "PC3", "(", paste( as.character(varPC3), "%", ")", sep=" ") )
PC4 <- paste( "PC4", "(", paste( as.character(varPC4), "%", ")", sep=" ") )

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

# Colour palette with black:
cbp2 <- c("#009E73", "#000000", "#D55E00", "#56B4E9", "#E69F00",
          "#F0E442", "#0072B2", "#CC79A7")


# Open a pdf file
pdf(args[3])

## plot PCA with IDs
pca = ggplot(data = df, aes(x=V2, y=V3, label=pop)) + geom_point(pch=16, size=4.5) + ggtitle("PC1 vs PC2 (allSNPs)") + geom_text()
pca = pca + xlab(PC1) + ylab(PC2) + theme + scale_color_manual(values = cbp2)
pca

## PC1vsPC2
pca = ggplot(data = df2, aes(x=V2, y=V3, color=groups)) + geom_point(pch=16, size=4.5) + ggtitle("PC1 vs PC2 (groups, allSNPs)")
pca = pca + xlab(PC1) + ylab(PC2) + theme + scale_color_manual(values = cbp2)
pca

## PC1vsPC3
pca = ggplot(data = df2, aes(x=V2, y=V4, color=groups)) + geom_point(pch=16, size=4.5) + ggtitle("PC1 vs PC3 (groups, allSNPs)")
pca = pca + xlab(PC1) + ylab(PC3) + theme + scale_color_manual(values = cbp2)
pca

## PC2vsPC3

pca = ggplot(data = df2, aes(x=V3, y=V4, color=groups)) + geom_point(pch=16, size=4.5) + ggtitle("PC2 vs PC3 (groups, allSNPs)")
pca = pca + xlab(PC2) + ylab(PC3) + theme + scale_color_manual(values = cbp2)
pca

dev.off()
