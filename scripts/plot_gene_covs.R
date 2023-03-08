rm(list=c(ls()))
setwd("/home/uibk/c7701178/scratch/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0/")
getwd()

##  This R script creates a plot of read coverage of reads based on bedtools genomecov -bga output
# Fileformat is in bedgraph format: 4 columns (chromosome, start, stop, score)
# chromosom     start   end     coverage
# chrom1        0       1       81
# chrom1        1       2       85
# ...
# chrom3        18      19      103
# chrom33       19      20      108
# ...

args <- commandArgs(trailingOnly=TRUE)

# Code modified from https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md
#cov = read.table('bedtools/AMM66.bb.RAPID.dedup.genomecov.bed')
#cov <- read.table('snakemake@input[[1]]')
cov <- read.table(args[1])
head(cov)

names <- list(args[1])

# extract the genome-wide (i.e., no the per-chromosome)   histogram entries
gcov = cov[cov[,1] == 'genome',]
head(gcov)
gcov

# Open a pdf file
#pdf("bedtools/plots/AMM66.pdf")
#pdf("bedtools/plots/snakemake@output[[1]].pdf")
pdf(args[2]) 


# Open a jpeg file
#jpeg("bedtools/plots/AMM66.jpg", width = 350, height = 350)
#jpeg("bedtools/plots/snakemake@output[[1]].jpg", width = 350, height = "350")

# plot a density function for the genome-wide coverage
plot(gcov[1:51,2], gcov[1:51,5], type='h', col='darkgreen', lwd=3, xlab="Depth", ylab="Fraction of genome at depth", main=names)
axis(1,at=c(1,5,10,15,20,25,30,35,40,45,50))


# Create a cumulative distribution from the "raw" hist 
# (truncate at depth >=50)
gcov_cumul = 1 - cumsum(gcov[,5])

# Create a plot of the CDF
plot(gcov[2:51,2], gcov_cumul[1:50], col='darkgreen', type='l', lwd=3, xlab="Depth", ylab="Fraction of genome >= depth", ylim=c(0,1.0), main=names)
axis(1,at=c(1,5,10,15,20,25,30,35,40,45,50))
axis(2,at=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))

dev.off()

# extract chromosom covergae
#dgal1_cov <- subset(cov, V1 == 'dgal1')
#head(dgal1_cov)

# plot a density function for dgal1
#plot(dgal1_cov[1:51,2], dgal1_cov[1:51,5], type='h', col='darkgreen', lwd=3, xlab="Depth", ylab="Fraction of dgal1 at depth")
#axis(1,at=c(1,5,10,15,20,25,30,35,40,45,50))


