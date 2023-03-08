# This scriptcreates diverse depth plots and calculates min and max depth filters
library(tidyverse)

getwd()
setwd("/home/uibk/c7701178/scratch/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0/")

args <- commandArgs(trailingOnly=T)

stats <- read.table(file = args[1], sep = '\t', header = TRUE)

sets <- as.list(args[2])
#sets <- args[2]

# subset dataframe stats
df1 <- subset(stats, mean_depth > 1)
df1

dfu1 <- subset(stats, mean_depth < 1)
dfu1$bamfile

df10 <- subset(stats, mean_depth > 10)
df10$bamfile

df20 <- subset(stats, mean_depth > 20)
df20$bamfile

df3 <- subset(stats, mean_depth > 3)
df3$bamfile


## Bar plot of mean read depth per sample ##
pdf(paste0("depth/plots/", sets, "_depth_hist_dfAll.pdf"))
barplot(stats$mean_depth, names=stats$id,las=2,cex.names=0.40, ann = F)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Mean read depth", line = 3)
title("Mean read depth per sample (all samples)")
dev.off()


## Bar plot for samples with a mean depth greater than 1 ##
pdf(paste0("depth/plots/", sets, "_depth_hist_df1.pdf"))
barplot(df1$mean_depth, names=df1$id,las=2,cex.names=0.40, ann = F)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Mean read depth", line = 3)
title("Mean read depth per sample (df 1)")
dev.off()

## Bar plot for samples with a mean depth greater than 10 ##
pdf(paste0("depth/plots/", sets, "_depth_hist_df10.pdf"))
barplot(df1$mean_depth, names=df1$id,las=2,cex.names=0.40, ann = F)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Mean read depth", line = 3)
title("Mean read depth per sample (df 10)")
dev.off()


## Check distributions of all samples using boxplots
#R's boxplot function uses the standard rule to indicate an observation as a potential outlier
#if it falls more than 1.5 times the IQR (Inter-Quartile Range, calculated as Q3-Q1) below Q1 or above Q3.
#The potential outliers are plotted with circles and the Whiskers (lines that extend from Q1 and Q3 typically to the minimum and maximum) are shortened to only go as far as observations that are within 1.5*IQR of the upper and lower quartiles. 
#The box part of the boxplot is a box that goes from Q1 to Q3 and the median is displayed as a line somewhere inside the box

# Distribution of mean read depth for all samples 
pdf(paste0("depth/plots/", sets, "_depth_boxplot_dfAll.pdf"))
boxplot(stats$mean_depth)
title("Distribution of mean read depth (all samples)")
dev.off()

# Distribution of mean read depth for depth > 1
pdf(paste0("depth/plots/", sets, "_depth_boxplot_df1.pdf"))
boxplot(df1$mean_depth)
title("Distribution of mean read depth")
dev.off()

# Distribution of mean read depth for depth > 10
#pdf(paste0("depth/plots/", sets, "_depth_boxplot_df10.pdf"))
#boxplot(df10$mean_depth)
#title("Distribution of mean read depth")
#dev.off()

# Summary of good samples
df1_sum <- summary(df1$mean_depth)
df1_sum
#mean
meanDepth_allInd <- mean(df1$mean_depth)
meanDepth_allInd
#median
medianDepth_allInd <- median(df1$mean_depth)
medianDepth_allInd
#IQR
IQRT<-IQR(df1$mean_depth)
IQRT

# number of individuals
N = nrow(df1)
N

#1.5 times the interquartile range from the mean (median) were excluded as they are expected to be enriched for paralogs
## mean
MaxDepth1 = (meanDepth_allInd + IQRT*1.5)*N
MaxDepth1
MinDepth1  = (meanDepth_allInd - IQRT*1.5)*N
MinDepth1
## mediam
MaxDepth2 = (medianDepth_allInd + IQRT*1.5)*N
MaxDepth2
MinDepth2  = (medianDepth_allInd - IQRT*1.5)*N
MinDepth2

# d + 3*sqrt(d), where d=mean coverage
HengLi1_max <- (meanDepth_allInd + 3*sqrt(meanDepth_allInd))*N
HengLi1_max
HengLi2_max <- (medianDepth_allInd + 3*sqrt(meanDepth_allInd))*N
HengLi2_max

HengLi1_min <- (meanDepth_allInd - 3*sqrt(meanDepth_allInd))*N
HengLi1_min
HengLi2_min <- (medianDepth_allInd - 3*sqrt(meanDepth_allInd))*N
HengLi2_min

# alternative approach
## sum mean_depth for all individuals
total <- mean(df1$mean_depth)
total
# sum the square root (SQRT) of the mean depth
SQRT <- sqrt(df1$mean_depth)
SQRT
# sum IQRT for all individuals
iqr <- IQR(df1$mean_depth)
iqr

## Filter1 - Heng Li
maxFilter1 <- mean((total + 3*SQRT))
minFilter1 <- mean((total - 3*SQRT))
## Filter2 - 1.5 times the interquartile range from the mean
maxFilter2 <- (total + iqr*1.5)
minFilter2 <- (total - iqr*1.5)
## Filter3 - 2 times the mean
maxFilter3 <- (2*total)


output <- data.frame(maxFilter1, minFilter1, maxFilter2, minFilter2, maxFilter3, N)
colnames(output) <- c("HengLi_max", "HengLi_min", "IQR_max", "IQR_min", "2x_mean", "Number_of_samples")
output <- rbind(output)
write.table(output, file = paste0("depth/stats/", sets, "_depthFilter.list"), sep = "\t", quote = F, row.names = F)

############################################################
#non-zero

dfNonZero <- summary(df1$mean_depth_nonzero)
dfNonZero

# Bar lot of mean read depth per sample
## samples having zero depth excluded
pdf(paste0("depth/plots/", sets, "_depth_NonZero_hist_df1.pdf"))
barplot(df1$mean_depth_nonzero, names=df1$id,las=2,cex.names=0.40, ann = F)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Mean read depth", line = 3)
title("Mean read depth (non-zero) per sample (df1)")
dev.off()

##########################################################
# write samples having a mean read depth higher than 1 to a text file, without using quotes
#realigned <- paste0(df3$bamfile,".minq20.realigned.bam")
#Extract Characters Before Pattern using sub()
## first pattern, then replacement of the pattern, then string
realignedBAM_df1 <-sub(".depth.gz", "", paste0('realigned/',df1$bamfile))
write.table(realignedBAM_df1, file = paste0("depth/stats/", sets, "_realignedBAM_df1.list"), sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)

# write samples having a mean read depth smaller than 1 (<10) to a text file, without using quotes
# Extract Characters Before Pattern using sub()
## first pattern, then replacement of the pattern, then string
realignedBAM_dfu1 <-sub(".depth.gz", "", paste0('realigned/', dfu1$bamfile))
write.table(realignedBAM_dfu1, file = paste0("depth/stats/", sets, "_realignedBAM_depth_dfu1.list"), sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
