##### Tania Holtzem
##### February 2022


#------------------------------------------------------------------------------------------------
#     1) PLOTTING ALL SITES QUALITY MEASUREMENTS
#------------------------------------------------------------------------------------------------
#Enter the path to the folder with the data file
#setwd("/mnt/data/snakemake_mach2/daphnia_analysis_smk")
getwd()

args <- commandArgs(trailingOnly=TRUE)

#Read the file into R, interprete "." as missing data
#values <- read.table("ancestry/parents_hybrids_INFODP.txt", na.strings=c("."))
values <- read.table(args[1], na.strings=c("."))
values$V1

#Add a header line to the data frame
colnames(values) <- c("DP")

#Check if the column headers were added correctly by printing the first lines
head(values)

summary(values)

meanDP <- round(mean(values$DP))
HengLimax <- round(meanDP + 3*sqrt(meanDP))
IQRmax <- round(1.5*IQR(values$DP))
twoxmean <- round(2*meanDP)


#Plot histograms of all parameters of interest for "LG05"
#Before saving all plots in one PDF file, execute the plotting commands individually
#for each quality measurement to adjust the plotting area. Replace the "XXX" for xlim and ylim
#with meaningful values for plotting.

pdf(args[2], height=25, width=10) #save them in one PDF file
par(mfrow=c(2,2)) #make a multi-paneled plotting window with 7 rows and 2 columns
hist(values$DP, xlab = "DP", main=paste0("LC", "_DP org."))
hist(values$DP, xlab = "DP", xlim = c(0,40000), ylim = c(0,100000), breaks=1000,  main=paste0("LC", "_DP scaled"))
dev.off() #close the plot

#To decide on a depth filtering theshold, you may be interested in the percentage of variants with an DP over a certain treshold 
proportion_DP1 <- ((sum((values$DP) > 1))*100)/nrow(values) #95.7% of sites have a DP > 1
proportion_DPu1 <- ((sum((values$DP) < 1))*100)/nrow(values) #3.1% of sites have a DP < 1
proportion_DP10 <- ((sum((values$DP) > 10))*100)/nrow(values) #91.3% of variants have a DP > 10
proportion_DP_HengLi <- ((sum((values$DP) > HengLimax))*100)/nrow(values) #53.1% of sites have a DP > maxDepth_Filter_HengLi
proportion_DP_IQR1.5 <- ((sum((values$DP) > IQRmax))*100)/nrow(values) #3.3% of sites have a DP > maxDepth_Filter_IQR1.5
proportion_DP_2xmean <- ((sum((values$DP) > twoxmean))*100)/nrow(values) #4.7% of sites have a DP > maxDepth_Filter_2xmean


# write to file
output <- data.frame(meanDP, HengLimax, IQRmax, twoxmean, proportion_DPu1, proportion_DP1, proportion_DP10, proportion_DP_HengLi, proportion_DP_IQR1.5, proportion_DP_2xmean)
colnames(output) <- c("meanDP", "HengLimax", "IQRmax", "twoxmean", "proportion_DPu1", "proportion_DP1", "proportion_DP10", "proportion_DP_HengLi", "proportion_DP_IQR1.5", "proportion_DP_2xmean")
output <- rbind(output)
write.table(output, file = args[3], sep = "\t",
            row.names = FALSE, quote = FALSE)
