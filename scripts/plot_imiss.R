args <- commandArgs(trailingOnly=TRUE)

setwd("/home/uibk/c7701178/scratch/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0/")

#imiss <- read.table("analyses/scikit_allel/data/daphnia_9pops_vars_alleles.imiss", header=FALSE)
imiss <- read.table(args[1], header=FALSE)

#Add a third column containing the percentage of missing data per individual
imiss <- cbind(imiss, V3=((imiss$V2*100)/args[2])) #with cbind we combine "imiss" with a third column "V3=", 
#containing the value in the second column "imiss$V2", devided by the total amount of reads.

#Sort the individuals according to the amount of missing data
imiss <- imiss[order(imiss$V3),]

#Make a barplot of missing data
#pdf("analyses/scikit_allel/data/daphnia_9pops_vars_alleles.imiss.pdf", height=25, width=10) #save the plot in a PDF file
pdf(args[3], height=25, width=10) #save the plot in a PDF file

par(mar=c(3, 5, 3, 1)) #enlargen figure margins (bottom, left, top, and right)
barplot(imiss$V3, main="Missingness per individual (total sites 1297802)", horiz=T, xlim=c(0,100), names.arg=imiss$V1,
        cex.names=0.5, las=1, space=5) #main is the title, we plot the bars horizontally, the x-axis maximum is 100%,
#column names are in column 1 of imiss, the names should be in small font (cex), the orientation of lables (las), 
#and the space and width of the bars is reduced (space).
abline(v=70, col="red") #add a line at the cut-off
dev.off() #close the plot

#Inspect the plots by opening the PDF through Rstudio (File > Open file...) and return to the activity web page.
#We set a cut-off at 70% missing sites, meaning that we will remove the three individuals with more than 
#70% missing data.
#Return to the activity web page
