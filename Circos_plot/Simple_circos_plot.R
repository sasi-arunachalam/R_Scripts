#************************************************
#SASI ARUNACHALAM PhD : St Jude Childern research hospital 
# Circos plot for fusion genes /Structural variants
#************************************************
library(readr)
library(data.table)
library(readxl)
library("circlize")
rm(list=ls())
#************************************************
##Read data
fusions<-read.table(file="fusions_reordered.txt",header=T) # Same list reordered
SVbed1Fusion <- read_excel("SVbed1Fusion1.xlsx")
SVbed2Fusion <-read_excel("SVbed1Fusion2.xlsx")

# Circos-Plot
pdf(file="Circos_plot.pdf")
circos.initializeWithIdeogram(species = "hg19",plotType = c("ideogram","labels"),ideogram.height=0.075)
circos.genomicLabels(fusions, labels.column = 4, side = "outside",cex=0.6)
circos.genomicLink(SVbed1Fusion,SVbed2Fusion,col="darkgrey")
dev.off()


