#************************************************
#SASI ARUNACHALAM : St Jude Childern research hospital 
# Fish plot
#************************************************
#library(devtools)
#install_github("chrisamiller/fishplot")
#************************************************
library(fishplot)
rm(list=ls())
#************************************************
#### Patient Sample1
#************************************************
fishTable <- read.table("Sample1_FishPlot_input.txt", header=FALSE, row.names=1, stringsAsFactors=FALSE, sep="\t")
timepoints <- as.vector(unlist(fishTable[nrow(fishTable),]))
cloneTable <- 200*as.matrix(fishTable[1:(nrow(fishTable)-1), ]) ## convert vaf to ccf percentage
parents <- c(0, 1, 1, 2, 2)
pdf("Sample1_fishplot.pdf")
#create a fish object
fish = createFishObject(cloneTable, parents, timepoints=timepoints, fix.missing.clones=TRUE)

#calculate the layout of the drawing
fish = layoutClones(fish)

#draw the plot, using the splining method (recommended)
#and providing both timepoints to label and a plot title

fishPlot(fish,shape="spline",title.btm="Sample1",
         cex.title=0.5, vlines=c(0,2000), 
         vlab=c("day 0","day 2000"))
dev.off()

#************************************************
#### Patient Sample2
#************************************************
fishTable <- read.table("Sample2_FishPlot_input.txt", header=FALSE, row.names=1, stringsAsFactors=FALSE, sep="\t")
timepoints <- as.vector(unlist(fishTable[nrow(fishTable),]))
cloneTable <- 200*as.matrix(fishTable[1:(nrow(fishTable)-1), ]) ## convert vaf to ccf percentage
parents <- c(0, 1, 1, 3)

#create a fish object
pdf("Sample2_fishplot.pdf")
fish = createFishObject(cloneTable, parents, timepoints=timepoints) #, fix.missing.clones=TRUE)

#calculate the layout of the drawing
fish = layoutClones(fish)

#draw the plot, using the splining method (recommended)
#and providing both timepoints to label and a plot title

fishPlot(fish,shape="spline",title.btm="Sample1",
         cex.title=0.5, vlines=c(0,2000), 
         vlab=c("day 0","day 2000"))
dev.off()


