#************************************************
#SASI ARUNACHALAM PhD : St Jude Childern research hospital 
#************************************************
library(data.table)
rm(list=ls())

##Read all inputs

inMatrix <- as.data.frame(read.table("Sample1_Copy_input.txt", sep="\t",
                                     header=TRUE, row.names=NULL, stringsAsFactors=FALSE))
chromMatrix <- as.data.frame(read.table("ChromMax.txt", sep="\t",
                                        header=FALSE, row.names=NULL, stringsAsFactors=FALSE))
cutMatrix1 <-as.data.frame(read.table("Sample1_BAF_input.txt", sep="\t",
                                     header=TRUE, row.names=NULL, stringsAsFactors=FALSE))


chromList <- seq(1,22)
head(chromMatrix)
#Get chromosome position on the graph
chromPos <- NULL
chromPos <- c(chromPos, 0)
namePos <- NULL
namePos <- c(namePos, chromMatrix[2, 2] / 2)


runningTotal <- 0
noLohNum <- 0

for (rowIndex in 2:nrow(chromMatrix))
{
  prevChromLength <- chromMatrix[rowIndex - 1, 2]
  thisChromLength <- chromMatrix[rowIndex, 2]
  
  runningTotal <- runningTotal + prevChromLength
  namePosThisOne <- runningTotal + (thisChromLength / 2)
  
  chromPos <- c(chromPos, runningTotal)
  namePos <- c(namePos, namePosThisOne)
}

chromMatrix <- cbind(chromMatrix, chromPos, namePos)
head(chromMatrix)
colnames(chromMatrix) <- c("Chrom", "Length", "StartPos", "NamePos")



#Get absolute copy
absCopy <- 2^(inMatrix[, "log2"] + 1)


cutMatrix <- cbind(inMatrix, absCopy)


# Begin PNG :plot
png(paste("Sample1_CNV.png", sep="_"), height=500, width=1000)
par(mfrow=c(2,1), oma=c(1,1,1,1), mar=c(3.8,3.8,1,1))

totalLength <- sum(chromMatrix[which(chromMatrix[, "Chrom"] %in% chromList), "Length"])



#For Copy number
# create empty blank plot and later add on the CNV data
plot(c(1,2,3), col="white", xlim=c(0, totalLength), xlab="",
     ylim=c(0,6), xaxs="i", yaxs="i", ylab="abs. copy", xaxt="n", font.main=1, main="Sample1")

abline(v=chromMatrix[,"StartPos"])
abline(h=c(1,2,3,4,5), lty=2, col="gray77")

for (chrom in 1:22)
{
  thisNamePos <- chromMatrix[which(chromMatrix[,"Chrom"] == chrom), "NamePos"]
  
  text(x=thisNamePos, y=5.5, label=chrom)
}
colnames(cutMatrix)

#Add copy data to plot
for (rowIndex in 1:nrow(cutMatrix))
{
  chrom <- cutMatrix[rowIndex,"chromosome"]
  startPos <- cutMatrix[rowIndex,"start"]
  endPos <- cutMatrix[rowIndex,"end"]
  absCopy <- cutMatrix[rowIndex,"absCopy"]
  chromBegin <- chromMatrix[which(chromMatrix[,"Chrom"] == chrom), "StartPos"]
  startPos <- startPos + chromBegin
  endPos <- endPos + chromBegin
  segments(x0=startPos, y0=absCopy, x1=endPos, y1=absCopy, lwd=3, col="dodgerblue3")
 
}
###


totalLength <- sum(chromMatrix[which(chromMatrix[, "Chrom"] %in% chromList), "Length"])


# AI:BAF number
plot(c(1,2,3), col="white", xlim=c(0, totalLength), xlab="",
     ylim=c(0,1), xaxs="i", yaxs="i", ylab="AI", xaxt="n", font.main=1)#, main=thisPatient)

abline(v=chromMatrix[,"StartPos"])
abline(h=c(1,2,3,4,5), lty=2, col="gray77")

for (chrom in 1:22)
{
  thisNamePos <- chromMatrix[which(chromMatrix[,"Chrom"] == chrom), "NamePos"]
  
  text(x=thisNamePos, y=5.5, label=chrom)
}
colnames(cutMatrix1)
# Add BAF data to plot
for (rowIndex in 1:nrow(cutMatrix1))
{
  chrom <- cutMatrix1[rowIndex,"chromosome"]
  startPos <- cutMatrix1[rowIndex,"start"]
  endPos <- cutMatrix1[rowIndex,"end"]
  #absCopy <- cutMatrix[rowIndex,"absCopy"]
  absCopy <- cutMatrix1[rowIndex,"TumorVaf"]
  
  chromBegin <- chromMatrix[which(chromMatrix[,"Chrom"] == chrom), "StartPos"]
  
  startPos <- startPos + chromBegin
  endPos <- endPos + chromBegin
  
  segments(x0=startPos, y0=absCopy, x1=endPos, y1=absCopy, lwd=3, col="dodgerblue3")
  #segments(x0=startPos, y0=TumorVaf, x1=endPos, y1=TumorVaf, lwd=3, col="dodgerblue3")
}

dev.off()


#************************************************
##Germline and tumor BAF 

# Begin PNG :plot
png(paste("Sample1_BAF_germline_tumor_CNV.png", sep="_"), height=500, width=1000)
par(mfrow=c(2,1), oma=c(1,1,1,1), mar=c(3.8,3.8,1,1))

totalLength <- sum(chromMatrix[which(chromMatrix[, "Chrom"] %in% chromList), "Length"])

plot(c(1,2,3), col="white", xlim=c(0, totalLength), xlab="",
     ylim=c(0,1), xaxs="i", yaxs="i", ylab="AI", xaxt="n", font.main=1)#, main=thisPatient)

abline(v=chromMatrix[,"StartPos"])
abline(h=c(1,2,3,4,5), lty=2, col="gray77")

for (chrom in 1:22)
{
  thisNamePos <- chromMatrix[which(chromMatrix[,"Chrom"] == chrom), "NamePos"]
  
  text(x=thisNamePos, y=5.5, label=chrom)
}
colnames(cutMatrix1)
# Add Tumor BAF data to plot
for (rowIndex in 1:nrow(cutMatrix1))
{
  chrom <- cutMatrix1[rowIndex,"chromosome"]
  startPos <- cutMatrix1[rowIndex,"start"]
  endPos <- cutMatrix1[rowIndex,"end"]
  #absCopy <- cutMatrix[rowIndex,"absCopy"]
  absCopy <- cutMatrix1[rowIndex,"TumorVaf"]
  
  chromBegin <- chromMatrix[which(chromMatrix[,"Chrom"] == chrom), "StartPos"]
  
  startPos <- startPos + chromBegin
  endPos <- endPos + chromBegin
  
  segments(x0=startPos, y0=absCopy, x1=endPos, y1=absCopy, lwd=3, col="dodgerblue3")
 
}

# Add Germline BAF data to plot
for (rowIndex in 1:nrow(cutMatrix1))
{
  chrom <- cutMatrix1[rowIndex,"chromosome"]
  startPos <- cutMatrix1[rowIndex,"start"]
  endPos <- cutMatrix1[rowIndex,"end"]
  #absCopy <- cutMatrix[rowIndex,"absCopy"]
  absCopy <- cutMatrix1[rowIndex,"GermlineVaf"]
  
  chromBegin <- chromMatrix[which(chromMatrix[,"Chrom"] == chrom), "StartPos"]
  
  startPos <- startPos + chromBegin
  endPos <- endPos + chromBegin
  
  segments(x0=startPos, y0=absCopy, x1=endPos, y1=absCopy, lwd=3, col="red")
  
}

dev.off()



