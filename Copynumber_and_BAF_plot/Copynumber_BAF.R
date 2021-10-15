#************************************************
#SASI ARUNACHALAM PhD : St Jude Childern research hospital 
#************************************************

setwd("~/Desktop/311_cnv/snp")
thisPatient <- "SJHGG016311_A1"
Patient <-paste(thisPatient,"txt", sep=".")

inMatrix <- as.data.frame(read.table(Patient, sep="\t",
                                     header=TRUE, row.names=NULL, stringsAsFactors=FALSE))
chromMatrix <- as.data.frame(read.table("ChromMax.txt", sep="\t",
                                        header=FALSE, row.names=NULL, stringsAsFactors=FALSE))

setwd("~/Desktop/325cnv/bplot")

cutMatrix1 <-as.data.frame(read.table("SJHGG016311_A1_G1_bplotinput.txt", sep="\t",
                                     header=TRUE, row.names=NULL, stringsAsFactors=FALSE))


chromList <- seq(1,22)
head(chromMatrix)
# get chromosome position on the graph
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



# get absolute copy
absCopy <- 2^(inMatrix[, "log2"] + 1)


cutMatrix <- cbind(inMatrix, absCopy)

colnames(cutMatrix)
#write.table(cutMatrix, file = "A1.txt", sep = "\t",
           # row.names = TRUE, col.names = NA)

#cutMatrix <- inMatrix

#paste("Hello", "world", sep=" ")

# begin PNG
png(paste(thisPatient,"CNV.png", sep="_"), height=500, width=1000)
par(mfrow=c(2,1), oma=c(1,1,1,1), mar=c(3.8,3.8,1,1))

totalLength <- sum(chromMatrix[which(chromMatrix[, "Chrom"] %in% chromList), "Length"])



#for copy number
# create empty blank plot and later add on the CNV data
plot(c(1,2,3), col="white", xlim=c(0, totalLength), xlab="",
     ylim=c(0,6), xaxs="i", yaxs="i", ylab="abs. copy", xaxt="n", font.main=1, main=thisPatient)

#for AI number
#plot(c(1,2,3), col="white", xlim=c(0, totalLength), xlab="",
    # ylim=c(0,1), xaxs="i", yaxs="i", ylab="AI", xaxt="n", font.main=1, main=thisPatient)

#axis(1, labels=c("Chr"), at=1400000000, col.ticks="white")
#axis(1,  at=1400000000, col.ticks="white")
abline(v=chromMatrix[,"StartPos"])
abline(h=c(1,2,3,4,5), lty=2, col="gray77")

for (chrom in 1:22)
{
  thisNamePos <- chromMatrix[which(chromMatrix[,"Chrom"] == chrom), "NamePos"]
  
  text(x=thisNamePos, y=5.5, label=chrom)
}
colnames(cutMatrix)
# add copy data to plot
for (rowIndex in 1:nrow(cutMatrix))
{
  chrom <- cutMatrix[rowIndex,"chromosome"]
  startPos <- cutMatrix[rowIndex,"start"]
  endPos <- cutMatrix[rowIndex,"end"]
  absCopy <- cutMatrix[rowIndex,"absCopy"]
  #absCopy <- cutMatrix[rowIndex,"TumorVaf"]
  
  chromBegin <- chromMatrix[which(chromMatrix[,"Chrom"] == chrom), "StartPos"]
  
  startPos <- startPos + chromBegin
  endPos <- endPos + chromBegin
  
  segments(x0=startPos, y0=absCopy, x1=endPos, y1=absCopy, lwd=3, col="dodgerblue3")
  #segments(x0=startPos, y0=TumorVaf, x1=endPos, y1=TumorVaf, lwd=3, col="dodgerblue3")
}
###


totalLength <- sum(chromMatrix[which(chromMatrix[, "Chrom"] %in% chromList), "Length"])


#for copy number
# create empty blank plot and later add on the CNV data
#plot(c(1,2,3), col="white", xlim=c(0, totalLength), xlab="",
    # ylim=c(0,6), xaxs="i", yaxs="i", ylab="abs. copy", xaxt="n", font.main=1, main=thisPatient)



#for AI number
plot(c(1,2,3), col="white", xlim=c(0, totalLength), xlab="",
     ylim=c(0,1), xaxs="i", yaxs="i", ylab="AI", xaxt="n", font.main=1)#, main=thisPatient)

#axis(1, labels=c("Chr"), at=1400000000, col.ticks="white")
#axis(1,  at=1400000000, col.ticks="white")
abline(v=chromMatrix[,"StartPos"])
abline(h=c(1,2,3,4,5), lty=2, col="gray77")

for (chrom in 1:22)
{
  thisNamePos <- chromMatrix[which(chromMatrix[,"Chrom"] == chrom), "NamePos"]
  
  text(x=thisNamePos, y=5.5, label=chrom)
}
colnames(cutMatrix1)
# add copy data to plot
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

