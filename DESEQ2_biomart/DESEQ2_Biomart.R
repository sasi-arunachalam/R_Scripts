#-----------------------------------------------------
#SASI ARUNACHALAM PhD : St Jude Childern Research Hospital
# DESEQ2 and Biomart
#-----------------------------------------------------
rm(list=ls())
#************************************************
#install.packages("htmltools") ## typed n for selection 
#library(htmltools)
#source("https://bioconductor.org/biocLite.R") 
#biocLite("DESeq2") # typed n for selection
#install.packages("memoise") # typed a
#library( "memoise" )
library( "DESeq2" )
library(ggplot2)
library( "biomaRt" )

#************************************************
### METHOD : 1
#************************************************
#DESEQ2: given data: ensembleID

countData <- read.csv('airway_rawcounts.csv', header = TRUE, sep = ",")
metaData <- read.csv('airway_metadata.csv', header = TRUE, sep = ",")

#Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~dex, tidy = TRUE)

#Run DESEQ function
dds <- DESeq(dds) 
#Results
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res)
res <- res[order(res$padj),]
head(res)

#**********
##Plots
#**********
#Countplots
par(mfrow=c(2,2))

plotCounts(dds, gene="ENSG00000179094", intgroup="dex")
plotCounts(dds, gene="ENSG00000116584", intgroup="dex")
plotCounts(dds, gene="ENSG00000189221", intgroup="dex")
plotCounts(dds, gene="ENSG00000120129", intgroup="dex")

## Volcanoplots
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


## PCA
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="dex")


#************************************************
### METHOD: 2
#************************************************
#DESEQ2: covert ensemble id to genenames 

#For plotting with Gene names rather than ensemble 
#id replace first column (ensemble id ) with gene name 
#ofcourse you will lose dimension(len) in the dataset as some ensemble id has no gene name and some have duplciates
#remove NA, duplicates and blanks in gene name colums(first column)

##Step 1
countData <- read.csv('airway_rawcounts.csv', header = TRUE, sep = ",")
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
#listAttributes(mart = ensembl)
genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = countData$ensgene,
                  mart = ensembl )
idx <- match( countData$ensgene, genemap$ensembl_gene_id )
countData$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
head(countData)

write.csv(as.data.frame(countData), file="ensemble_id_vs_gene.csv" ) #use this file for future referncefor any ensemble id and gene name

# Go to this  output ,"ensemble_id_vs_gene.csv", remove duplicates , NA and replace gene name in first column
# i renamed them to airway_rawcounts_geneadded.csv

##Step 2
#DESEQ2 starts here for METHOD 2

countData <- read.csv('airway_rawcounts_geneadded.csv', header = TRUE, sep = ",")
which(is.na(countData)) #this will tell us location of missing rows  (NA)
countData <- na.omit(countData) # omit the NA
metaData <- read.csv('airway_metadata.csv', header = TRUE, sep = ",")


#Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~dex, tidy = TRUE)

#Run DESEQ function
dds <- DESeq(dds) 
#Results
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res)
res <- res[order(res$padj),]
head(res)


#**********
##Plots
#**********
#Countplots
par(mfrow=c(2,2))

plotCounts(dds, gene="STOM", intgroup="dex")
plotCounts(dds, gene="ABCC3", intgroup="dex")
plotCounts(dds, gene="ABHD15", intgroup="dex")
plotCounts(dds, gene="ABHD10", intgroup="dex")



## Volcanoplots
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

## PCA
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="dex")




#************************************************
### Normalized counst for plotting
#************************************************
#"""DESeq2 doesn’t actually use normalized counts, rather it uses the raw counts and 
#models the normalization inside the Generalized Linear Model (GLM). 
#These normalized counts will be useful for downstream visualization of results, such as heatmap,clustering 
#but cannot be used as input to DESeq2 or any other tools that peform differential 
#expression analysis which use the negative binomial model."""
counts <- read.csv('airway_rawcounts.csv', header = TRUE, sep = ",")
geoMeans <- exp(rowMeans(log(counts(dds))))
dds <- estimateSizeFactors(dds,geoMeans=geoMeans)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
#head(counts)
#head(normalized_counts)
write.csv(as.data.frame(normalized_counts), file="normalized_counts.csv" ) 








