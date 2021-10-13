#************************************************
#SASI ARUNACHALAM PhD: St Jude Childern research hospital 
# This program calculate PCA, Plots PCA and writes out, Top 100 features
#************************************************
#install.packages("devtools")
#devtools::install_github("HuntsmanCancerInstitute/hciR")
library(hciR)
library(readr)
library(SummarizedExperiment)
library("DESeq2")
library("Biobase")
library(ggrepel)
library(ggplot2)
library(data.table)
rm(list=ls())
#************************************************
#Read data
myfpkm <- readr::read_tsv(file = "PCA_Top_100_features_input.txt")
sidx <- !colnames(myfpkm) %in% c("gene")
mat_fpkm <- as.matrix(myfpkm[,sidx])  #except gene column
rownames(mat_fpkm) <- myfpkm$gene # gene column as row names

#Plot PCA
pdf("PCA.pdf")
p <- plot_pca(log2(mat_fpkm + 1), ggplot=TRUE)
q<- p + geom_text_repel(aes(label=COLNAMES), cex=3, box.padding=.1)
q
dev.off()

#Write out PC1 and PC2 
p$data  ## gives PC1 and PC2 values
fwrite(p$data,"PC1vsPC2.csv",row.names = TRUE)

#Top 100 features
x<-mat_fpkm
n <- apply(x, 1, var) # apply var function
top100 <- head(x[order(n, decreasing = TRUE),  ], 100) ## choose top 100 features 
top100<-as.data.frame(top100)

#Write out
fwrite(top100,"Top100_Features.csv",row.names = TRUE)

####END####

