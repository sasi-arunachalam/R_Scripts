#************************************************
#SASI ARUNACHALAM PhD : St Jude Childern research hospital 
# t-SNE plot for expression data
#************************************************
# This program  craetes t-SNE plot for expression data 
# Two different disease comdition namely Dis1 and Dis2 .  For eg: Dis1 could AML , Dis2 could be melanoma 
##############################
library(readr)
library(Rtsne)
library(data.table)
library(tidyverse) #ggplot2,tibble,tidyr,readr,purrr,dplyr,stringr,forcats
library(data.table)
library(readxl)
library(openxlsx)
rm(list=ls())
#************************************************

# Read data
raw_data<-read.delim("t_SNE_input.txt")

#Data preprocessing
sidx <- !colnames(raw_data) %in% c("gene")
mat_fpkm <- as.matrix(raw_data[,sidx])
rownames(mat_fpkm) <- raw_data$gene
#Filter out low-expression features.
keep_genes <- apply(mat_fpkm, 1, max) > 3
mat_fpkm <- mat_fpkm[keep_genes,]
mat_fpkm1 <-log2(mat_fpkm + 1)

#t-SNE plot 


training <- mat_fpkm1
#training_set<- t(training)
training_set<- as.data.frame(training_set)
training_set <-training_set %>%rownames_to_column( var="gene") 
training_set$label <- ifelse(grepl("Dis1", training_set$gene),"Dis1","Dis2")
training_set$label <- as.factor(training_set$label)
train<-training_set
##Performing t-SNE
tsne <- Rtsne(train[,-1], dims = 2, perplexity=50, verbose=TRUE, max_iter = 500,check_duplicates = FALSE)
colors = rainbow(length(unique(train$label)))
names(colors) = unique(train$label)

##Exploratory plotting
#plot on R studio
par(mgp=c(2.5,1,0))
plot(tsne$Y, t='n', main="tSNE", xlab="tSNE dimension 1", ylab="tSNE dimension 2", "cex.main"=2, "cex.lab"=1.5)
text(tsne$Y, labels=train$label, col=colors[train$label])


#Iterate different perplexity_values and choose one
tsne_plot <- function(perpl=30,iterations=500,learning=200){
  set.seed(1) # for reproducibility
  tsne <- Rtsne(train[,-1], dims = 2, perplexity=perpl, verbose=TRUE, max_iter=iterations, eta=learning)
  plot(tsne$Y, t='n', main = print(paste0("perplexity = ",perpl, ", max_iter = ",iterations, ", learning rate = ",learning)), xlab="tSNE dimension 1", ylab="tSNE dimension 2", "cex.main"=1, "cex.lab"=1.5)
  text(tsne$Y, labels=train$label, col=colors[train$label])
}
perplexity_values <- c(2,5,30,50,100)
sapply(perplexity_values,function(i){tsne_plot(perpl=i)})

#************************************************
#Plot t-SNE with choosen perplexity_values
pdf("t-SNE_perplexity30.pdf")
perplexity_values <- c(30)

sapply(perplexity_values,function(i){tsne_plot(perpl=i)})
dev.off()


##END

