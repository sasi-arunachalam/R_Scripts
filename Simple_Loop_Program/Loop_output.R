#************************************************
#SASI ARUNACHALAM : St Jude Childern research hospital 
# Simple loop 
#************************************************
# This program is for creating input for matrix counts 
# choose only "Chr","Position","Position.1","Reference_Allele","Variant_Allele" .  
# The write out unique patient mutation by looping through

##############################

library(tidyverse) #ggplot2,tibble,tidyr,readr,purrr,dplyr,stringr,forcats
library(data.table)
library(readxl)
library(openxlsx)

rm(list=ls())
##############################
##read data 
pre <- read_excel("Loop_output_input.xlsx",sheet="Two_Platform")
pre <- pre %>% filter(Data=="Both")
pre<-pre%>%select("Sample", "Chr" ,"WU_HG19_Pos", "ReferenceAllele", "MutantAllele") #select only columns you need
colnames(pre) <- c("Sample","Chr","Position","Reference_Allele","Variant_Allele") # remane the columns

##############################
## loop
samples <- unique(pre$Sample)
for(i in 1:length(samples)){
  pre_out <- pre %>%filter(Sample== samples[i])
  pre_out <- pre_out%>%select("Chr","Position","Reference_Allele","Variant_Allele"  )
   write.table(pre_out, paste( "Loop_",samples[i], "_output.txt", sep=""), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}

########END########