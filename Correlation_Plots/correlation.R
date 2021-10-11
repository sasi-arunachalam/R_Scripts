#-----------------------------------------------------
#SASI ARUNACHALAM PhD : St Jude Childern Research Hospital
# Correlation_stats
#-----------------------------------------------------

#install.packages("ggpubr")
library(tidyverse) #ggplot2,tibble,tidyr,readr,purrr,dplyr,stringr,forcats
library(data.table)
library(readxl)
library(openxlsx)
library(scales)
library("ggpubr")

rm(list=ls())
###################
#pre
###################
pre <- read_excel("Correlation_Input1.xlsx", 
                                   sheet = "For_python_correlation_pre")
pdf("Pre_plots.pdf")
ggscatter(pre, x = "Ori_vaf", y = "Redux_vaf", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Pipeline_VAF", ylab = "Redux_VAF")
ggscatter(pre, x = "MinD", y = "Mut_in_Tumor", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Pipeline_ALT", ylab = "Redux_ALT")
ggscatter(pre, x = "TinD", y = "Total_in_Tumor", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Pipeline_TotalTumor", ylab = "Redux_TotalTumor")
ggscatter(pre, x = "TinN", y = "Total_in_normal", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Pipeline_TotalNormal", ylab = "Redux_TotalNormal")
dev.off()

 

###################
#post
###################

post <- read_excel("Correlation_Input1", 
                  sheet = "For_python_correlation_post")
pdf("Post_plots.pdf")
ggscatter(post, x = "Ori_vaf", y = "Redux_vaf", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Pipeline_VAF", ylab = "Redux_VAF")
ggscatter(post, x = "MinR", y = "Mut_in_Tumor", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Pipeline_ALT", ylab = "Redux_ALT")
ggscatter(post, x = "TinR", y = "Total_in_Tumor", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Pipeline_TotalTumor", ylab = "Redux_TotalTumor")
ggscatter(post, x = "TinN", y = "Total_in_normal", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Pipeline_TotalNormal", ylab = "Redux_TotalNormal")
dev.off()
###################
#Counts
###################
patient_counts <- read_excel("Correlation_Input2.xlsx", 
                             sheet = "Sheet3")
pdf("Counts.pdf")

ggscatter(patient_counts, x = "hg18_p_vaf", y = "hg19_p_vaf", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Pipeline_P_HG18", ylab = "Pipeline_P_HG19")


ggscatter(patient_counts, x = "hg18_x_vaf", y = "hg19_x_vaf", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Pipeline_X_HG18", ylab = "PipelineX_HG19")



ggscatter(patient_counts, x = "matrix_P_vaf", y = "hg19_p_vaf", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Matrix_P_HG19", ylab = "Pipeline_P_HG19")


ggscatter(patient_counts, x = "matrix_X_vaf", y = "hg19_x_vaf", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Matrix_X_HG19", ylab = "Pipeline_X_HG19")



dev.off()







