#-----------------------------------------------------
#SASI ARUNACHALAM PhD : St Jude Childern Research Hospital
#NGS_stas
#-----------------------------------------------------
library(tidyverse) #ggplot2,tibble,tidyr,readr,purrr,dplyr,stringr,forcats
library(data.table)
library(readxl)
library(openxlsx)
library(scales)
rm(list=ls())
#############


stat <- read_excel("NGS_Stats_Input_2.xlsx")

###### Qc+Fusion
##1


##ggplot for Duplicate reads
ggplot(data=stat, aes(x=fct_reorder(Sample,Fusion),y=Dup_per,fill=factor(Fusion))) +
  geom_bar(stat="identity") + 
  ggtitle("Dup_percentage")+
  theme_classic()+
  geom_hline(yintercept=40, linetype="dashed")+
  theme(
    legend.title=element_blank(),  
    #legend.position=c(.73,.7),
    axis.title.y=element_blank(), 
    text=element_text(family="serif",size=20),
    plot.title=element_text(face="bold",hjust=c(0,0))
  )
ggsave("~/Desktop/R_Scripts_compilation/NGS_QC_Stats_R/NGS_Stats_output_2/Dup+Fusion.pdf",width = 50, height = 10, units = "cm")


##2
##ggplot for Total reads
ggplot(data=stat, aes(x=fct_reorder(Sample,Fusion),y=Total_reads,fill=factor(Fusion))) +
  geom_bar(stat="identity") + 
  ggtitle("Total_reads")+
  theme_classic()+
  #scale_y_continuous(labels=function(n){format(n, scientific = FALSE)})+
  geom_hline(yintercept=100000000, linetype="dashed")+
  theme(
    legend.title=element_blank(),  
    #legend.position=c(.73,.7),
    axis.title.y=element_blank(), 
    text=element_text(family="serif",size=20),
    plot.title=element_text(face="bold",hjust=c(0,0))
  )
ggsave("~/Desktop/R_Scripts_compilation/NGS_QC_Stats_R/NGS_Stats_output_2/Totalreads+Fusion.pdf",width = 50, height = 10, units = "cm")

##3
##ggplot for Mapped reads
ggplot(data=stat, aes(x=fct_reorder(Sample,Fusion),y=Mapped_reads,fill=factor(Fusion))) +
  geom_bar(stat="identity") + 
  ggtitle("Mapped_reads")+
  theme_classic()+
  #scale_y_continuous(labels=function(n){format(n, scientific = FALSE)})+
  geom_hline(yintercept=100000000, linetype="dashed")+
  theme(
    legend.title=element_blank(),  
    #legend.position=c(.73,.7),
    axis.title.y=element_blank(), 
    text=element_text(family="serif",size=20),
    plot.title=element_text(face="bold",hjust=c(0,0))
  )
ggsave("~/Desktop/R_Scripts_compilation/NGS_QC_Stats_R/NGS_Stats_output_2/Mapped_reads+Fusion.pdf",width = 50, height = 10, units = "cm")


##4
##ggplot for NonDupMpd_reads
ggplot(data=stat, aes(x=fct_reorder(Sample,Fusion),y=NonDupMpd_reads,fill=factor(Fusion))) +
  geom_bar(stat="identity") + 
  ggtitle("NonDupMpd_reads")+
  theme_classic()+
  #scale_y_continuous(labels=function(n){format(n, scientific = FALSE)})+
  geom_hline(yintercept=100000000, linetype="dashed")+
  theme(
    legend.title=element_blank(),  
    #legend.position=c(.73,.7),
    axis.title.y=element_blank(), 
    text=element_text(family="serif",size=20),
    plot.title=element_text(face="bold",hjust=c(0,0))
  )
ggsave("~/Desktop/R_Scripts_compilation/NGS_QC_Stats_R/NGS_Stats_output_2/NonDupMpd_reads+Fusion.pdf",width = 50, height = 10, units = "cm")




##5
##ggplot for mapped_per
ggplot(data=stat, aes(x=fct_reorder(Sample,Fusion),y=Mapped_per,fill=factor(Fusion))) +
  geom_bar(stat="identity") + 
  ggtitle("Mapped_percentage")+
  theme_classic()+
  coord_cartesian(ylim=c(80, 100))+
  geom_hline(yintercept=96, linetype="dashed")+
  theme(
    legend.title=element_blank(),  
    #legend.position=c(.73,.7),
    axis.title.y=element_blank(), 
    text=element_text(family="serif",size=20),
    plot.title=element_text(face="bold",hjust=c(0,0))
  )
ggsave("~/Desktop/R_Scripts_compilation/NGS_QC_Stats_R/NGS_Stats_output_2/Mapped_per+Fusion.pdf",width = 50, height = 10, units = "cm")




##6
##ggplot for coverage reads
ggplot(data=stat, aes(x=fct_reorder(Sample,Fusion),y=Coverage_20X,fill=factor(Fusion))) +
  geom_bar(stat="identity") + 
  ggtitle("Coverage min. 20X(%) ")+
  theme_classic()+
  geom_hline(yintercept=30, linetype="dashed")+
  theme(
    legend.title=element_blank(),  
    #legend.position=c(.73,.7),
    axis.title.y=element_blank(), 
    text=element_text(family="serif",size=20),
    plot.title=element_text(face="bold",hjust=c(0,0))
  )
ggsave("~/Desktop/R_Scripts_compilation/NGS_QC_Stats_R/NGS_Stats_output_2/Coverage+Fusion.pdf",width = 50, height = 10, units = "cm")



##7
##ggplot for bamsize
ggplot(data=stat, aes(x=fct_reorder(Sample,Fusion),y=Bamsize_Gb,fill=factor(Fusion))) +
  geom_bar(stat="identity") + 
  ggtitle("Bamsize_Gb")+
  theme_classic()+
  #geom_hline(yintercept=30, linetype="dashed")+
  theme(
    legend.title=element_blank(),  
    #legend.position=c(.73,.7),
    axis.title.y=element_blank(), 
    text=element_text(family="serif",size=20),
    plot.title=element_text(face="bold",hjust=c(0,0))
  )
ggsave("~/Desktop/R_Scripts_compilation/NGS_QC_Stats_R/NGS_Stats_output_2/Bamsize+Fusion.pdf",width = 50, height = 10, units = "cm")


##8
ggplot(data=stat, aes(x=Sample)) +
  geom_bar(aes(y=Total_reads), stat="identity", position ="identity", alpha=.3, fill='lightblue', color='lightblue4') +
  geom_bar(aes(y=Mapped_reads), stat="identity", position="identity", alpha=.8, fill='pink', color='red')

##8
ggplot(data=stat, aes(x=fct_reorder(Sample,Fusion),fill=factor(Fusion))) +
  geom_bar(aes(y=Total_reads), stat="identity", position ="identity", alpha=.3, fill='blue', color='lightblue4') +
  geom_bar(aes(y=Mapped_reads), stat="identity", position="identity", alpha=.8, fill='lavender', color='red')+
  geom_bar(aes(y=NonDupMpd_reads), stat="identity", position="identity", alpha=.8, fill='pink', color='red')+
  theme_classic()+
  #scale_y_continuous(labels=function(n){format(n, scientific = FALSE)})+
  geom_hline(yintercept=100000000, linetype="dashed")+
  theme(
    legend.title=element_blank(),  
    #legend.position=c(.73,.7),
    axis.title.y=element_blank(), 
    text=element_text(family="serif",size=20),
    plot.title=element_text(face="bold",hjust=c(0,0))
  )
ggsave("~/Desktop/R_Scripts_compilation/NGS_QC_Stats_R/NGS_Stats_output_2/Allreads+Fusion.pdf",width = 50, height = 10, units = "cm")





