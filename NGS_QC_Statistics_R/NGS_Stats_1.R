#-----------------------------------------------------
#SASI ARUNACHALAM PhD : St Jude Childern Research Hospital
# NGS_Stats
#-----------------------------------------------------
#############
library(tidyverse) #ggplot2,tibble,tidyr,readr,purrr,dplyr,stringr,forcats
library(data.table)
library(readxl)
library(openxlsx)
library(scales)
rm(list=ls())
#############
stat <- read_excel("NGS_Stats_Input_1.xlsx")

###### ALL_reads

##1
##ggplot for Allreads

ggplot(data=stat, aes(x=fct_reorder(Sample,Fusion),fill=factor(Fusion))) +
  geom_bar(aes(y=Totalreads), stat="identity", position ="identity", alpha=.3, fill='blue', color='lightblue4') +
  geom_bar(aes(y=MappedReads), stat="identity", position="identity", alpha=.8, fill='lavender', color='red')+
  geom_bar(aes(y=Non_dupped), stat="identity", position="identity", alpha=.8, fill='pink', color='red')+
  ggtitle("All_Reads")+
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

  
ggsave("~/Desktop/R_Scripts_compilation/NGS_QC_Stats_R/NGS_Stats_output_1/Allreads.pdf",width = 50, height = 10, units = "cm")




##2
##ggplot for mapped_per
ggplot(data=stat, aes(x=Sample, y=Map_perc)) +
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
ggsave("~/Desktop/R_Scripts_compilation/NGS_QC_Stats_R/NGS_Stats_output_1/mapped_per.pdf",width = 50, height = 10, units = "cm")


##3
##ggplot for Dup_per
ggplot(data=stat, aes(x=Sample, y=Dup_perc)) +
  geom_bar(stat="identity") + 
  ggtitle("Dup_percentage")+
  theme_classic()+
  coord_cartesian(ylim=c(0, 50))+
  geom_hline(yintercept=40, linetype="dashed")+
  theme(
    legend.title=element_blank(),  
    #legend.position=c(.73,.7),
    axis.title.y=element_blank(), 
    text=element_text(family="serif",size=20),
    plot.title=element_text(face="bold",hjust=c(0,0))
    
  )
ggsave("~/Desktop/R_Scripts_compilation/NGS_QC_Stats_R/NGS_Stats_output_1/dup_prec.pdf",width = 50, height = 10, units = "cm")


##4
##ggplot for OnTargetDepth
ggplot(data=stat, aes(x=Sample, y=OnTargetDepth)) +
  geom_bar(stat="identity") + 
  ggtitle("Target_depth")+
  theme_classic()+
  coord_cartesian(ylim=c(0, 100))+
  geom_hline(yintercept=75, linetype="dashed")+
  theme(
    legend.title=element_blank(),  
    #legend.position=c(.73,.7),
    axis.title.y=element_blank(), 
    text=element_text(family="serif",size=20),
    plot.title=element_text(face="bold",hjust=c(0,0))
    
  )

ggsave("~/Desktop/R_Scripts_compilation/NGS_QC_Stats_R/NGS_Stats_output_1/Target.pdf",width = 50, height = 10, units = "cm")
