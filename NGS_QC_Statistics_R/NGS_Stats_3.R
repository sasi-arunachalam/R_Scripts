#-----------------------------------------------------
#SASI ARUNACHALAM PhD : St Jude Childern Research Hospital
# NGS_Stats
#-----------------------------------------------------
library(tidyverse) #ggplot2,tibble,tidyr,readr,purrr,dplyr,stringr,forcats
library(data.table)
library(readxl)
library(openxlsx)
library(scales)
rm(list=ls())
##############



stat <- read_excel("NGS_Stats_Input_3.xlsx")
                  

##ggplot for Total percent read map 
ggplot(data=stat, aes(x=fct_reorder(Sample,Type),y=Percent_readmapped,fill=factor(Type))) +
  geom_bar(stat="identity") + 
  ggtitle("Reads mapped %")+
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

ggsave("~working_directory/NGS_Stats_output_3/readmapped.pdf",width = 50, height = 10, units = "cm")
  
##ggplot for Totalreads read map 
ggplot(data=stat, aes(x=fct_reorder(Sample,Type),y=Total_reads,fill=factor(Type))) +
  geom_bar(stat="identity") + 
  ggtitle("Total reads (million) ")+
  theme_classic()+
  scale_y_continuous(labels=function(n){format(n, scientific = FALSE)})+
geom_hline(yintercept=100000000, linetype="dashed")+
  theme(
    legend.title=element_blank(),  
    #legend.position=c(.73,.7),
    axis.title.y=element_blank(), 
    text=element_text(family="serif",size=20),
    plot.title=element_text(face="bold",hjust=c(0,0))
  )
ggsave("~working_directory/NGS_Stats_output_3/Totalreads.pdf",width = 50, height = 10, units = "cm")

##ggplot for duplicate rate 

ggplot(data=stat, aes(x=fct_reorder(Sample,Type),y=Duplicaterate,fill=factor(Type))) +
  geom_bar(stat="identity") + 
  ggtitle("Duplicate reads % ")+
  theme_classic()+
   geom_hline(yintercept=40, linetype="dashed")+
  theme(
    legend.title=element_blank(),  
    #legend.position=c(.73,.7),
    axis.title.y=element_blank(), 
    text=element_text(family="serif",size=20),
    plot.title=element_text(face="bold",hjust=c(0,0))
  )
ggsave("~working_directory/NGS_Stats_output_3/Duplicate.pdf",width = 50, height = 10, units = "cm")



##ggplot for bamsize
ggplot(data=stat, aes(x=fct_reorder(Sample,Type),y=BamSize_G,fill=factor(Type))) +
  geom_bar(stat="identity") + 
  ggtitle("Data output (Gb) ")+
  theme_classic()+
  geom_hline(yintercept=10, linetype="dashed")+
  theme(
    legend.title=element_blank(),  
    #legend.position=c(.73,.7),
    axis.title.y=element_blank(), 
    text=element_text(family="serif",size=20),
    plot.title=element_text(face="bold",hjust=c(0,0))
  )
ggsave("~working_directory/NGS_Stats_output_3/Dataoutput.pdf",width = 50, height = 10, units = "cm")



##ggplot for coverage
ggplot(data=stat, aes(x=fct_reorder(Sample,Type),y= coverage_20X,fill=factor(Type))) +
  geom_bar(stat="identity") + 
  ggtitle("Coverage min. 20X (%) ")+
  theme_classic()+
  geom_hline(yintercept=40, linetype="dashed")+
  theme(
    legend.title=element_blank(),  
    #legend.position=c(.73,.7),
    axis.title.y=element_blank(), 
    text=element_text(family="serif",size=20),
    plot.title=element_text(face="bold",hjust=c(0,0))
  )
ggsave("~working_directory/NGS_Stats_output_3/Coverage.pdf",width = 50, height = 10, units = "cm")



##ggplot for coverage
ggplot(data=stat, aes(x=fct_reorder(Sample,coverage_20X),y= coverage_20X,fill=factor(Type))) +
  geom_bar(stat="identity") + 
  ggtitle("Coverage min. 20X (%) ")+
  theme_classic()+
  geom_hline(yintercept=40, linetype="dashed")+
  theme(
    legend.title=element_blank(),  
    #legend.position=c(.73,.7),
    axis.title.y=element_blank(), 
    text=element_text(family="serif",size=20),
    plot.title=element_text(face="bold",hjust=c(0,0))
  )
ggsave("~working_directory/NGS_Stats_output_3/Coverage_sorted.pdf",width = 50, height = 10, units = "cm")







