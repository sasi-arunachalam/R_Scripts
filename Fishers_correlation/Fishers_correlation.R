#-----------------------------------------------------
#SASI ARUNACHALAM PhD : St Jude Childern Research Hospital
# Fishers_Correlation_calculation
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
##Fishers
###################

# read the data
Supplementary_new_final <- read_excel("Fishers_input.xlsx", sheet = "Tier3")

#Calculate VAF
Supplementary_new_final <- Supplementary_new_final %>% mutate (wgs_p_vaf = Mutant_In_Tumor/ Total_In_Tumor )
Supplementary_new_final <- Supplementary_new_final %>% mutate (wgs_x_vaf = Mutant_In_Xeno/ Total_In_Xeno )
Supplementary_new_final <- Supplementary_new_final %>% mutate (Val_p_vaf = Mutant_In_Tumor_Validation/ Total_In_Tumor_Validation )
Supplementary_new_final <- Supplementary_new_final %>% mutate (Val_x_vaf = Mutant_In_Xeno_Validation/ Total_In_Xeno_Validation )


###################
# Plot the correlation 
###################



# Plot the corelation for validated ( non-resue)
WGS <- Supplementary_new_final %>% filter( Validation == "shared" | Validation == "Xeno-only" | Validation == "patient-only")
pdf("Correlation_exploratory_validated.pdf")
ggscatter(WGS, x = "wgs_p_vaf", y = "wgs_x_vaf", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "WGS_P", ylab = "WGS_XenoCP")
ggscatter(WGS, x = "Val_p_vaf", y = "Val_x_vaf", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Val_P", ylab = "Val_XenoCP")
dev.off()
# Plot the corelation for rescue
Rescue <- Supplementary_new_final %>% filter( Validation == "not_in_validation")
pdf("Correlation_exploratory_rescue.pdf")
ggscatter(Rescue, x = "wgs_p_vaf", y = "wgs_x_vaf", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Wgs_P_rescue", ylab = "Wgs_XenoCP_rescue")
dev.off()


###################
##Fishers correlation calculation
###################
options(scipen=999)
dt <- data.table(WGS) # for validated , non-resue


dt[, p.val_Tumor := fisher.test(matrix(c(Total_In_Tumor, Total_In_Normal, Mutant_In_Tumor, Mutant_In_Normal), ncol=2), workspace=1e4)$p.value, by=Chr.Pos.Ref.Alt]
dt[, p.val_Xeno := fisher.test(matrix(c(Total_In_Xeno, Total_In_Normal, Mutant_In_Xeno, Mutant_In_Normal), ncol=2), workspace=1e4)$p.value, by=Chr.Pos.Ref.Alt]
dt[, p.val_Val_tumor := fisher.test(matrix(c(Total_In_Tumor_Validation, Total_In_Normal_Validation, Mutant_In_Tumor_Validation, Mutant_In_Normal_Validation), ncol=2), workspace=1e4)$p.value, by=Chr.Pos.Ref.Alt]
dt[, p.val_Xeno_tumor := fisher.test(matrix(c(Total_In_Xeno_Validation, Total_In_Normal_Validation, Mutant_In_Xeno_Validation, Mutant_In_Normal_Validation), ncol=2), workspace=1e4)$p.value, by=Chr.Pos.Ref.Alt]
df <- as.data.frame(dt)


dt1 <- data.table(Rescue)# for  non-resue
dt1[, p.val_Tumor := fisher.test(matrix(c(Total_In_Tumor, Total_In_Normal, Mutant_In_Tumor, Mutant_In_Normal), ncol=2), workspace=1e4)$p.value, by=Chr.Pos.Ref.Alt]
dt1[, p.val_Xeno := fisher.test(matrix(c(Total_In_Xeno, Total_In_Normal, Mutant_In_Xeno, Mutant_In_Normal), ncol=2), workspace=1e4)$p.value, by=Chr.Pos.Ref.Alt]
df1 <- as.data.frame(dt1)
df1 <- df1 %>% mutate (p.val_Val_tumor= "NA") %>%mutate (p.val_Xeno_tumor= "NA")  # this column is set to NA as no data to coorelated, so its easy for r bind

#r bind the above two data frame
final <- rbind(df,df1)

#Write out the file
write.xlsx(final, file = "Fishers_output.xlsx")

### END####

