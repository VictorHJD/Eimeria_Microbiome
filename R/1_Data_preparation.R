###Mice infection experiment summary
library("lifecycle", lib.loc="/usr/local/lib/R/site-library") 
library("ggplot2")
library("data.table")
library("tidyverse")
require("ggpubr")
library("dplyr")
library("plyr")
library("vegan")
library("gridExtra")
library("grid")
library("lattice")
library("pheatmap")
library("viridisLite")
#library("rcompanion")
#library("FSA")
library("phyloseq")

##Load data
sample.data <- read.csv("/SAN/Victors_playground/Eimeria_microbiome/sample_data_infb_Exp005.csv")
exp.des <- read.csv("/SAN/Victors_playground/Eimeria_microbiome/Inf1b_Exp005.DESIGN.csv")
exp.des$InfectionStrain <- NULL
exp.des$Departure <- NULL
exp.des$NOTE <- NULL
  
  ###Use function from Alice to calculate OPG
  calculateOPG <- function(sample.data){
    sample.data$mean_Neubauer <- 
      (sample.data$oocyst_sq1 + sample.data$oocyst_sq2 + sample.data$oocyst_sq3 + sample.data$oocyst_sq4) / 4
    # NB! Limit of detection = 1 oocysts
    sample.data$mean_Neubauer[sample.data$oocyst_sq1 + sample.data$oocyst_sq2 + sample.data$oocyst_sq3 + sample.data$oocyst_sq4 == 1] <- 0
    sample.data$oocysts.per.tube <- sample.data$mean_Neubauer * 10000 * sample.data$dilution
    sample.data$OPG <- sample.data$oocysts.per.tube / sample.data$fecweight_flot
    ## If we don't have the fecal weight BUT we counted in Neubauer chamber 0, then OPG = 0
    sample.data$oocysts.per.tube[sample.data$fecweight_flot == 0 & sample.data$mean_Neubauer == 0] <- 0
    sample.data$OPG[sample.data$fecweight_flot == 0 & sample.data$mean_Neubauer == 0] <- 0
    return(sample.data)
  }
  
sample.data<- calculateOPG(sample.data = sample.data)
  
##Check for spaces
sample.data$EH_ID <- gsub(pattern = " ", replacement = "", x = sample.data$EH_ID)
exp.des$EH_ID <- gsub(pattern = " ", replacement = "", x = exp.des$EH_ID)
exp.des%>%
  filter(!(EH_ID== "LM0205"))->exp.des ##Eliminate mouse not included in the experiment

sample.data <- merge(sample.data, exp.des, by= "EH_ID", all= TRUE) ##merge sample data with genotype of mice 
rownames(sample.data) <- sample.data$labels
  
###Estimate microbiota density (MD) from Contijoch et al. 2019
### MD = total DNA per sample (µg)/ mg of fresh feces
##considering 30µL of elution volume 
sample.data %>%
  mutate(Total_DNA = (sample.data$Conc_DNA*30)*0.001) -> sample.data ### Add a new variable that will contain total DNA extracted per sample in µg
  
sample.data %>%
  mutate(Microbial_density = sample.data$Total_DNA/(sample.data$fecweight_DNA*1000)) -> sample.data ### Total DNA extracted per sample in µg by feces weight in mg
  
sample.data$labels<- as.vector(sample.data$labels)
rownames(sample.data) <- make.unique(sample.data$labels)

### Table 1:Infection experiment design
sample.data%>%
  dplyr::group_by(EH_ID)%>%
  dplyr::summarise(Weight= mean(weight_dpi0), Sex= Sex, Age= mean(ageAtdpi0expe1a), Strain= Strain)%>%
  dplyr::distinct()-> tab1

min(tab1$Weight)
min(tab1$Age)

max(tab1$Weight)
max(tab1$Age)

mean(tab1$Weight)
mean(tab1$Age)

##7 females, 15 males

#write.csv(tab1, "/SAN/Victors_playground/Eimeria_microbiome/Table_1_Experiment_Setting_Summary.csv", row.names = F) 
##For now later change it to the repo

rm(tab1)

