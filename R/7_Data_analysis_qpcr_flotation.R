### Code to analyse
## 1) Correlation among qPCR and oocyst flotation quantification

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
library("phyloseq")
library("microbiome")
library("grid")
library("ggsci")
library("knitr")
library("kableExtra")


if(!exists("sample.data")){
  source("~/GitProjects/Eimeria_Microbiome/R/1_Data_preparation.R")
}

#if(!exists("data.inf.exp")){
#  source("~/GitProjects/Eimeria_Microbiome/R/6_qPCR_data_preparation.R")
#}
if(!exists("data.inf.exp")){
  data.inf.exp<- read.csv(file="/SAN/Victors_playground/Eimeria_microbiome/qPCR/sample_data_qPCR.csv")
}

setdiff(sample.data$labels, data.inf.exp$labels)

##Get unique labels from qPCR data
data.inf.exp%>%
  select(labels, Qty_mean, Genome_copies_mean,Tm_mean, Infection)%>%
  distinct(labels, .keep_all = TRUE)-> data.inf.exp

### Join all the data in the same dataframe 
sdt<- join(sample.data, data.inf.exp, by="labels") ## Add qPCR data

sdt$dpi<- as.factor(sdt$dpi)

###Let's start plotting and analysing the data!
### 1) Correlation among Eimeria quantification methods

####OPG vs qPCR 
sdt%>%
  ggplot(aes(OPG, Qty_mean))+
  geom_smooth(method = lm)+
  scale_x_log10(name = "log10 Oocyst per gram feces (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 Number of Eimeria Oocysts (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_cor(label.x = 5.5, label.y = 1.5, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 5.5, label.y = 2)+
  stat_cor(label.x = 5.5,  label.y = 1,method = "spearman")+
  annotation_logticks()#-> opgqpcr

sdt%>%
  ggplot(aes(OPG, Genome_copies_mean))+
  geom_smooth(method = lm)+
  scale_x_log10(name = "log10 Oocyst per gram feces (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 Eimeria genome copies (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_cor(label.x = 5.5, label.y = 2.5, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 5.5, label.y = 3)+
  stat_cor(label.x = 5.5,  label.y = 2,method = "spearman")+
  annotation_logticks()-> opgqpcr

###Course of infection 
## Oocysts
sdt%>%
  ggplot(aes(dpi, OPG))+
  geom_boxplot()+
  xlab("Day post infection")+
  scale_y_log10("log10 Oocyst per gram feces (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "black")+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks(sides = "l")-> a

sdt%>%
  dplyr::arrange(dpi)%>%
  ggplot(aes(as.numeric(as.character(dpi)), OPG, colour= EH_ID))+
  xlab("Day post infection")+
  scale_y_log10("log10 Oocyst per gram feces (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10))+
  geom_line(alpha= 0.5)+
  geom_jitter(shape=21, position=position_jitter(0.0), size=2.5, aes(fill= EH_ID), color= "black", alpha= 0.5)+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_smooth(color= "black", method = "loess")+
  annotation_logticks(sides = "l")-> a2

##qPCR
sdt%>%
  ggplot(aes(dpi, Genome_copies_mean))+
  geom_boxplot()+
  xlab("Day post infection")+
  scale_y_log10(name = "log10 Number of Eimeria genome copies (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
  geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "black")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks(sides = "l")-> c

sdt%>%
  dplyr::arrange(dpi)%>%
  ggplot(aes(as.numeric(as.character(dpi)), Genome_copies_mean, colour= EH_ID))+
  xlab("Day post infection")+
  scale_y_log10("log10 Eimeria genome copies (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10))+
  geom_line(alpha=0.5)+
  geom_jitter(shape=21, position=position_jitter(0.0), size=2.5, aes(fill= EH_ID), color= "black", alpha= 0.5)+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_smooth(color= "black", method = "loess")+
  annotation_logticks(sides = "l")-> c2

pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_2.pdf", width = 10, height = 20)
grid.arrange(a,c,opgqpcr, ncol= 1, nrow= 3)
dev.off()

### Predictions Early DNA vs Late Oocysts
sdt%>%
  dplyr::group_by(EH_ID)%>%
  dplyr::arrange(as.factor(as.character(dpi)))%>%
  select(dpi,OPG, Genome_copies_mean)%>%
  filter(EH_ID%in%c("LM0204"))

sdt%>%
  filter(dpi%in%c("3", "7"))%>%
  select(EH_ID, dpi,OPG, Genome_copies_mean)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)->comp2 ##for comparison 

ggpaired(comp2, x= "dpi", y= "Genome_copies_mean",line.color= "gray", line.size= 0.4, color= "dpi")+
  scale_y_log10("log10 Eimeria genome copies (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  xlab("Day post infection")+
  scale_color_npg()+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks(sides = "l")-> dna37

ggpaired(comp2, x= "dpi", y= "OPG",line.color= "gray", line.size= 0.4, color= "dpi")+
  scale_y_log10("log10 Oocysts per gram feces (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  xlab("Day post infection")+
  scale_color_npg()+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks(sides = "l")-> oocyst37

sdt%>%
  dplyr::arrange(dpi)%>%
  filter(dpi%in%c("3"))%>%
  select(EH_ID,dpi,Genome_copies_mean)-> earlyDNA

sdt%>%
  dplyr::arrange(dpi)%>%
  filter(dpi%in%c("7"))%>%
  select(EH_ID,dpi,OPG)-> lateOoc

elod<- join(earlyDNA, lateOoc, by= "EH_ID")

set.seed(2020)
elod%>%
  select(EH_ID,Genome_copies_mean,OPG)%>%
  ggplot(aes(Genome_copies_mean, OPG))+
  geom_smooth(method = lm)+
  scale_x_log10(name = "log10 Eimeria genome copies dpi 3 (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 Oocyst per gram feces dpi 7 (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= EH_ID), color= "black")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_cor(label.x = 3.5, label.y = 6.15, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 3.5, label.y = 6.20)+
  stat_cor(label.x = 3.5,  label.y = 6.085,method = "spearman")+
  annotation_logticks()-> elodp

summary(lm(OPG~Genome_copies_mean,  data = elod)) ## Using just genome copies as predictor
anova(lm(OPG~Genome_copies_mean,  data = elod))

elod<- join(earlyDNA, lateOoc, by= "EH_ID")

##Using linear models for time-series analysis
summary(lm(OPG~Genome_copies_mean,  data = sdt)) ## Using just genome copies as predictor
anova(lm(OPG~Genome_copies_mean,  data = sdt))
summary(lm(OPG~Genome_copies_mean+dpi,  data = sdt)) ## Genome copies and dpi as predictors
anova(lm(OPG~Genome_copies_mean+dpi,  data = sdt))
summary(lm(OPG~Genome_copies_mean*dpi,  data = sdt)) ## Genome copies and dpi and interactions as predictors
anova(lm(OPG~Genome_copies_mean*dpi,  data = sdt))
summary(lm(OPG~Genome_copies_mean+dpi+EH_ID,  data = sdt)) ## Genome copies, dpi and individual as predictors 
anova(lm(OPG~Genome_copies_mean+dpi+EH_ID,  data = sdt))

##Using local regression models (non-parametric approach that fits multiple regressions in local neighborhood)
summary(loess(OPG~Genome_copies_mean,  data = sdt, span = 0.1)) ## Using just genome copies as predictor
predict(loess(OPG~Genome_copies_mean,  data = sdt, span = 0.1))




##Save plots
pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_2.pdf", width = 10, height = 20)
grid.arrange(a,c,opgqpcr, ncol= 1, nrow= 3)
dev.off()

pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_3.pdf", width = 10, height = 10)
grid.arrange(dna37,oocyst37,elodp, widths = c(1, 1), layout_matrix = rbind(c(1, 2), c(3, 3)))
dev.off()
