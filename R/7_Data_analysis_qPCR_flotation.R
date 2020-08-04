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
  select(labels, Genome_copies_mean,Tm_mean, Infection)%>%
  distinct(labels, .keep_all = TRUE)-> data.inf.exp

### Join all the data in the same dataframe 
sdt<- join(sample.data, data.inf.exp, by="labels") ## Add qPCR data

sdt$dpi<- as.factor(sdt$dpi)

###Let's start plotting and analysing the data!
### 1) Correlation among Eimeria quantification methods

####OPG vs qPCR 
sdt%>%
  ggplot(aes(OPG, Genome_copies_mean))+
  geom_smooth(method = lm)+
  scale_x_log10(name = "log10 Oocyst per gram feces (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 Genome copies/µL gDNA (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_cor(label.x = 5.5, label.y = 1.5, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 5.5, label.y = 2)+
  stat_cor(label.x = 5.5,  label.y = 1,method = "spearman")+
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
  scale_y_log10(name = "log10 Genome copies/µL gDNA (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
  geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "black")+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks(sides = "l")-> b

sdt%>%
  dplyr::arrange(dpi)%>%
  ggplot(aes(as.numeric(as.character(dpi)), Genome_copies_mean, colour= EH_ID))+
  xlab("Day post infection")+
  scale_y_log10("log10 Genome copies/µL gDNA (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10))+
  geom_line(alpha=0.5)+
  geom_jitter(shape=21, position=position_jitter(0.0), size=2.5, aes(fill= EH_ID), color= "black", alpha= 0.5)+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_smooth(color= "black", method = "loess")+
  annotation_logticks(sides = "l")-> b2

### Predictions Early DNA vs Late Oocysts
sdt%>%
  filter(dpi%in%c("3","4", "6", "7", "8", "9", "10"))%>%
  select(EH_ID, dpi,OPG, Genome_copies_mean)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  ggplot(aes(x= dpi, y= Genome_copies_mean))+
  scale_y_log10("log10 Genome copies/µL gDNA (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot(aes(color= dpi))+
  geom_point(aes(color=dpi))+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray")+
  scale_color_npg()+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks(sides = "l")+
  #stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "3", paired = TRUE, na.rm = TRUE)+
  stat_compare_means(method =  "anova")-> dna3_10

  compare_means(formula = Genome_copies_mean~dpi, method = "wilcox.test", ref.group = "3",paired = TRUE, data = sdt)

  sdt%>%
    filter(dpi%in%c("3","4", "6", "7", "8", "9", "10"))%>%
    select(EH_ID, dpi,OPG, Genome_copies_mean)%>%
    dplyr::arrange(EH_ID)%>%
    dplyr::arrange(dpi)%>% ##for comparison 
    ggplot(aes(x= dpi, y= OPG))+
    scale_y_log10("log10 Oocysts per gram feces (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    geom_boxplot(aes(color= dpi))+
    geom_point(aes(color=dpi))+
    xlab("Day post infection")+
    geom_line(aes(group = EH_ID), color= "gray")+
    scale_color_npg()+
    labs(tag= "B)")+
    theme_bw()+
    theme(text = element_text(size=16))+
    annotation_logticks(sides = "l")+
    stat_compare_means("anova")-> oocyst3_10

compare_means(formula = OPG~dpi, method = "wilcox.test", ref.group = "3",paired = TRUE, data = sdt)

### Do Eimeria DNA in early dpi can predict the oocyst at pick day?
##for comparison between early DNA and Pick oocysts
sdt%>%
  filter(dpi%in%c("3","6"))%>%
  select(EH_ID, dpi,OPG, Genome_copies_mean)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% 
  ggplot(aes(x= dpi, y= OPG))+
  scale_y_log10("log10 Oocysts per gram feces (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot(aes(color= dpi))+
  geom_point(aes(color=dpi))+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray")+
  scale_color_npg()+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks(sides = "l")+
  coord_cartesian(ylim = c(100, 10000000))-> oocyst36

sdt%>%
  filter(dpi%in%c("3","6"))%>%
  select(EH_ID, dpi,OPG, Genome_copies_mean)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% 
  ggplot(aes(x= dpi, y= Genome_copies_mean))+
  scale_y_log10("log10 Genome copies/µL gDNA (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot(aes(color= dpi))+
  geom_point(aes(color=dpi))+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray")+
  scale_color_npg()+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks(sides = "l")+
  coord_cartesian(ylim = c(100, 10000000))-> dna36

sdt%>%
  dplyr::arrange(dpi)%>%
  filter(dpi%in%c("3"))%>%
  select(EH_ID,dpi,Genome_copies_mean)-> earlyDNA

sdt%>%
  dplyr::arrange(dpi)%>%
  filter(dpi%in%c("6"))%>%
  select(EH_ID,dpi,OPG)-> pickOoc

elop<- join(earlyDNA, pickOoc, by= "EH_ID")
rm(pickOoc)

set.seed(2020)
elop%>%
  select(EH_ID,Genome_copies_mean,OPG)%>%
  ggplot(aes(Genome_copies_mean, OPG))+
  geom_smooth(method = lm)+
  scale_x_log10(name = "log10 Genome copies/µL gDNA dpi 3 (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 Oocyst per gram feces dpi 6 (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= EH_ID), color= "black")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_cor(label.x = 2.25, label.y = 4.5, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 2.25, label.y = 5.0)+
  stat_cor(label.x = 2.25,  label.y = 4.75,method = "spearman")+
  annotation_logticks()+
  coord_cartesian(ylim = c(10000, 10000000))-> elopp

summary(lm(OPG~Genome_copies_mean,  data = elop)) ## Using just genome copies as predictor
anova(lm(OPG~Genome_copies_mean,  data = elod))

### Do Eimeria DNA at early dpi can predict the oocyst at late dpi?
sdt%>%
  filter(dpi%in%c("3","7"))%>%
  select(EH_ID, dpi,OPG, Genome_copies_mean)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>%
  ggplot(aes(x= dpi, y= OPG))+
  scale_y_log10("log10 Oocysts per gram feces (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot(aes(color= dpi))+
  geom_point(aes(color=dpi))+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray")+
  scale_color_npg()+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks(sides = "l")+
  coord_cartesian(ylim = c(100, 10000000))-> oocyst37

sdt%>%
  filter(dpi%in%c("3","7"))%>%
  select(EH_ID, dpi,OPG, Genome_copies_mean)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>%
  ggplot(aes(x= dpi, y= Genome_copies_mean))+
  scale_y_log10("log10 Genome copies/µL gDNA (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot(aes(color= dpi))+
  geom_point(aes(color=dpi))+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray")+
  scale_color_npg()+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks(sides = "l")+
  coord_cartesian(ylim = c(100, 10000000))-> dna37

sdt%>%
  dplyr::arrange(dpi)%>%
  filter(dpi%in%c("7"))%>%
  select(EH_ID,dpi,OPG)-> lateOoc ## DPI where most of samples still are excreting oocyst after the pick 

elod<- join(earlyDNA, lateOoc, by= "EH_ID")

set.seed(2020)
elod%>%
  select(EH_ID,Genome_copies_mean,OPG)%>%
  ggplot(aes(Genome_copies_mean, OPG))+
  geom_smooth(method = lm)+
  scale_x_log10(name = "log10 Genome copies/µL gDNA dpi 3 (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 Oocyst per gram feces dpi 7 (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= EH_ID), color= "black")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_cor(label.x = 2.25, label.y = 4.75, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 2.25, label.y = 5.0)+
  stat_cor(label.x = 2.25,  label.y = 4.5,method = "spearman")+
  annotation_logticks()+
  coord_cartesian(ylim = c(10000, 10000000))-> elodp

summary(lm(OPG~Genome_copies_mean,  data = elod)) ## Using just genome copies as predictor
anova(lm(OPG~Genome_copies_mean,  data = elod))

### Do Eimeria DNA at early dpi can predict the oocyst before pick dpi?
sdt%>%
  filter(dpi%in%c("3","5"))%>%
  select(EH_ID, dpi,OPG, Genome_copies_mean)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>%
  ggplot(aes(x= dpi, y= OPG))+
  scale_y_log10("log10 Oocysts per gram feces (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot(aes(color= dpi))+
  geom_point(aes(color=dpi))+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray")+
  scale_color_npg()+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks(sides = "l")+
  coord_cartesian(ylim = c(100, 10000000))-> oocyst35

sdt%>%
  filter(dpi%in%c("3","5"))%>%
  select(EH_ID, dpi,OPG, Genome_copies_mean)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>%
  ggplot(aes(x= dpi, y= Genome_copies_mean))+
  scale_y_log10("log10 Genome copies/µL gDNA (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot(aes(color= dpi))+
  geom_point(aes(color=dpi))+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray")+
  scale_color_npg()+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks(sides = "l")+
  coord_cartesian(ylim = c(100, 10000000))-> dna35

sdt%>%
  dplyr::arrange(dpi)%>%
  filter(dpi%in%c("5"))%>%
  select(EH_ID,dpi,OPG)-> beforeOoc ## DPI where most of samples still are excreting oocyst before the pick 

elob<- join(earlyDNA, beforeOoc, by= "EH_ID")

set.seed(2020)
elob%>%
  select(EH_ID,Genome_copies_mean,OPG)%>%
  ggplot(aes(Genome_copies_mean, OPG))+
  geom_smooth(method = lm)+
  scale_x_log10(name = "log10 Genome copies/µL gDNA dpi 3 (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log10 Oocyst per gram feces dpi 5 (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= EH_ID), color= "black")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_cor(label.x = 2.25, label.y = 6.25, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 2.25, label.y = 6.5)+
  stat_cor(label.x = 2.25,  label.y = 6.0,method = "spearman")+
  annotation_logticks()+
  coord_cartesian(ylim = c(10000, 10000000))-> elobp

summary(lm(OPG~Genome_copies_mean,  data = elob)) ## Using just genome copies as predictor
anova(lm(OPG~Genome_copies_mean,  data = elob))

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
pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_3.pdf", width = 10, height = 20)
grid.arrange(a,b,opgqpcr, ncol= 1, nrow= 3)
dev.off()

pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_4.pdf", width = 20, height = 10)
grid.arrange(oocyst36, dna36, elopp, oocyst37, dna37, elodp, ncol= 3, nrow= 2)
dev.off()

pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_4.1.pdf", width = 10, height = 10)
grid.arrange(oocyst36, dna36, elopp, widths = c(1, 1), layout_matrix = rbind(c(1, 2), c(3, 3)))
dev.off()

pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_4.2.pdf", width = 10, height = 10)
grid.arrange(oocyst37, dna37, elodp, widths = c(1, 1), layout_matrix = rbind(c(1, 2), c(3, 3)))
dev.off()

pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_4.3.pdf", width = 10, height = 10)
grid.arrange(oocyst35, dna35, elobp, widths = c(1, 1), layout_matrix = rbind(c(1, 2), c(3, 3)))
dev.off()

pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_4.4.pdf", width = 25, height = 8)
grid.arrange(elobp, elopp, elodp, ncol= 3, nrow= 1)
dev.off()