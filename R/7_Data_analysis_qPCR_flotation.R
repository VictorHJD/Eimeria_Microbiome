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
  geom_smooth(method = lm, col= "black")+
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

summary(lm(OPG~Genome_copies_mean,  data = sdt, na.action = na.exclude)) ## Using just genome copies as predictor

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
  filter(dpi%in%c("0","1","2","3","4", "6", "7", "8", "9", "10"))%>%
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
    filter(dpi%in%c("0","1","2","3","4", "6", "7", "8", "9", "10"))%>%
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
  dplyr::select(EH_ID, dpi,OPG, Genome_copies_mean)%>%
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
  dplyr::select(EH_ID,Genome_copies_mean,OPG)%>%
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
  dplyr::select(EH_ID,Genome_copies_mean,OPG)%>%
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
  dplyr::select(EH_ID,Genome_copies_mean,OPG)%>%
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

##Using lm for time-series analysis

summary(glm(OPG~Genome_copies_mean*dpi,  data = sdt, na.action = na.exclude)) ## Using just genome copies as predictor
anova(glm(OPG~Genome_copies_mean*dpi,  data = sdt, na.action = na.exclude))

## Create Time-Series
require("reshape")
sdt%>%
  dplyr::select(EH_ID, dpi, Genome_copies_mean)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>%
  na.exclude()-> dna
dna<- reshape(dna, idvar = "EH_ID", timevar = "dpi", direction = "wide")

sdt%>%
  dplyr::select(EH_ID, dpi, OPG)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>%
  na.exclude()-> oocysts
oocysts<- reshape(oocysts, idvar = "EH_ID", timevar = "dpi", direction = "wide")

ts.data<- join(oocysts, dna, by= "EH_ID")

ts.data%>%
  dplyr::rowwise()%>%
  dplyr::mutate(Sum_Oocysts= sum(c(OPG.0,OPG.3,OPG.4,OPG.5,OPG.6,OPG.7,OPG.8,OPG.9,OPG.10)))-> ts.data

### Models
##Check distribution of data 
require("fitdistrplus")
require("logspline")

x <- ts.data$Sum_Oocysts
x <- ts.data$OPG.6
x <- round(na.omit(ts.data$Genome_copies_mean.0))
x <- ts.data$Genome_copies_mean.1
x <- ts.data$Genome_copies_mean.2
x <- ts.data$Genome_copies_mean.3
x <- ts.data$Genome_copies_mean.4
x <- ts.data$Genome_copies_mean.5
x <- ts.data$Genome_copies_mean.6
x <- ts.data$Genome_copies_mean.7
x <- ts.data$Genome_copies_mean.8
x <- ts.data$Genome_copies_mean.9
x <- ts.data$Genome_copies_mean.10

plotdist(x, histo = TRUE, demp = TRUE)
descdist(x, boot = 1000)
fit.norm <- fitdist(x, "norm")
#fit.nbinom <- fitdist(x, "nbinom")
plot(fit.norm)
#plot(fit.nbinom)
## Genome copies have a negative binomial distribution use glm.nb
##Total OPGs during infection are predicted by DNA at different dpi? Genome copies per dpi as individual predictors
sum.opg <- glm.nb(formula = Sum_Oocysts~ Genome_copies_mean.0+
      Genome_copies_mean.1+
      Genome_copies_mean.2+
      Genome_copies_mean.3+
      Genome_copies_mean.4+
      Genome_copies_mean.5+
      Genome_copies_mean.6+
      Genome_copies_mean.7+
      Genome_copies_mean.8+
      Genome_copies_mean.9+
      Genome_copies_mean.10, data = ts.data, na.action = na.exclude)

summary(sum.opg)
plot(sum.opg)

library(sjPlot)
library(sjmisc)
library(sjlabelled)
tab_model(sum.opg)

##extract p values for bonferroni correction
p.sum.opg<- as.data.frame(coef(summary(sum.opg))[,'Pr(>|z|)'])
colnames(p.sum.opg)<- "P_unadjusted"
p.sum.opg$P_adjusted<-p.adjust(p.sum.opg$`P_unadjusted`, method = "bonferroni")

##OPG at pick of infection are predicted by DNA at different dpi? Genome copies per dpi as individual predictors 
dpi6.opg <- glm.nb(formula = OPG.6~ Genome_copies_mean.0+
                 Genome_copies_mean.1+
                 Genome_copies_mean.2+
                 Genome_copies_mean.3+
                 Genome_copies_mean.4+
                 Genome_copies_mean.5+
                 Genome_copies_mean.6+
                 Genome_copies_mean.7+
                 Genome_copies_mean.8+
                 Genome_copies_mean.9+
                 Genome_copies_mean.10, data = ts.data, na.action = na.exclude)

summary(dpi6.opg)
plot(dpi6.opg)
tab_model(dpi6.opg)

##extract p values for bonferroni correction
p.dpi6.opg<- as.data.frame(coef(summary(dpi6.opg))[,'Pr(>|z|)'])
colnames(p.dpi6.opg)<- "P_unadjusted"
p.dpi6.opg$P_adjusted<-p.adjust(p.dpi6.opg$`P_unadjusted`, method = "bonferroni")

##Using local regression models (non-parametric approach that fits multiple regressions in local neighborhood)
summary(loess(OPG~Genome_copies_mean,  data = sdt, span = 0.1)) ## Using just genome copies as predictor
predict(loess(OPG~Genome_copies_mean,  data = sdt, span = 0.1))

## DNA as a predictor of weightloss
sdt%>%
  ggplot(aes(Genome_copies_mean, (weight/weight_dpi0)*100))+
  geom_smooth(method = lm)+
  scale_y_continuous(name = "Relative weight loss to 0 dpi (%)")+
  scale_x_log10(name = "log10 Genome copies/µL gDNA (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  #annotation_logticks("b")+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_cor(label.x = 2.0, label.y = 84, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 2.0, label.y = 88)+
  stat_cor(label.x = 2.0,  label.y = 86,method = "spearman")-> wlqpcr

## Weight loss
sdt%>%
  ggplot(aes(dpi, (weight/weight_dpi0)*100))+
  geom_boxplot()+
  xlab("Day post infection")+
  scale_y_continuous("Relative weight loss at 0 dpi (%)")+
  geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "black")+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))-> ai

### Predictions Early DNA vs Late Oocysts
sdt%>%
  filter(dpi%in%c("3","4", "6", "7", "8", "9", "10"))%>%
  select(EH_ID, dpi,OPG, Genome_copies_mean, weight, weight_dpi0)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% ##for comparison 
  ggplot(aes(x= dpi, y= (weight/weight_dpi0)*100))+
  scale_y_continuous("Relative weight loss at 0 dpi (%)")+
  geom_boxplot(aes(color= dpi))+
  geom_point(aes(color=dpi))+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray")+
  scale_color_npg()+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  #stat_compare_means(label= "p.signif", method = "wilcox.test", ref.group = "3", paired = TRUE, na.rm = TRUE)+
  stat_compare_means(method =  "anova")-> weight3_10

compare_means(formula = Genome_copies_mean~dpi, method = "wilcox.test", ref.group = "3",paired = TRUE, data = sdt)

##for comparison between early DNA and Pick oocysts
sdt%>%
  filter(dpi%in%c("3","6"))%>%
  select(EH_ID, dpi,OPG, Genome_copies_mean, weight, weight_dpi0)%>%
  dplyr::arrange(EH_ID)%>%
  dplyr::arrange(dpi)%>% 
  ggplot(aes(x= dpi, y= (weight/weight_dpi0)*100))+
  scale_y_continuous("Relative weight loss at 0 dpi (%)")+
  geom_boxplot(aes(color= dpi))+
  geom_point(aes(color=dpi))+
  xlab("Day post infection")+
  geom_line(aes(group = EH_ID), color= "gray")+
  scale_color_npg()+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))-> weight36

sdt%>%
  dplyr::arrange(dpi)%>%
  filter(dpi%in%c("6"))%>%
  select(EH_ID,dpi, weight, weight_dpi0, weightloss)-> pickweight

elopw<- join(earlyDNA, pickweight, by= "EH_ID")
rm(pickweight)

set.seed(2020)
elopw%>%
  select(EH_ID,Genome_copies_mean, weight, weight_dpi0, weightloss)%>%
  ggplot(aes(Genome_copies_mean, (weight/weight_dpi0)*100))+
  geom_smooth(method = lm)+
  scale_x_log10(name = "log10 Genome copies/µL gDNA dpi 3 (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous("Relative Weight loss at 6 to 0 dpi")+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= EH_ID), color= "black")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_cor(label.x = 2.25, label.y = 83.5, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 2.25, label.y = 87.5)+
  stat_cor(label.x = 2.25,  label.y = 85,method = "spearman")+
  annotation_logticks("b")+
  coord_cartesian(ylim = c(10000, 10000000))-> elopp

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