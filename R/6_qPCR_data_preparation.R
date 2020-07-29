## Code for 
## 1) Standard curve of qPCR Eimeria
## 2) Determine Eimeria amount for infection experiment samples
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


Standards<- T
Unknowns<- F
Mock<- F
Infexp<- T

##Load data
if(Standards){
  data.std<- read.csv("/SAN/Victors_playground/Eimeria_microbiome/qPCR/Eimeria_quantification_Std_Curve_data.csv")
  data.std%>%
    dplyr::mutate(Genome_copies= Oocyst_count*8)-> data.std
  
  ##Define numeric and factor variables 
  num.vars <- c("Ct", "Ct_mean", "Sd_Ct", "Qty", "Qty_mean", "Sd_Qty", "Oocyst_count", "Feces_weight", "Qubit", "NanoDrop", "Beads_weight", "Tm", "Genome_copies")
  fac.vars <- c("Well", "Sample.Name", "Detector", "Task",  "Std_series","Date", "Operator", "Cycler", "Parasite", "Sample_type", "Extraction")  
  data.std[, num.vars] <- apply(data.std[, num.vars], 2, as.numeric)
  data.std[, fac.vars] <- apply(data.std[, fac.vars], 2, as.factor)
}

if(Unknowns){
  data.unk<-read.csv("/SAN/Victors_playground/Eimeria_microbiome/qPCR/Eimeria_quantification_Sample_data.csv")
  
  ##Define numeric and factor variables 
  num.vars <- c("Ct", "Ct_mean", "Sd_Ct", "Qty", "Qty_mean", "Sd_Qty", "Oocyst_count", "Feces_weight", "Qubit", "NanoDrop", "Beads_weight", "Tm", 
                "Oocyst_1", "Oocyst_2", "Oocyst_3", "Oocyst_4", "Oocyst_5", "Oocyst_6", "Oocyst_7", "Oocyst_8", "Dilution_factor", "Volume")
  fac.vars <- c("Well", "Sample.Name", "Detector", "Task", "Date", "Operator", "Cycler", "Parasite", "Sample_type", "Extraction")  
  data.unk[, num.vars] <- apply(data.unk[, num.vars], 2, as.numeric)
  data.unk[, fac.vars] <- apply(data.unk[, fac.vars], 2, as.factor)
}

if(Infexp){
  data.inf<-read.csv("/SAN/Victors_playground/Eimeria_microbiome/qPCR/Eimeria_quantification_Inf_exp_data.csv")
  data.inf%>%
    select(Content, Sample, Plate_number, Cq, Melt_Temperature)%>%
    dplyr::rename(Ct= Cq, labels= Sample, Task= Content, Tm= Melt_Temperature)-> data.inf
  
  ##Define numeric and factor variables 
  num.vars <- c("Ct", "Tm")
  fac.vars <- c("labels", "Task", "Plate_number")  
  data.inf[, num.vars] <- apply(data.inf[, num.vars], 2, as.numeric)
  data.inf[, fac.vars] <- apply(data.inf[, fac.vars], 2, as.factor)
  
  ### Add sample information
  if(!exists("sample.data")){
    source("~/GitProjects/Eimeria_Microbiome/R/1_Data_preparation.R")
  }
}

rm(fac.vars, num.vars)
####### Standard curves #######
set.seed(2020)

data.std%>%
  select(Sample.Name,Task,Std_series,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm)%>%
  filter(Task=="Standard" & Cycler=="Eppendorf")%>%
  dplyr::group_by(Parasite)%>%
  ggplot(aes(Qty, Ct))+
  scale_x_log10("log 10 Number of Eimeria Oocyst (Flotation)")+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Std_series), color= "black")+
  geom_smooth(aes(color= Std_series, fill= Std_series), method = "lm")+            
  stat_cor(aes(color = Std_series), label.x = 4,  label.y = c(30, 31),method = "spearman")+
  stat_cor(label.x = 4, label.y = c(28, 30),aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Std_series))+        # Add correlation coefficient
  stat_regline_equation(aes(color = Std_series), label.x = 3, label.y = c(30, 31))+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  labs(tag = "A)")+
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none") #-> A

summary(lm(formula = log10(Oocyst_count)~Ct, data = subset(data.std, Cycler=="Eppendorf")))

data.std%>%
  select(Sample.Name,Task,Std_series,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm)%>%
  filter(Task=="Standard" & Cycler=="ABI")%>%
  dplyr::group_by(Parasite)%>%
  ggplot(aes(Qty, Ct))+
  scale_x_log10("log 10 Number of Eimeria Oocyst (Flotation)")+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Std_series), color= "black")+
  geom_smooth(aes(color= Std_series, fill= Std_series), method = "lm")+            
  stat_cor(aes(color = Std_series), label.x = 3,  label.y = c(30, 30, 30),method = "spearman")+
  stat_cor(label.x = 3, label.y = c(28, 28, 28),aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Std_series))+        # Add correlation coefficient
  stat_regline_equation(aes(color = Std_series), label.x = 3, label.y = c(32, 32,32))+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  labs(tag = "B)")+
  theme_bw() +
  facet_grid(cols = vars(Std_series))+
  theme(text = element_text(size=20), legend.position = "none") #-> B

summary(lm(formula = log10(Oocyst_count)~Ct, data = subset(data.std, Cycler=="ABI" & Task== "Standard")))

data.std%>%
  select(Sample.Name,Task,Std_series,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm)%>%
  filter(Task=="Standard" & Parasite=="E_falciformis" & Cycler=="BioRad")%>%
  dplyr::group_by(Sample.Name)%>%
  ggplot(aes(Qty, Ct))+
  scale_x_log10("log 10 Number of Eimeria Oocyst (Flotation)")+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Std_series), color= "black")+
  geom_smooth(aes(color= Std_series, fill= Std_series), method = "lm")+            
  stat_cor(aes(color = Std_series), label.x = 3,  label.y = c(34, 34),method = "spearman")+
  stat_cor(label.x = 3, label.y = c(32, 32),aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Std_series))+        # Add correlation coefficient
  stat_regline_equation(aes(color = Std_series), label.x = 3, label.y = c(36, 36))+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  labs(tag = "C)")+
  theme_bw() +
  facet_grid(cols = vars(Std_series))+
  theme(text = element_text(size=20), legend.position = "none") #-> C

summary(lm(formula = log10(Oocyst_count)~Ct, data = subset(data.std, Cycler=="BioRad" & Task== "Standard")))

#pdf(file = "~/AA_Microbiome/qPCR_Eimeria_Quantification/figures/Standard_curves.pdf", width = 10, height = 15)
#grid.arrange(A, B, C)
#dev.off()
#rm(A,B,C)

##Final standard curve for Inf. experiment samples 
## Eimeria Oocysts
data.std%>%
  select(Sample.Name,Task,Std_series,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm)%>%
  filter(Task=="Standard" & Parasite=="E_falciformis" & Cycler=="BioRad")%>%
  dplyr::group_by(Sample.Name)%>%
  #filter(Sample.Name!="Eimeria_10_0")%>%
  ggplot(aes(Qty, Ct))+
  scale_x_log10("log 10 Number of Eimeria Oocyst (Flotation)")+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Task), color= "black")+
  geom_smooth(aes(color= Task, fill= Task),method = "lm")+            
  stat_cor(aes(color= Task), label.x = 5,  label.y = 34, method = "spearman")+
  stat_cor(label.x = 5, label.y = 32, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Task))+        # Add correlation coefficient
  stat_regline_equation(aes(color = Task), label.x = 5, label.y = 36)+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  labs(tag = "A)")+
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none")-> D

## Eimeria genome copies
data.std%>%
  select(Sample.Name,Task,Std_series,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm, Genome_copies)%>%
  filter(Task=="Standard" & Parasite=="E_falciformis" & Cycler=="BioRad")%>%
  dplyr::group_by(Sample.Name)%>%
  #filter(Sample.Name!="Eimeria_10_0")%>%
  ggplot(aes(Genome_copies, Ct))+
  scale_x_log10("log 10 Number of Eimeria genome copies")+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Task), color= "black")+
  geom_smooth(aes(color= Task, fill= Task),method = "lm")+            
  stat_cor(aes(color= Task), label.x = 5,  label.y = 34, method = "spearman")+
  stat_cor(label.x = 5, label.y = 32, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Task))+        # Add correlation coefficient
  stat_regline_equation(aes(color = Task), label.x = 5, label.y = 36)+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  labs(tag = "B)")+
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none")-> D2

### Figure 1 Final Standard curves 
pdf(file = "~/AA_Microbiome/Figures/Figure_1.pdf", width = 10, height = 10)
grid.arrange(D, D2)
dev.off()

rm(D,D2)
###### Intersample variation #####

data.unk%>%
  select(Sample.Name, Task, Ct,Qty,Cycler,Parasite, Sample_type, Tm, Oocyst_1, Oocyst_2, Oocyst_3, Oocyst_4, Oocyst_5,Oocyst_6, Oocyst_7, Oocyst_8, Dilution_factor, Volume)%>%
  filter(Sample_type=="Oocysts" & Task=="Unknown")%>%
  dplyr::group_by(Sample.Name)%>%
  dplyr::mutate(N= n())%>%
  dplyr::mutate(Total_oocysts= (sum(Oocyst_1, Oocyst_2, Oocyst_3, Oocyst_4, Oocyst_5,Oocyst_6, Oocyst_7, Oocyst_8))/N)%>%
  dplyr::mutate(Oocyst_count= (((Total_oocysts*10000)/8)*Dilution_factor)*Volume)%>%
  ggplot(aes(Oocyst_count,Qty), geom=c("point", "smooth"))+
  scale_x_log10(name = "log10 Number of Eimeria Oocysts (Flotation)")+
  scale_y_log10(name = "log10 Number of Eimeria Oocysts (qPCR)")+ 
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Task), color= "black")+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  theme_bw() +
  theme(legend.text=element_text(size=20)) +
  theme(legend.key.size = unit(3,"line")) +
  geom_smooth(aes(color= Task, fill= Task), method = "lm")+            
  stat_cor(aes(color = Task), label.x = 5.5,  label.y = 4,method = "spearman")+
  stat_cor(label.x = 5.5, label.y = 3.75,aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Task))+        # Add correlation coefficient
  stat_regline_equation(aes(color = Task), label.x = 5, label.y = 4)+
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(text = element_text(size=20),legend.position = "none")+
  labs(tag = "A)") #-> E

summary(lm(formula = log10(Oocyst_count)~log10(Qty), data = subset(data.unk, Sample_type== "Oocysts")))

#pdf(file = "~/AA_Microbiome/qPCR_Eimeria_Quantification/figures/Intersample_Variation.pdf", width = 10, height = 10)
#E
#dev.off()
#rm(E)
########## Mock samples Experiment #########

if(Mock){
  set.seed(2020)

  ##Standard curve for this experiment 
  data.std%>%
    select(Sample.Name,Task,Std_series,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm, Date)%>%
    filter(Task=="Standard" & Cycler=="ABI" & Date=="300620")%>%
    ggplot(aes(Qty, Ct), geom=c("point", "smooth"))+
    scale_x_log10(name = "Number of Eimeria oocysts")+
    geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Std_series), color= "black")+
    geom_smooth(aes(color= Std_series, fill= Std_series), method = "lm")+            
    stat_cor(aes(color = Std_series), label.x = 4,  label.y = c(35, 31),method = "spearman")+
    stat_cor(label.x = 4, label.y = c(34, 30),aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Std_series))+        # Add correlation coefficient
    stat_regline_equation(aes(color = Std_series), label.x = 3, label.y = c(35, 31))+
    stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
    labs(tag = "A)")+
    theme_bw() +
    theme(text = element_text(size=20), legend.position = "none")
  
  set.seed(2020)
  data.std%>%
    select(Sample.Name,Task,Std_series,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm, Date)%>%
    filter(Task=="Standard" & Cycler=="ABI" & Date=="300620")%>%
    ggplot(aes(x = Qty, y = Ct), geom=c("point", "smooth")) +
    scale_x_log10(name = "Number of Eimeria oocysts") +
    scale_y_continuous(name = "Ct")+ 
    geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Std_series), color= "black")+
    stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
    theme_bw() +
    theme(legend.text=element_text(size=20)) +
    theme(legend.key.size = unit(3,"line")) +
    geom_smooth(method = "lm", se = T, col = "red") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme(text = element_text(size=20),legend.position = "none")+
    labs(tag = "B)")+
    stat_cor(label.x = log10(100000), label.y = 31, method = "spearman")+
    stat_cor(label.x = log10(100000), label.y = 30, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
    stat_regline_equation(label.x = log10(100000), label.y = 32)
  
  #pdf(file = "~/AA_Microbiome/qPCR_Eimeria_Quantification/figures/ABI_Cycler/Fig_1_Mock.pdf", width = 15, height = 8)
  #grid.arrange(A, B, widths = c(2, 2), layout_matrix = rbind(c(1, 2)))
  #dev.off()
  
  ##Mean comparison standars against NTC
  #comp<- compare_means(Ct ~ Sample.Name,  data = subset(data.mock, Task%in%c("Standard", "NTC")), method = "t.test", p.adjust.method = "bonferroni")
  
  set.seed(2020)
  data.std%>%
    select(Sample.Name,Task,Std_series,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm, Date)%>%
    filter(Task%in%c("Standard", "NTC") & Cycler=="ABI" & Date=="300620")%>%
    ggplot(aes(x = Sample.Name, y = Ct)) +
    scale_x_discrete(name = "Standard", labels= c("Eimeria_10_0"= "Oocysts 10⁰", "Eimeria_10_1"= "Oocysts 10¹",
                                                  "Eimeria_10_2"= "Oocysts 10²", "Eimeria_10_3"= "Oocysts 10³",
                                                  "Eimeria_10_4"= "Oocysts 10⁴", "Eimeria_10_5"= "Oocysts 10⁵",
                                                  "Eimeria_10_6"= "Oocysts 10⁶", "H2O"= "NTC")) +
    scale_y_continuous(name = "Ct")+ 
    geom_jitter(shape=21, position=position_jitter(0.2), color= "black", aes(size= 25, fill= Std_series))+
    theme_bw() +
    theme(text = element_text(size=16),legend.position = "none")+
    theme(axis.text.x = element_text(angle=90))+
    stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
    labs(tag = "A)")+
    geom_hline(yintercept = 30, linetype = 2)+
    stat_compare_means(method = "anova",
                       aes(label = paste0(..method.., ",\n","p=",..p.format..)), label.y= 34, label.x = 7)+
    stat_compare_means(label = "p.signif", method = "t.test",ref.group = "H2O", 
                       label.y = c(34, 31, 27, 24, 21, 17, 15, 0)) ##Determine that 10^0 meassurments are basically like NTC
  
  ##Compair mock samples qPCR estimation with real oocyst count
  set.seed(2020)
  data.unk%>%
    select(Sample.Name, Task, Ct,Qty,Cycler,Parasite, Sample_type, Feces_weight, Extraction, Oocyst_count)%>%
    filter(Sample_type=="Feces" & Task=="Unknown")%>%
    ggplot(aes(x = Oocyst_count, y = Qty, color=Extraction), geom=c("point", "smooth")) +
    scale_x_log10(name = "log10 Number of Eimeria Oocysts (Flotation)")+
    scale_y_log10(name = "log10 Number of Eimeria Oocysts (qPCR)")+ 
    geom_jitter(shape=21, position=position_jitter(0.2), color= "black", aes(size= 25, fill= Extraction))+
    theme_bw() +
    geom_smooth(aes(color= Extraction, fill= Extraction), method = "lm")+
    facet_grid(cols = vars(Extraction))+
    stat_cor(aes(color = Extraction), label.x = log10(1000),  label.y = log10(100000),method = "spearman")+
    stat_cor(label.x = log10(1000), label.y = log10(50000), aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Extraction))+        # Add correlation coefficient
    stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
    #stat_regline_equation(aes(color = Std_series), label.x = 3, label.y = c(35, 31))
    theme(text = element_text(size=20),legend.position = "none")+
    labs(tag = "A)")
  
  #pdf(file = "~/AA_Microbiome/qPCR_Eimeria_Quantification/figures/ABI_Cycler/Fig_3_Mock.pdf", width = 15, height = 8)
  #E
  #dev.off()
  
  #summary(lm(formula = log10(Amount_SYBR_Copies)~Ct_SYBR, data = data.epp))
  
  ##Standard curve and ceramic beads data
  data.std%>%
    select(Sample.Name,Task,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm, Sample_type,Feces_weight,Extraction,Date)%>%
    filter(Task%in%c("Standard", "NTC") & Cycler=="ABI" & Date=="300620")-> Std.mock
  
  data.unk%>%
    select(Sample.Name,Task,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm,Sample_type,Feces_weight,Extraction,Date)%>%
    filter(Sample_type=="Feces" & Task=="Unknown" & Extraction!="Glass_beads")%>%
    bind_rows(Std.mock)-> data.mock
  
  rm(Std.mock)
  
  data.mock%>%
    select(Sample.Name,Task,Qty,Ct,Oocyst_count)%>%  
    filter(Task%in%c("Standard", "Unknown"))%>%
    ggplot(aes(x = Oocyst_count, y = Ct), geom=c("point", "smooth")) +
    scale_x_log10(name = "Number of Eimeria oocysts") +
    scale_y_continuous(name = "Ct")+ 
    geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Task), color= "black")+
    stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
    theme_bw() +
    theme(legend.text=element_text(size=20)) +
    theme(legend.key.size = unit(3,"line")) +
    geom_smooth(aes(color= Task, fill= Task), method = "lm")+            
    stat_cor(aes(color = Task), label.x = 4,  label.y = c(35, 31),method = "spearman")+
    stat_cor(label.x = 4, label.y = c(34, 30),aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Task))+        # Add correlation coefficient
    stat_regline_equation(aes(color = Task), label.x = 3, label.y = c(35, 31))+
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme(text = element_text(size=20),legend.position = "none")+
    labs(tag = "A)")
  
  data.mock%>%
    select(Sample.Name,Task,Qty,Ct,Oocyst_count)%>%  
    filter(Task%in%c("Standard", "Unknown"))%>%
    ggplot(aes(x = Oocyst_count, y = Qty), geom=c("point", "smooth")) +
    scale_x_log10(name = "log10 Number of Eimeria Oocysts (Flotation)")+
    scale_y_log10(name = "log10 Number of Eimeria Oocysts (qPCR)")+ 
    geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Task), color= "black")+
    stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
    theme_bw() +
    theme(legend.text=element_text(size=20)) +
    theme(legend.key.size = unit(3,"line")) +
    geom_smooth(aes(color= Task, fill= Task), method = "lm")+            
    stat_cor(aes(color = Task), label.x = 2,  label.y = c(6, 5),method = "spearman")+
    stat_cor(label.x = 2, label.y = c(5.5, 4.5),aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Task))+        # Add correlation coefficient
    stat_regline_equation(aes(color = Task), label.x = 1, label.y = c(6, 5))+
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme(text = element_text(size=20),legend.position = "none")+
    labs(tag = "B)")
  
  #pdf(file = "~/AA_Microbiome/qPCR_Eimeria_Quantification/figures/ABI_Cycler/Fig_4_Mock.pdf", width = 10, height = 15)
  #grid.arrange(f, G)
  #dev.off()
  
  ##Mean comparison unknown samples against NTC
  compare_means(Ct ~ Sample.Name,  data = subset(data.mock, Task%in%c("Unknown", "NTC")), method = "t.test", p.adjust.method = "bonferroni")
  
  set.seed(2020)
  data.mock%>%
    select(Sample.Name,Task,Qty,Ct,Oocyst_count,Sample_type)%>%  
    filter(Task%in%c("Unknown", "NTC"))%>%
    ggplot(aes(x = Sample.Name, y = Ct)) +
    scale_x_discrete(name = "No. of oocysts", labels= c("E_falciformis_GB_01"= "1.115x10⁶", "E_falciformis_GB_02"= "1.115x10⁵",
                                                        "E_falciformis_GB_03"= "1.115x10⁴", "E_falciformis_GB_04"= "1.115x10³",
                                                        "E_falciformis_GB_05"= "5.575x10³", "E_falciformis_GB_06"= "5.58x10²",
                                                        "E_falciformis_GB_07"= "5.6x10¹", "E_falciformis_GB_08"= "0", 
                                                        "E_falciformis_PB_01"= "1.115x10⁶", "E_falciformis_PB_02"= "1.115x10⁵",
                                                        "E_falciformis_PB_03"= "1.115x10⁴", "E_falciformis_PB_04"= "1.115x10³",
                                                        "E_falciformis_PB_05"= "5.575x10³", "E_falciformis_PB_06"= "5.58x10²",
                                                        "E_falciformis_PB_07"= "5.6x10¹", "E_falciformis_PB_08"= "0","H2O"= "NTC")) +
    scale_y_continuous(name = "Ct")+ 
    geom_jitter(shape=21, position=position_jitter(0.2), color= "black", aes(size= 25, fill= Sample_type))+
    theme_bw() +
    #facet_grid(cols = vars(Extraction))+
    theme(text = element_text(size=16),legend.position = "none")+
    theme(axis.text.x = element_text(angle=90))+
    stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
    labs(tag = "B)")+
    geom_hline(yintercept = 30, linetype = 2)+
    stat_compare_means(method = "anova",
                       aes(label = paste0(..method.., ",\n","p ",..p.format..)), label.y= 17, label.x = 8)+
    stat_compare_means(label = "p.signif", method = "t.test",ref.group = "H2O", 
                       label.y = c(16.5, 21, 25, 27.5, 26, 27.5, 31, 32.5, 30)) ##Determine that samples with <=10^1 oocysts are basically like NTC
  
  #pdf(file = "~/AA_Microbiome/qPCR_Eimeria_Quantification/figures/ABI_Cycler/Fig_2_Mock.pdf", width = 10, height = 15)
  #grid.arrange(C, D)
  #dev.off()
}

###### Infection experiment samples ########
## Considering the standard curve generated with the data from the BioRad Cycler
#Ct = 38 -4(log10Number of Oocysts)
# Number of oocysts = 10^((Ct-38)/-4)

##Estimate number of Oocysts with qPCR Ct value 

##Define real positive based on Tm 
data.inf %>% 
  dplyr::mutate(Infection = case_when(is.na(Tm)  ~ "Negative",
                                      Tm >= 80   ~ "Negative", Tm < 80 ~ "Positive"))%>%
  dplyr::mutate(Qty= 10^((Ct-38)/-4)) -> data.inf

data.inf %>%
  select(Tm, Qty, labels) %>% # select variables to summarise
  na.omit()%>%
  dplyr::group_by(labels)%>%
  summarise_each(funs(min = min, q25 = quantile(., 0.25), median = median, q75 = quantile(., 0.75), 
                      max = max, mean = mean, sd = sd)) -> Sum.inf

data.inf<- join(data.inf, Sum.inf, by= "labels")
data.inf<- join(data.inf, sample.data, by= "labels")

data.inf%>%
  select(labels, Qty_mean, Tm_mean, Infection)%>%
  filter(!labels%in%c("Pos_Ctrl","Neg_Ctrl","FML"))-> data.inf.exp

write.csv(data.inf.exp, "/SAN/Victors_playground/Eimeria_microbiome/sample_data_qPCR.csv")

###Plot Qty by DPI
data.inf%>%
  filter(!labels%in%c("Pos_Ctrl","Neg_Ctrl","FML"))%>%
  ggplot(aes(dpi, Qty))+
  geom_boxplot()+
  xlab("Day post infection")+
  scale_y_log10(name = "log10 Number of Eimeria Oocysts (qPCR)")+ 
  geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "black")+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16)) -> f

data.inf%>%
  filter(!labels%in%c("Pos_Ctrl","Neg_Ctrl","FML"))%>%
  ggplot(aes(OPG, Qty))+
  geom_smooth(method = lm)+
  scale_x_log10(name = "log10 Oocyst per gram feces (Flotation)")+
  scale_y_log10(name = "log10 Number of Eimeria Oocysts (qPCR)")+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  #stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  #facet_wrap(~dpi) +
  stat_cor(label.x = 5.5, label.y = 1, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 5.5, label.y = 1.5) +
  stat_cor(label.x = 5.5,  label.y = 0.5,method = "spearman")-> G

#summary(lm(formula = log10(OPG)~log10(Qty), data = data.inf))

pdf(file = "~/AA_Microbiome/qPCR_Eimeria_Quantification/figures/Track_infection.pdf", width = 10, height = 15)
grid.arrange(f, G)
dev.off()
