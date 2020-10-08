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
                "Oocyst_1", "Oocyst_2", "Oocyst_3", "Oocyst_4", "Oocyst_5", "Oocyst_6", "Oocyst_7", "Oocyst_8", "Dilution_factor", "Volume", "Sporulated")
  fac.vars <- c("Well", "Sample.Name", "Detector", "Task", "Date", "Operator", "Cycler", "Parasite", "Sample_type", "Extraction", "Strain")  
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
# comput simple linear models from standards
# "Genome copies modeled by Ct"
data.std.lm<- subset(data.std, Task== "Standard") ## Select just data from standards 
data.std.lm %>% 
  select(Sample.Name, Task, Ct, Cycler, Oocyst_count, Parasite, Genome_copies)-> data.std.lm ## Select useful data

##Ct modeled by Oocyst counts; data from different Cyclers
data.std.lm%>%
  ggplot(aes(x = Oocyst_count, y = Ct, color= Cycler)) +
  geom_smooth(method = "lm", se = T) +
  guides(color = FALSE, size = FALSE) +  # Size legend also removed
  scale_x_log10("log 10 Eimeria Oocysts Count", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 20, fill= Cycler), color= "black", alpha= 0.5)+
  stat_cor(label.x = 5, label.y = c(35,30,25),aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+# Add correlation coefficient
  stat_regline_equation(label.x = 5, label.y = c(36.5,31.5,26.5))+
  #stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color= "black")+
  labs(tag = "A)")+
  theme_bw() +
  theme(text = element_text(size=20))+
  annotation_logticks(sides = "b")-> A1

###Model by cycler 
lm.CtABI1 <- lm(Ct~log10(Oocyst_count), subset(data.std.lm, Cycler=="ABI"))
lm.CtABI2 <- lm(Ct~log10(Oocyst_count)+Parasite, subset(data.std.lm, Cycler=="ABI"))
lm.CtABI3 <- lm(Ct~log10(Oocyst_count), subset(data.std.lm, Cycler=="ABI"&Parasite=="E_falciformis"))
lm.CtABI4 <- lm(Ct~log10(Oocyst_count), subset(data.std.lm, Cycler=="ABI"&Parasite=="E_ferrisi"))
lm.CtEpp <- lm(Ct~log10(Oocyst_count), subset(data.std.lm, Cycler=="Eppendorf"))
lm.CtBR <- lm(Ct~log10(Oocyst_count), subset(data.std.lm, Cycler=="BioRad"))

require("rcompanion")
compareLM(lm.CtABI1, lm.CtEpp, lm.CtBR)
compareLM(lm.CtABI1, lm.CtABI2, lm.CtABI3, lm.CtABI4)

lm.CtAll<- lm(Ct~log10(Oocyst_count)+Parasite+Cycler, data.std.lm)
lm.CtPar<- lm(Ct~log10(Oocyst_count)+Parasite, data.std.lm)
lm.CtCyc<- lm(Ct~log10(Oocyst_count)+Parasite+Cycler, data.std.lm)
lm.CtInt<- lm(Ct~log10(Oocyst_count)+Parasite*Cycler, data.std.lm)

compareLM(lm.CtAll, lm.CtPar, lm.CtCyc, lm.CtInt)
summary(lm.CtAll)
summary(lm.CtInt)

#require("ggeffects")
#ggpredict(lm.CtEpp)

##Linear model Genome copies modeled by Oocyst count 
data.std.lm%>%
  ggplot(aes(x = Oocyst_count, y = Genome_copies)) +
  geom_smooth(method = "lm", se = F, color= "black") +
  guides(color = FALSE, size = FALSE) +  # Size legend also removed
  scale_x_log10("log 10 Eimeria Oocysts Count", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10("log 10 Eimeria genome copies/µL gDNA", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 20, fill= Cycler), color= "black", alpha= 0.5)+
    stat_cor(label.x = 4, label.y = c(3,2,1),aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+# Add correlation coefficient
  stat_regline_equation(label.x = 4, label.y = c(3.5,2.5,1.5))+
  labs(tag = "B)")+
  theme_bw() +
  theme(text = element_text(size=20))+
  annotation_logticks(sides = "bl")-> A2

lm.GC <- lm(log10(Genome_copies)~log10(Oocyst_count), data.std.lm)

##Linear model (Standard curve for the rest of experiments)
data.std.lm%>%
  ggplot(aes(x = Ct, y = Genome_copies)) +
  geom_smooth(method = "lm", se = T, color="black") +
  guides(color = FALSE, size = FALSE) +  # Size legend also removed
  scale_y_log10("log 10 Eimeria genome copies/µL gDNA", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 20, fill= Cycler), color= "black", alpha= 0.5)+
  stat_cor(label.x = 5, label.y = 4,aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+# Add correlation coefficient
  stat_regline_equation(label.x = 5, label.y = 4.5)+
  labs(tag = "A)")+
  theme_bw() +
  theme(text = element_text(size=20))+
  annotation_logticks(sides = "l")-> A3

lm.GC1 <- lm(log10(Genome_copies)~Ct, data.std.lm)
summary(lm.GC1)
data.std.lm$predicted<- 10^predict(lm.GC1)
data.std.lm$residuals<- 10^residuals(lm.GC1)

# check 
data.std.lm %>% 
  select(Genome_copies, predicted, residuals) %>%
  head()

#Linear model perfect fit Genome copies ~ Oocyst count
lm.GC2<- lm(Genome_copies~ Oocyst_count, data = data.std.lm)
summary(lm.GC2)

ggplot(data.std.lm, aes(x = Oocyst_count, y = Ct, color= Parasite)) +
  geom_smooth(method = "lm", se = T) +
  guides(color = FALSE, size = FALSE) +  # Size legend also removed
  #geom_point(aes(y = predicted), shape = 21) +
  scale_x_log10("log 10 Eimeria Oocysts Count", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 20, fill= Parasite), color= "black", alpha= 0.5)+
  #stat_cor(label.x = 4,  label.y = 3,method = "spearman")+
  stat_cor(label.x = 5, label.y = c(35,30,25),aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+# Add correlation coefficient
  stat_regline_equation(label.x = 5, label.y = c(36.5,31.5,26.5))+
  #stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  labs(tag = "B)")+
  theme_bw() +
  theme(text = element_text(size=20))+
  annotation_logticks(sides = "b")-> B1

ggplot(data.std.lm, aes(x = Oocyst_count, y = Genome_copies)) +
  geom_smooth(method = "lm", se = FALSE, color= "darkgray") +
  guides(color = FALSE, size = FALSE) +  # Size legend also removed
  #geom_point(aes(y = predicted), shape = 21) +
  scale_x_log10("log 10 Oocysts count", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10("log 10 Eimeria genome copies/µL gDNA", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Cycler), color= "black", alpha= 0.5)+
  #stat_cor(label.x = 4,  label.y = 3,method = "spearman")+
  stat_cor(label.x = 4, label.y = 3.5,aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+# Add correlation coefficient
  stat_regline_equation(label.x = 4, label.y = 3)+
  #stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  labs(tag = "A)")+
  theme_bw() +
  theme(text = element_text(size=20))+
  annotation_logticks(sides = "b")->B2

##Oocysts
data.std%>%
  dplyr::select(Sample.Name,Task,Std_series,Ct,Cycler,Oocyst_count,Parasite,Tm,Date)%>%
  filter(Task=="Standard" & Cycler=="ABI" & Std_series%in%c("A","B"))%>%
  dplyr::group_by(Parasite)%>%
  ggplot(aes(Oocyst_count, Ct))+
  scale_x_log10("log 10 Eimeria Oocysts", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Cycler), color= "black", alpha= 0.5)+
  geom_smooth(color= "black", method = "lm")+            
  stat_cor(label.x = 5,  label.y = 34,method = "spearman")+
  stat_cor(label.x = 5, label.y = 32,aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+        # Add correlation coefficient
  stat_regline_equation(label.x = 5, label.y = 36)+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  labs(tag = "A)")+
  theme_bw() +
  #facet_grid(cols = vars(Std_series))+
  theme(text = element_text(size=20), legend.position = "none")+
  annotation_logticks(sides = "b")+
  coord_cartesian(ylim = c(10, 40)) -> A

std_abi<- lm(formula = log10(Oocyst_count)~Ct, data = subset(data.std, Cycler=="ABI" & Task== "Standard" & Std_series%in%c("A","B")))
summary(std_abi)

##Final standard curve for Inf. experiment samples 
data.std%>%
  dplyr::select(Sample.Name,Task,Std_series,Ct,Cycler,Oocyst_count,Parasite,Tm)%>%
  filter(Task=="Standard" & Parasite=="E_falciformis" & Cycler=="BioRad")%>%
  dplyr::group_by(Sample.Name)%>%
  ggplot(aes(Oocyst_count*8, Ct))+
  scale_x_log10("log 10 Eimeria genome copies/µL gDNA", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Cycler), color= "black", alpha= 0.5)+
  geom_smooth(color= "black",method = "lm")+            
  stat_cor(label.x = 5,  label.y = 34, method = "spearman")+
  stat_cor(label.x = 5, label.y = 32, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+        # Add correlation coefficient
  stat_regline_equation(label.x = 5, label.y = 36)+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  labs(tag = "B)")+
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none")+
  annotation_logticks(sides = "b")+
  coord_cartesian(ylim = c(10, 40))-> B

##Oocysts
data.std%>%
  dplyr::select(Sample.Name,Task,Std_series,Ct,Cycler,Oocyst_count,Parasite,Tm)%>%
  filter(Task=="Standard" & Parasite=="E_falciformis" & Cycler=="BioRad")%>%
  dplyr::group_by(Sample.Name)%>%
  ggplot(aes(Oocyst_count, Ct))+
  scale_x_log10("log 10 Eimeria Oocysts", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Cycler), color= "black", alpha= 0.5)+
  geom_smooth(color= "black",method = "lm")+            
  stat_cor(label.x = 5,  label.y = 34, method = "spearman")+
  stat_cor(label.x = 5, label.y = 32, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+        # Add correlation coefficient
  stat_regline_equation(label.x = 5, label.y = 36)+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  labs(tag = "B)")+
  theme_bw() +
  theme(text = element_text(size=20), legend.position = "none")+
  annotation_logticks(sides = "b")+
  coord_cartesian(ylim = c(10, 40))-> B

std_br<- lm(formula = log10(Oocyst_count)~Ct, data = subset(data.std, Task=="Standard" & Parasite=="E_falciformis" & Cycler=="BioRad"))
summary(std_br)

### Figure 1 Final Standard curves 
pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_1.pdf", width = 8, height = 10)
grid.arrange(A, B)
dev.off()

pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_1.1.pdf", width = 10, height = 10)
grid.arrange(A, B, C)
dev.off()
rm(A,B,C)

pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_1.2.pdf", width = 10, height = 8)
grid.arrange(A1)
dev.off()

pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_1.3.pdf", width = 10, height = 8)
grid.arrange(B1)
dev.off()

pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_1.4.pdf", width = 10, height = 8)
grid.arrange(A2)
dev.off()

pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_1.5.pdf", width = 10, height = 8)
grid.arrange(A3)
dev.off()

## First test on Eppendorf cycler
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
  theme(text = element_text(size=20), legend.position = "none")

summary(lm(formula = log10(Oocyst_count)~Ct, data = subset(data.std, Cycler=="Eppendorf")))

### Comparison between series in BioRad cycler
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
  theme(text = element_text(size=20), legend.position = "none")

summary(lm(formula = log10(Oocyst_count)~Ct, data = subset(data.std, Cycler=="BioRad" & Task== "Standard")))

###### Intersample variation experiment #####
if(Unknowns){
    
data.unk%>%
  dplyr::select(Sample.Name, Task, Ct,Qty,Cycler,Parasite, Sample_type, Extraction, Tm, Oocyst_1, Oocyst_2, Oocyst_3, Oocyst_4, Oocyst_5,Oocyst_6, Oocyst_7, Oocyst_8, Sporulated, Dilution_factor, Volume, Strain)%>%
  dplyr::mutate(Qty= 10^((Ct-36)/-3.1))%>% ## Considering ABI std curve (Fig.1A)
  filter(Sample_type=="Oocysts" & Task=="Unknown" & Cycler== "ABI", Extraction== "Ceramic_beads")%>%
  dplyr::group_by(Sample.Name)%>%
  dplyr::mutate(N= n())%>%
  dplyr::mutate(Total_oocysts= (sum(Oocyst_1, Oocyst_2, Oocyst_3, Oocyst_4, Oocyst_5,Oocyst_6, Oocyst_7, Oocyst_8))/N)%>%
  dplyr::mutate(Oocyst_count= (((Total_oocysts*10000)/8)*Dilution_factor))%>% ##Concentration of Oocyst in the solution
  dplyr::mutate(Sporulated_count= (((Sporulated*10000)/8)*Dilution_factor))%>% ##Concentration of sporulated oocyst in the solution
  dplyr::mutate(Sporulation_rate= (Sporulated_count/Oocyst_count)*100)%>% ##Sporulation rate
  ggplot(aes(Oocyst_count,Qty), geom=c("point", "smooth"))+
  scale_x_log10(name = "log10 Eimeria Oocysts count (Flotation)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10(name = "log 10 Eimeria genome copies/µL gDNA \n (qPCR)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
  geom_jitter(shape=21, position=position_jitter(0.1), aes(size= 25, fill= Task), color= "black", alpha= 0.5)+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  theme_bw() +
  theme(legend.text=element_text(size=20)) +
  theme(legend.key.size = unit(3,"line")) +
  geom_smooth(color= "black", method = "lm")+            
  stat_cor(label.x = 5.5,  label.y = 5,method = "spearman")+
  stat_cor(label.x = 5.5, label.y = 4.75,aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+        # Add correlation coefficient
  stat_regline_equation(label.x = 5.5, label.y = 5.25)+
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(text = element_text(size=20),legend.position = "none")+
  labs(tag = "A)")+
  annotation_logticks(sides = "bl") -> D

data.unk%>%
  dplyr::select(Sample.Name, Task, Ct,Qty,Cycler,Parasite, Sample_type, Extraction, Tm, Oocyst_1, Oocyst_2, Oocyst_3, Oocyst_4, Oocyst_5,Oocyst_6, Oocyst_7, Oocyst_8, Sporulated, Dilution_factor, Volume, Strain)%>%
  dplyr::mutate(Qty= 10^((Ct-36)/-3.1))%>% ## Considering ABI std curve (Fig.1A)
  filter(Sample_type=="Oocysts" & Task=="Unknown")%>%
  dplyr::group_by(Sample.Name)%>%
  dplyr::mutate(N= n())%>%
  dplyr::mutate(Total_oocysts= (sum(Oocyst_1, Oocyst_2, Oocyst_3, Oocyst_4, Oocyst_5,Oocyst_6, Oocyst_7, Oocyst_8))/N)%>%
  dplyr::mutate(Oocyst_count= (((Total_oocysts*10000)/8)*Dilution_factor))%>% ##Concentration of Oocyst in the solution
  dplyr::mutate(Sporulated_count= (((Sporulated*10000)/8)*Dilution_factor))%>% ##Concentration of sporulated oocyst in the solution
  dplyr::mutate(Sporulation_rate= (Sporulated_count/Oocyst_count)*100)%>%
  dplyr::mutate(Sporulation_rate= as.numeric(Sporulation_rate))-> data.unk.lm

data.unk.lm$predicted.Gc<- 10^predict(lm.GC1, data.unk.lm)
data.unk.lm$residuals.Gc<- 10^residuals(lm.GC1, data.unk.lm)

ggplot(data.unk.lm, aes(x = Oocyst_count, y = predicted.Gc)) +
  geom_smooth(method = "lm", se = F,color= "darkgrey") +
  guides(color = FALSE, size = FALSE) +  # Size legend also removed
  #geom_point(aes(y = predicted), shape = 21) +
  scale_x_log10("log 10 Eimeria Oocysts Count", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_log10("log 10 Eimeria genome copies/µL gDNA", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 20, fill= Strain), color= "black", alpha= 0.5)+
  #stat_cor(label.x = 4,  label.y = 3,method = "spearman")+
  stat_cor(label.x = 4, label.y = 7.75,aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+# Add correlation coefficient
  stat_regline_equation(label.x = 4, label.y = 8)+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
  labs(tag = "A)")+
  theme_bw() +
  theme(text = element_text(size=20))+
  annotation_logticks(sides = "bl")->D1

summary(glm(formula = log10(predicted.Gc)~log10(Oocyst_count)+Parasite+Strain+Cycler, data = data.unk.lm))

pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_2.1.pdf", width = 10, height = 8)
grid.arrange(D1)
dev.off()

}
########## Mock samples Experiment #########
if(Mock){
  ##Mean comparison standars against NTC
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
                       aes(label = paste0(..method.., ",\n","p=",..p.format..)), label.y= 33, label.x = 7)+
    stat_compare_means(label = "p.signif", method = "t.test",ref.group = "H2O", 
                       label.y = c(34, 31, 27, 24, 21, 17, 15, 0)) ##Determine that 10^0 meassurments are basically like NTC
  
  ##Compair mock samples qPCR estimation with real oocyst count by two extraction methods
  set.seed(2020)
  data.unk%>%
    select(Sample.Name, Task, Ct,Qty,Cycler,Parasite, Sample_type, Feces_weight, Extraction, Oocyst_count)%>%
    filter(Sample_type=="Feces" & Task=="Unknown")%>%
    dplyr::mutate(Qty= 10^((Ct-36)/-3.1))%>%
    ggplot(aes(x = Oocyst_count, y = Qty, color=Extraction), geom=c("point", "smooth")) +
    scale_x_log10(name = "log10 Eimeria Oocysts (Flotation)", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_y_log10(name = "log10 Eimeria genome copies/µL gDNA", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
    geom_jitter(shape=21, position=position_jitter(0.2), color= "black", aes(size= 25, fill= Extraction), alpha= 0.5)+
    theme_bw() +
    geom_smooth(aes(color= Extraction, fill= Extraction), method = "lm")+
    facet_grid(cols = vars(Extraction))+
    stat_cor(aes(color = Extraction), label.x = log10(100),  label.y = log10(100000),method = "spearman")+
    stat_cor(label.x = log10(100), label.y = log10(50000), aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Extraction))+        # Add correlation coefficient
    stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
    #stat_regline_equation(aes(color = Std_series), label.x = 3, label.y = c(35, 31))
    theme(text = element_text(size=20),legend.position = "none")+
    labs(tag = "A)")+
    annotation_logticks(sides = "bl")
  
  ##Standard curve and ceramic beads data
  data.std%>%
    dplyr::select(Sample.Name,Task,Std_series,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm,Sample_type,Feces_weight,Extraction,Date)%>%
    filter(Task%in%c("Standard", "NTC") & Cycler=="ABI" & Std_series%in%c("A","B"))%>%
    dplyr::mutate(Qty= Qty*8)%>% ##Transform to Genome copies per uL gDNA qPCR 
    dplyr::group_by(Parasite)-> Std.mock
  
  data.unk%>%
    dplyr::select(Sample.Name,Task,Ct,Qty,Cycler,Oocyst_count,Parasite,Tm,Sample_type,Feces_weight,Extraction,Date, NanoDrop, Strain)%>%
    filter(Sample_type=="Feces" & Task=="Unknown" & Extraction!="Glass_beads")%>%
    dplyr::mutate(Qty= 10^((Ct-36)/-3.1), ## Transform Ct to Genome copies per uL gDNA qPCR
                  #GC_ngDNA= Genome_copies/NanoDrop, ## Estimate Genome copies by ng of fecal DNA
                  GC_ngDNA= Qty/50, ## Estimate Genome copies by ng of fecal DNA
                  DNA_sample= NanoDrop*40, ## Estimate total gDNA of sample
                  DNA_g_feces= DNA_sample/Feces_weight, ## Transform it to ng fecal DNA by g of feces
                  GC_gfeces= GC_ngDNA*DNA_g_feces, ## Estimate genome copies by g of feces
                  OPG=Oocyst_count/Feces_weight)-> data.mock ## Estimate oocyst per g of feces for mock samples
  
  data.mock$predicted.Gc<- 10^predict(lm.GC1, data.mock)
  data.mock$residuals.Gc<- 10^residuals(lm.GC1, data.mock) 
  
  data.mock%>%
  bind_rows(data.std.lm)-> data.mock
  
  #rm(Std.mock)
  
  ##
  set.seed(2020)
  data.mock%>%
    dplyr::select(Sample.Name,Task,Qty,Ct,Oocyst_count, predicted.Gc)%>%  
    filter(Task%in%c("Standard", "Unknown"))%>%
    ggplot(aes(x = Oocyst_count, y = predicted.Gc), geom=c("point", "smooth")) +
    scale_x_log10(name = "log10 Eimeria Oocyst count (Flotation)", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_y_log10(name = "log10 Eimeria genome copies/µL gDNA \n (qPCR)", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
    geom_jitter(shape=21, position=position_jitter(0.2), aes(fill= Task), size= 5, color= "black", alpha= 0.5)+
    stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
    theme_bw() +
    theme(legend.text=element_text(size=20)) +
    theme(legend.key.size = unit(3,"line")) +
    geom_smooth(aes(color= Task, fill= Task), method = "lm")+            
    #stat_cor(aes(color = Task), label.x = 2,  label.y = c(7, 6),method = "spearman")+
    stat_cor(label.x = 0.75, label.y =  5.5,aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"), color = Task))+        # Add correlation coefficient
    stat_regline_equation(aes(color = Task), label.x = 0.75, label.y = 6)+
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme(text = element_text(size=20),legend.position = "none")+
    labs(tag = "B)")+
    annotation_logticks(sides = "bl")->E#+
    #theme(legend.position = c(0.85, 0.25), legend.direction = "vertical",
          # Change legend key size and key width
          #legend.key.size = unit(0.25, "cm"),
          #legend.key.width = unit(0.15,"cm"))
  
  pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_2.2.pdf", width = 10, height = 8)
  grid.arrange(E)
  dev.off()
  
  summary(lm(formula = log10(Qty)~log10(Oocyst_count), data = subset(data.mock, Task== "Standard")))
  modelstd<- lm(formula = log10(Qty)~log10(Oocyst_count), data = subset(data.mock, Task== "Standard"))
  summary(lm(formula = log10(Qty)~log10(Oocyst_count), data = subset(data.mock, Task== "Unknown" & Oocyst_count >0)))
  modelmock<- lm(formula = log10(Qty)~log10(Oocyst_count), data = subset(data.mock, Task== "Unknown" & Oocyst_count >0))
  
  data.mock%>%
    dplyr::select(Sample.Name, Qty, Oocyst_count, Task)%>%
    filter(Task== "Unknown"& Oocyst_count >0)%>%
    dplyr::mutate(Qty_estimated= 10^(0.9+log10(Oocyst_count)), Percent_error= ((Qty_estimated- Qty)/Qty_estimated)*100)->Error_mock
  
  mean_ci(Error_mock$Percent_error)
  mean_sd(Error_mock$Percent_error)
  
  data.mock%>%
    select(Sample.Name,Task,Qty,Ct,Oocyst_count, GC_gfeces, OPG)%>%  
    filter(Task== "Unknown")%>%
    ggplot(aes(x = OPG, y = Qty), geom=c("point", "smooth")) +
    scale_x_log10(name = "log10 Eimeria Oocysts per gram of feces (Flotation)", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_y_log10(name = "log10 Eimeria genome copies/µL gDNA \n (qPCR)", 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+ 
    geom_jitter(shape=21, position=position_jitter(0.2), aes(size= 25, fill= Task), color= "black", alpha= 0.5)+
    stat_summary(fun.data=mean_cl_boot, geom="pointrange", shape=16, size=0.5, color="black")+
    theme_bw() +
    theme(legend.text=element_text(size=20)) +
    theme(legend.key.size = unit(3,"line")) +
    geom_smooth(color= "black", method = "lm")+            
    stat_cor(label.x = 2,  label.y = 6,method = "spearman")+
    stat_cor(label.x = 2, label.y = 5.75,aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+        # Add correlation coefficient
    stat_regline_equation(label.x = 2, label.y = 6.25)+
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme(text = element_text(size=20),legend.position = "none")+
    labs(tag = "C)")+
    annotation_logticks(sides = "bl")-> G
  
  pdf(file = "~/AA_Microbiome/Figures/Oocysts_qPCR_Manuscript/Figure_2.pdf", width = 20, height = 15)
  grid.arrange(D,E,G, widths = c(1, 1), layout_matrix = rbind(c(1, 2), c(3, 3)))
  dev.off()
  rm(data.mock, data.unk, data.std, D, E, G)
}

###### Infection experiment samples ########
if(Infexp){
## Considering the standard curve generated with the data from the BioRad Cycler
#Ct = 42x -4(log10Number of genome number per uL gDNA) Figure 1.1B
# Number of genome copies = 10^((Ct-42)/-4)

##Estimate number of genome copies with qPCR Ct value 

##Define real positive and negatives based on Tm 
data.inf %>% 
  dplyr::mutate(Infection = case_when(is.na(Tm)  ~ "Negative",
                                      Tm >= 80   ~ "Negative", Tm < 80 ~ "Positive"))-> data.inf#%>%
  #dplyr::mutate(Qty= 10^((Ct-42)/-4), Genome_copies= 10^((Ct-42)/-4)) -> data.inf

  data.inf$Genome_copies<- 10^predict(lm.GC1, data.inf)
  data.inf$residuals<- 10^residuals(lm.GC1, data.inf)

  data.inf %>%
  select(Tm, Genome_copies,labels) %>% # select variables to summarise
  na.omit()%>%
  dplyr::group_by(labels)%>%
  dplyr::summarise_each(funs(min = min, q25 = quantile(., 0.25), median = median, q75 = quantile(., 0.75), 
                      max = max, mean = mean, sd = sd)) -> Sum.inf

data.inf<- join(data.inf, Sum.inf, by= "labels")

data.inf%>%
  select(labels, Genome_copies_mean, Tm_mean, Infection)%>%
  filter(!labels%in%c("Pos_Ctrl","Neg_Ctrl","FML"))%>% ## Replace NAs in real negative samples to 0 
  dplyr::mutate(Genome_copies_mean= replace_na(Genome_copies_mean, 0))-> data.inf.exp

#write.csv(data.inf.exp, "/SAN/Victors_playground/Eimeria_microbiome/qPCR/sample_data_qPCR.csv", row.names = FALSE)
}

rm(data.inf, data.std, Sum.inf)
