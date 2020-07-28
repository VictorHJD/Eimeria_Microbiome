###Alpha diversity graphics
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
library("microbiome")
library("grid")
library("ggsci")
library("knitr")
library("kableExtra")

if(!exists("PS")){
  source("2_MA_Multimarker.R")
}else{
  ##############Load data######################
  PS<- readRDS(file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/PhyloSeqData_All.Rds")  
}

############# General overview ##################
summarize_phyloseq(PS) 
rank_names(PS)
tax_table(PS)

############# Filtering ####################### 
###Step 1) Eliminate samples with 0 read counts
PS <- prune_samples(sample_sums(PS)>0, PS)

###Step 2) Eliminate reads Unassigned at Superkingdom level
###Check taxa
table(tax_table(PS)[, "superkingdom"], exclude = NULL)
### Relative low proportion just 269 ASV's are not assigned as Eukaryotes or Bacteria
#kable(table(tax_table(PS)[, "superkingdom"], exclude = NULL), format = "html") 
PS <- subset_taxa(PS, !is.na(superkingdom) & !superkingdom %in% c("", "uncharacterized"))

###Step 3) Eliminate reads from host DNA 
##Check for host DNA reads
## Host read numbers
hist(rowSums(otu_table(subset_taxa(PS, genus%in%"Mus"))))
sum(otu_table(subset_taxa(PS, genus%in%"Mus")))/sum(otu_table(PS))
###Just 1.22% of Mus DNA 
###Eliminate reads assigned as "Mus"
#PS <- subset_taxa(PS, genus!= "Mus") ##Taxa en OTU table decreased dramatically, check and then apply 

###Step 4) Taxonomic filter (eliminate spurius taxa)
##Taxonomic filtering
# Create table, number of features for each phyla
table(tax_table(PS)[, "phylum"], exclude = NULL)

# Compute prevalence of each feature, store as data.frame
PrevAll <- apply(X = otu_table(PS),
                 MARGIN = ifelse(taxa_are_rows(PS), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
PrevAll <- data.frame(Prevalence = PrevAll,
                      TotalAbundance = taxa_sums(PS),
                      tax_table(PS))

plyr::ddply(PrevAll, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Define phyla to filter (present in less than 5 samples)
filterPhyla <- c("Tardigrada", "Microsporidia", "Fusobacteria", "Equinodermata", "Cryptomycota", "Candidatus Saccharibacteria") ###For the full dataset

# Filter entries with unidentified Phylum.
PS <- subset_taxa(PS, !phylum %in% filterPhyla)
rm(PrevAll, filterPhyla)

###Check how many Reads have every superkingdom
###Read counts 
counts_multi <- data.frame(rowSums(otu_table(PS)))
counts_multi[,2] <- rownames(counts_multi)
###Bacterial and Eukaryotic counts
counts_multi[,3] <- as.data.frame(rowSums(otu_table(subset_taxa(PS, superkingdom%in%"Bacteria"))))
counts_multi[,4] <- as.data.frame(rowSums(otu_table(subset_taxa(PS, superkingdom%in%"Eukaryota"))))
colnames(counts_multi) <- c("Read_counts", "labels", "Bacteria_reads", "Eukaryota_reads")
rownames(counts_multi) <- c(1:nrow(counts_multi))
counts_multi <- data.frame(labels = counts_multi$labels, 
                             Read_counts = counts_multi$Read_counts,
                             Bacteria_reads= counts_multi$Bacteria_reads,
                             Eukaryota_reads= counts_multi$Eukaryota_reads)

hist(counts_multi$Read_counts)
summary(counts_multi$Read_counts)
sum(counts_multi$Read_counts)
sum(counts_multi$Bacteria_reads)
sum(counts_multi$Eukaryota_reads)

## Bacteria reads proportion
sum(otu_table(subset_taxa(PS, superkingdom %in% "Bacteria")))/sum(otu_table(PS))
## Eukaryote reads proportion 
sum(otu_table(subset_taxa(PS, superkingdom %in% "Eukaryota")))/sum(otu_table(PS))
## ASVs by superkingdom
table(tax_table(PS)[, "superkingdom"], exclude = NULL)

####Taxa detected
as.data.frame(table(tax_table(PS)[, "phylum"]))
as.data.frame(table(tax_table(PS)[, "genus"]))

## Merge ASVs that have the same taxonomy at a certain taxonomic rank (in this case phylum and genus)
PS.Gen <-  tax_glom(PSHigh, "genus", NArm = F)
summarize_phyloseq(PS.Gen)

PS.Fam<-  tax_glom(PSHigh, "family", NArm = F)
summarize_phyloseq(PS.Fam)

PS.Ord <-  tax_glom(PSHigh, "order", NArm = F)
summarize_phyloseq(PS.Ord)

PS.Phy <-  tax_glom(PSHigh, "phylum", NArm = TRUE)
summarize_phyloseq(PS.Phy)

###Seq depths
qplot(log10(rowSums(otu_table(PS))),binwidth=0.2) +
  xlab("Log10 counts-per-sample")+ ylab("Count")+
  theme_bw()+
  labs(tag = "A)")

###Summarize sequencing depths 
sample_data(PS)$dpi <- as.factor(sample_data(PS)$dpi)

sdt <- data.table(as(sample_data(PS), "data.frame"),
                  TotalReads= sample_sums(PS), keep.rownames = T)

setnames(sdt, "rn", "Sample_ID")

##Check for Eimeria
PS.eimeria <- subset_taxa(PS, genus%in%"Eimeria")

###Eimeria 
sdtEim <- data.table(as(sample_data(PS.eimeria), "data.frame"),
                     ReadsEim= sample_sums(PS.eimeria), keep.rownames = T)


sdtEim <- dplyr::select(sdtEim, 5,57)

sdt <- plyr::join(sdt, sdtEim, by= "labels")

sdt %>%
  mutate(Eimeria_abundance = sdt$ReadsEim/sdt$TotalReads) -> sdt ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

####For here one should go into a different code 



####Normalizing data :S 
PSlog <- transform_sample_counts(PS, function(x) log(1 + x))


###Rarefaction curve 
#rarcurv <- vegan::rarecurve(otu_table(PS),
#                            label = F)

## Eliminate samples with low counts
##At 6000 reads samples have reached the species saturation (from rarecurve2!)
PSHigh <- prune_samples(sample_sums(PS)>=6000, PS)
summarize_phyloseq(PSHigh)
nsamples(PSHigh)
ntaxa(PSHigh)
sum(otu_table(PSHigh))
sum(otu_table(subset_taxa(PSHigh, superkingdom %in% "Bacteria")))/sum(otu_table(PSHigh))
sum(otu_table(subset_taxa(PSHigh, superkingdom %in% "Eukaryota")))/sum(otu_table(PSHigh))
#sum(otu_table(subset_taxa(PSHigh, genus%in%"Mus")))/sum(otu_table(PSHigh))

table(tax_table(PSHigh)[, "superkingdom"], exclude = NULL)

#rarcurv2 <- vegan::rarecurve(otu_table(PSHigh),
#                             label = F)
## Merge ASVs that have the same taxonomy at a certain taxonomic rank (in this case phylum and genus)
PS.Gen <-  tax_glom(PSHigh, "genus", NArm = F)
summarize_phyloseq(PS.Gen)

PS.Fam<-  tax_glom(PSHigh, "family", NArm = F)
summarize_phyloseq(PS.Fam)

PS.Ord <-  tax_glom(PSHigh, "order", NArm = F)
summarize_phyloseq(PS.Ord)

PS.Phy <-  tax_glom(PSHigh, "phylum", NArm = TRUE)
summarize_phyloseq(PS.Phy)

plot_bar(PS, fill="phylum") + facet_wrap(~dpi, scales= "free_x", nrow=1)
plot_richness(PS, x="dpi", color="dpi", measures=c("Chao1", "Shannon", "Simpson"))+ 
  geom_boxplot()+
  theme_bw()

alphaDiv <- phyloseq::estimate_richness(PS)

#eveAll <- evenness(PS1, c("pielou", "simpson", "evar", "bulla"))

###Difference accross DPI?
#pairwise.wilcox.test(alphaDiv$Shannon, sample_data(PS1)$dpi)

#Beta diversity 
dis <- phyloseq::distance(PS, method="bray")
#ordi <- ordinate(PS1, method="PCoA", distance="bray")
#beta.plot <- plot_ordination(PS1, ordi, color="dpi")+
#  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")+
#  theme_classic(base_size = 15, base_family = "Helvetica")+
#  geom_point(size=5, alpha= 0.75)+
#  labs(tag= "A)")

##Permanova
#require("vegan")
#adonis(dis ~ sample_data(PS1)$dpi)

###Using Negative Bionomial in microbiome differential abundance testing (not working check it again :( )
#require("DESeq2")
#dipdds <- phyloseq_to_deseq2(PS1.16SA, ~ dpi)
#DESeq(dipdds)

# calculate geometric means prior to estimate size factors
#gm_mean <- function(x, na.rm=TRUE){
#  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
#}
#geoMeans <- apply(counts(dipdds), 1, gm_mean)
#dipdds <- estimateSizeFactors(dipdds, geoMeans)
#dipdds <- DESeq(dipdds, fitType="local")

###Ordination
#OrdPS16SA <- ordinate(PS16SAlog, method = "MDS", distance = "wunifrac")
#evals <- out.wuf.log$values$Eigenvalues
#plot_ordination(pslog, out.wuf.log, color = "age_binned") +
#  labs(col = "Binned Age") +
#  coord_fixed(sqrt(evals[2] / evals[1]))

###Composition
PS1c <- PS1 %>%
  aggregate_taxa(level = "phylum") %>%  
  microbiome::transform(transform = "compositional")

PS1c %>%
  plot_composition(sample.sort = "Eimeria", otu.sort = "abundance", group_by= "dpi") + theme(legend.position="none")

plot_composition(PS1c,
                 average_by = "dpi")+ theme(legend.position="bottom")

pdf(file = "~/AA_Microbiome/Figures/Composition_Infection.pdf", width = 15, height = 15)
grid.arrange(g,h,i, ncol= 1, nrow= 3)
dev.off()


###Incorporate qPCR information
if(!exists("data.inf.exp")){
  data.inf.exp<- read.csv(file="/SAN/Victors_playground/Eimeria_microbiome/sample_data_qPCR.csv")
}

data.inf.exp%>%
  select(labels, Qty_mean, Tm_mean, Infection)%>%
  distinct()-> data.inf.exp

sdt<- join(sdt, data.inf.exp, by="labels")

##Get single marker data (18S)
PS1.18S<- readRDS("/SAN/Victors_playground/Eimeria_microbiome/PS1_18S.Rds")

PS.18S.Eim <- subset_taxa(PS1.18S, genus%in%"Eimeria")
sdt18SEim <- data.table(as(sample_data(PS.18S.Eim), "data.frame"),
                        ReadsEim18S= sample_sums(PS.18S.Eim), keep.rownames = T)


sdt18SEim <- dplyr::select(sdt18SEim, 5,57)

sdt <- plyr::join(sdt, sdt18SEim, by= "labels")

####OPG vs reads Eimeria (Multiamplicon) 
sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  ggplot(aes(OPG, ReadsEim))+
  geom_smooth(method = lm)+
  scale_x_log10(name = "log10 Oocyst per gram feces (Flotation)")+
  scale_y_log10("log10 Sequence reads count Multiamplicon (Eimeria)")+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  #facet_wrap(~dpi) +
  stat_cor(label.x = 5.5, label.y = 1.25, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 5.5, label.y = 1.5)+
  stat_cor(label.x = 5.5,  label.y = 1,method = "spearman")-> opgrem

####OPG vs reads Eimeria (Single amplicon) 
sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  ggplot(aes(OPG, ReadsEim18S))+
  geom_smooth(method = lm)+
  scale_x_log10(name = "log10 Oocyst per gram feces (Flotation)")+
  scale_y_log10("log10 Sequence reads count 18S (Eimeria)")+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  #facet_wrap(~dpi) +
  stat_cor(label.x = 5.5, label.y = 1.25, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 5.5, label.y = 1.5)+
  stat_cor(label.x = 5.5,  label.y = 1,method = "spearman")-> opgres

####OPG vs qPCR 
sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  ggplot(aes(OPG, Qty_mean))+
  geom_smooth(method = lm)+
  scale_x_log10(name = "log10 Oocyst per gram feces (Flotation)")+
  scale_y_log10(name = "log10 Number of Eimeria Oocysts (qPCR)")+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_cor(label.x = 5.5, label.y = 1.5, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 5.5, label.y = 2)+
  stat_cor(label.x = 5.5,  label.y = 1,method = "spearman")-> opgqpcr

sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  ggplot(aes(Qty_mean, OPG))+
  geom_smooth(method = lm)+
  scale_y_log10(name = "log10 Oocyst per gram feces (Flotation)")+
  scale_x_log10(name = "log10 Number of Eimeria Oocysts (qPCR)")+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_cor(label.x = 4.5, label.y = 2.5, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 4.5, label.y = 3)+
  stat_cor(label.x = 4.5,  label.y = 2,method = "spearman")

summary(lm(OPG~Qty_mean*dpi,  data = sdt))
  
####qPCR vs Multiamplicon
sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  ggplot(aes(Qty_mean, ReadsEim))+
  geom_smooth(method = lm)+
  scale_x_log10(name = "log10 Number of Eimeria Oocysts (qPCR)")+
  scale_y_log10(name = "log10 Sequence reads count Multiamplicon (Eimeria)")+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "D)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_cor(label.x = 4.0, label.y = 0.75, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 4.0, label.y = 1)+
  stat_cor(label.x = 4.0,  label.y = 0.5,method = "spearman")-> qpcrrem

####qPCR vs Single amplicon
sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  ggplot(aes(Qty_mean, ReadsEim18S))+
  geom_smooth(method = lm)+
  scale_x_log10(name = "log10 Number of Eimeria Oocysts (qPCR)")+
  scale_y_log10(name = "log10 Sequence reads count 18S (Eimeria)")+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "E)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_cor(label.x = 4.0, label.y = 1.75, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 4.0, label.y = 2)+
  stat_cor(label.x = 4.0,  label.y = 1.5,method = "spearman")-> qpcrres

####Multiamplicon vs Single amplicon
sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  ggplot(aes(ReadsEim18S, ReadsEim))+
  geom_smooth(method = lm)+
  scale_x_log10(name = "log10 Sequence reads count Multiamplicon (Eimeria)")+
  scale_y_log10(name = "log10 Sequence reads count 18S (Eimeria)")+
  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
  labs(tag= "F)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_cor(label.x = 3.0, label.y = 0.75, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 3.0, label.y = 1)+
  stat_cor(label.x = 3.0,  label.y = 0.5,method = "spearman")-> remres

#opgea <- ggplot(sdt, aes(OPG, Eimeria_abundance))+
#  geom_smooth(method = lm)+
# scale_x_log10(name = "log10 Oocyst per gram feces (Flotation)")+
#  scale_y_continuous("Relative abundance of Eimeria (Sequencing reads)")+
#  geom_jitter(shape=21, position=position_jitter(0.2), size=5, aes(fill= dpi), color= "black")+
#  labs(tag= "B)")+
#  theme_bw()+
#  theme(text = element_text(size=16))+
#  stat_cor(label.x = 5.5, label.y = 0.1, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
#  stat_regline_equation(label.x = 5.5, label.y = 0.2)+
#  stat_cor(label.x = 5.5,  label.y = 0.15,method = "spearman")


pdf(file = "~/AA_Microbiome/Figures/OPG_qPCR_Eimeria_Reads_Multimarker.pdf", width = 15, height = 20)
grid.arrange(opgre,opgqpcr, qpcrre, ncol= 1, nrow= 3)
dev.off()

pdf(file = "~/AA_Microbiome/Figures/OPG_Eimeria_Reads_Multimarker.pdf", width = 15, height = 15)
grid.arrange(opgre,opgea, ncol= 1, nrow= 2)
dev.off()

pdf(file = "~/AA_Microbiome/Figures/Eimeria_quant_methods.pdf", width = 15, height = 20)
grid.arrange(opgrem, qpcrrem, opgres, qpcrres, opgqpcr, remres, ncol= 2, nrow= 3)
dev.off()


###Course of infection 
#compare_means(OPG ~ dpi,  data = sdt) #Adjust table to run it
sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  dplyr::arrange(dpi)%>%
  filter(dpi%in%c("4", "8"))->comp ##for comparison later 

sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  ggplot(aes(dpi, OPG))+
  geom_boxplot()+
  xlab("Day post infection")+
  scale_y_log10("log10 Oocyst per gram feces (Flotation)")+
  geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "black")+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))-> a

sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  dplyr::arrange(dpi)%>%
  ggplot(aes(as.numeric(as.character(dpi)), OPG, colour= EH_ID))+
  xlab("Day post infection")+
  scale_y_log10("log10 Oocyst per gram feces (Flotation)")+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10))+
  geom_line(alpha= 0.5)+
  geom_jitter(shape=21, position=position_jitter(0.0), size=2.5, aes(fill= EH_ID), color= "black", alpha= 0.5)+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_smooth(color= "black", method = "loess")-> a2

ggpaired(comp, x= "dpi", y= "OPG",line.color= "gray", line.size= 0.4, color= "dpi")+
  scale_y_log10("log10 Oocyst per gram feces (Flotation)")+
  xlab("Day post infection")+
  scale_color_npg()+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16)) -> a3

sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  ggplot(aes(dpi, ReadsEim))+
  geom_boxplot()+
  xlab("Day post infection")+
  scale_y_log10("log10 Sequence reads count \n Multiamplicon (Eimeria)")+
  geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "black")+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))-> b

sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  dplyr::arrange(dpi)%>%
  ggplot(aes(as.numeric(as.character(dpi)), ReadsEim, colour= EH_ID))+
  xlab("Day post infection")+
  scale_y_log10("log10 Sequence reads count \n Multiamplicon (Eimeria)")+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10))+
  geom_line(alpha=0.5)+
  geom_jitter(shape=21, position=position_jitter(0.0), size=2.5, aes(fill= EH_ID), color= "black", alpha= 0.5)+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_smooth(color= "black", method = "loess")-> b2

ggpaired(comp, x= "dpi", y= "ReadsEim",line.color= "gray", line.size= 0.4, color= "dpi")+
  scale_y_log10("log10 Sequence reads count \n Multiamplicon (Eimeria)")+
  xlab("Day post infection")+
  scale_color_npg()+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))-> b3

sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  ggplot(aes(dpi, Qty_mean))+
  geom_boxplot()+
  xlab("Day post infection")+
  scale_y_log10(name = "log10 Number of Eimeria Oocysts (qPCR)")+ 
  geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "black")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))-> c

sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  dplyr::arrange(dpi)%>%
  ggplot(aes(as.numeric(as.character(dpi)), Qty_mean, colour= EH_ID))+
  xlab("Day post infection")+
  scale_y_log10("log10 Number of Eimeria Oocysts (qPCR)")+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10))+
  geom_line(alpha=0.5)+
  geom_jitter(shape=21, position=position_jitter(0.0), size=2.5, aes(fill= EH_ID), color= "black", alpha= 0.5)+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_smooth(color= "black", method = "loess")-> c2

ggpaired(comp, x= "dpi", y= "Qty_mean",line.color= "gray", line.size= 0.4, color= "dpi")+
  scale_y_log10("log10 Number of Eimeria Oocysts (qPCR)")+
  xlab("Day post infection")+
  scale_color_npg()+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))-> c3

sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  ggplot(aes(dpi, ReadsEim18S))+
  geom_boxplot()+
  xlab("Day post infection")+
  scale_y_log10("log10 Sequence reads count 18S (Eimeria)")+
  geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "black")+
  labs(tag= "D)")+
  theme_bw()+
  theme(text = element_text(size=16))-> d

sdt%>%
  distinct(Sample_ID, .keep_all=T)%>%
  dplyr::arrange(dpi)%>%
  ggplot(aes(as.numeric(as.character(dpi)), ReadsEim18S, colour= EH_ID))+
  xlab("Day post infection")+
  scale_y_log10("log10 Sequence reads count 18S (Eimeria)")+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10))+
  geom_line(alpha= 0.5)+
  geom_jitter(shape=21, position=position_jitter(0.0), size=2.5, aes(fill= EH_ID), color= "black", alpha= 0.5)+
  labs(tag= "D)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_smooth(color= "black", method = "loess")-> d2

ggpaired(comp, x= "dpi", y= "ReadsEim18S",line.color= "gray", line.size= 0.4, color= "dpi")+
  scale_y_log10("log10 Sequence reads count 18S (Eimeria)")+
  xlab("Day post infection")+
  scale_color_npg()+
  labs(tag= "D)")+
  theme_bw()+
  theme(text = element_text(size=16))-> d3

pdf(file = "~/AA_Microbiome/Figures/Course_of_Eimeria_Infection_Multimethods.pdf", width = 15, height = 10)
grid.arrange(a,c,b,d, ncol= 2, nrow= 2)
dev.off()

pdf(file = "~/AA_Microbiome/Figures/Course_of_Eimeria_Infection_Multimethods_sample.pdf", width = 15, height = 10)
grid.arrange(a2,c2,b2,d2, ncol= 2, nrow= 2)
dev.off()

pdf(file = "~/AA_Microbiome/Figures/Course_of_Eimeria_Infection_dpi_comparison.pdf", width = 15, height = 10)
grid.arrange(a3,c3,b3,d3, ncol= 2, nrow= 2)
dev.off()
##For time series graph (under construction)
require("reshape")

sdt_reshape<- melt(sdt, id.vars = "EH_ID")
sdt_reshape%>%
  select(EH_ID)%>%
  plyr::join(sdt, by = "EH_ID")%>%
  group_by(EH_ID)%>%
  ggplot(aes(x= dpi, y= OPG))+
  geom_line()

library("TSrepr")
  
sdt%>%
  dplyr::group_by(EH_ID)%>%
  dplyr::arrange(as.factor(as.character(dpi)))%>%
  select(dpi,OPG,Qty_mean, ReadsEim, ReadsEim18S)%>%
  filter(EH_ID%in%c("LM0206"))

###Alpha diversity along infection
alphaDiv$labels<- rownames(alphaDiv)
eveAll$labels <- rownames(eveAll)

sdt<- plyr::join(sdt, alphaDiv, by= "labels")
sdt<- plyr::join(sdt, eveAll, by= "labels")

d<- ggplot(sdt, aes(dpi, Chao1, fill= dpi))+
  geom_boxplot()+
  xlab("Day post infection")+
  ylab("Richness (Chao1 Index)")+
  labs(tag= "A)")+
  theme_bw()

e<- ggplot(sdt, aes(dpi, Shannon, fill= dpi))+
  geom_boxplot()+
  xlab("Day post infection")+
  ylab("Diversity (Shannon Index)")+
  labs(tag= "B)")+
  theme_bw()

f<- ggplot(sdt, aes(dpi, pielou, fill= dpi))+
  geom_boxplot()+
  xlab("Day post infection")+
  ylab("Evenness (Pielou's Index)")+
  labs(tag= "C)")+
  theme_bw()

pdf(file = "~/AA_Microbiome/Figures/Alpha_div_Infection.pdf", width = 15, height = 15)
grid.arrange(d,e,f, ncol= 1, nrow= 3)
dev.off()


##Rarefaction curves
rarecurve(t(otu_table(PS)), step=50, cex=0.5)


sampleData(PS)$parental <- getVariable(PS, "HybridStatus") %in% "parental strains"
PSHigh <- prune_samples(sample_sums(PS)>6700, PS) ## Eliminate samples with low counts considering as limit Mean - sd 

###Alpha diversity
alphaDiv <- as.data.frame(estimate_richness(PSHigh, measures=c("Observed", "InvSimpson", "Shannon", "Chao1")))
####Add richness alpha diversity data to phyloseq 
sample_data(PSHigh)$Shannon <- alphaDiv$Shannon
sample_data(PSHigh)$InvSimpson <- alphaDiv$InvSimpson
sample_data(PSHigh)$Chao1 <- alphaDiv$Chao1
PSdata <- meta(PSHigh)


###Compair differences in Shannon idex between groups (e.g. time points)
ggplot(PSdata, aes(x= as.factor(dpi), y= Shannon, fill= as.factor(dpi)))+
  geom_boxplot()+
  geom_jitter()+
  facet_grid(HybridStatus~.)+ 
  theme_classic()

plot_regression(Shannon ~ as.factor(dpi), meta(PSHigh)) + xlab("DPI") + ylab("Diversity (Shannon)")

plot_landscape(PSHigh, "MDS", "bray", col = "dpi")


###Compositional 

PSComp <- microbiome::transform(PSHigh, "compositional")

plot_composition(PSHigh, average_by = "dpi", transform = "compositional")

plot_composition(PSComp, sample.sort = "Apicomplexa", otu.sort = "abundance") +
  # Set custom colors
  scale_fill_manual(values = default_colors("Phylum")[taxa(PSComp1)]) 


##Plot_richness
prune_species(speciesSums(PS)>0, PS)
plot_richness(PS, x= "dpi", color = "Strain" , measures = c("Chao1", "Shannon")) + geom_line()
plot_richness(PSHigh, x= as.character("dpi"), color = "Strain" , measures = c("Chao1", "Shannon"))+ geom_line()
plot_richness(PSHigh, x= "parental", color = "dpi" , measures = c("Chao1", "Shannon"))
##rarefy 
rarecurve(t(otu_table(PSHigh)), step=100, cex=0.5)
PS.rarefied <- rarefy_even_depth(PSHigh, rngseed=1, sample.size=0.9*min(sample_sums(PSHigh)), replace=F)


###Plot ASVs

plot_bar(PS, fill="phylum") + facet_wrap(~dpi, scales= "free_x", nrow=1)
plot_bar(PSHigh, fill="phylum") + facet_wrap(~dpi, scales= "free_x", nrow=1)
plot_bar(PS.rarefied, fill="phylum") + facet_wrap(~dpi, scales= "free_x", nrow=1)


##Beta diversity
# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist <- distance(PSHigh, method="unifrac", weighted=F)
ordination <- ordinate(PS, method="PCoA", distance=wunifrac_dist)
plot_ordination(PS, ordination, color="dpi") + theme(aspect.ratio=1)

###Emanuel's code
nameOtuByTax <- function(PS, taxon="genus"){
  otable <- otu_table(PS)
  rownames(otable) <- make.unique(as.character(tax_table(PS)[, "genus"]))
  otable
}

ordi <- ordinate(PS, method="MDS", distance="bray")
ordi <- ordinate(PS.rarefied, method="MDS", distance="bray")

plot_ordination(PS, ordi, color="dpi")
plot_ordination(PS.rarefied, ordi, color="dpi") ##add circles by dpi 

PSHigh <- prune_samples(sample_sums(PS)>6700, PS) ## Eliminate samples with low counts considering as limit Mean - sd 


ordi2 <- ordinate(PSHigh, method="MDS", distance="bray")

plot_ordination(PSHigh, ordi2, color="OPG", label="dpi")


ordi3 <- ordinate(subset_taxa(PSHigh, !phylum%in%"Apicomplexa"),
                  method="MDS", distance="bray")

plot_ordination(subset_taxa(PSHigh, !phylum%in%"Apicomplexa"),
                ordi3, color="dpi")


ordi4 <- ordinate(PSHigh, formula = ~dpi, "MDS", "bray")
ordi4plot <- plot_ordination(PSHigh, ordi4, "Sample", color="dpi",shape="HybridStatus")


#ordi4 <- ordinate(PSHigh, method="PCoA", distance="unifrac")
#plot_ordination(PSHigh, ordi4, color="dpi") + theme(aspect.ratio=1)



plot_richness(PSHigh, "dpi", color= "OPG", measures="Shannon") + stat_smooth()

plot_richness(PSHigh, "dpi", color= "OPG", measures="Chao1") + stat_smooth()

plot_richness(PSHigh, "dpi", color= "fecweight_DNA",measures="Chao1") + stat_smooth()

sexplot <- plot_richness(PSHigh, x= "Sex") + geom_boxplot(aes(x = Sex, color = dpi), alpha = 0.1)

dpiplot <- plot_richness(PSHigh, x= "dpi") + geom_boxplot(aes(x = dpi, color = dpi), alpha = 0.1)


####

# prune to just the top 100 most abundant OTUs across all samples (crude).
PS50 <- prune_taxa(names(sort(taxa_sums(PS), TRUE))[1:50], PS)
jg <-  make_network(PS100, "taxa", "jaccard", 0.3)
plot_network(jg, PS100, "taxa", color = "phylum", line_weight = 0.4, label = NULL)


TopOTUs <- names(sort(taxa_sums(PS), TRUE)[1:50])
PS10 <- prune_taxa(TopOTUs, PS)
plot_bar(PS10, "dpi", fill = "Strain", facet_grid = ~phylum)
