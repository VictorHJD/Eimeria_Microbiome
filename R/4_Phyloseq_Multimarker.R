### Code to 
## 1) Filter process of reads and 
## 2) Generate a working dataframe with sample information + total and Eimeria read counts from multimarker amplification
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
## Eimeria reads proportion 
sum(otu_table(subset_taxa(PS, genus %in% "Eimeria")))/sum(otu_table(PS))
## ASVs by superkingdom
table(tax_table(PS)[, "superkingdom"], exclude = NULL)
## ASVs by genus
table(tax_table(PS)[, "genus"], exclude = NULL) ## Eimeria has 398 ASVs

####Taxa detected
as.data.frame(table(tax_table(PS)[, "phylum"]))
as.data.frame(table(tax_table(PS)[, "genus"]))

## Merge ASVs that have the same taxonomy at a certain taxonomic rank
#PS.Gen <-  tax_glom(PSHigh, "genus", NArm = F)

#PS.Fam<-  tax_glom(PSHigh, "family", NArm = F)

#PS.Ord <-  tax_glom(PSHigh, "order", NArm = F)

#PS.Phy <-  tax_glom(PSHigh, "phylum", NArm = TRUE)

###Summarize total sequencing depths 
sample_data(PS)$dpi <- as.factor(sample_data(PS)$dpi)

sdt <- data.table(as(sample_data(PS), "data.frame"),
                  TotalReads= sample_sums(PS), keep.rownames = T)

setnames(sdt, "rn", "Sample_ID")

###Summarize Eimeria sequencing depths
PS.eimeria <- subset_taxa(PS, genus%in%"Eimeria")

sdtEim <- data.table(as(sample_data(PS.eimeria), "data.frame"),
                     ReadsEim= sample_sums(PS.eimeria), keep.rownames = T)

sdtEim <- dplyr::select(sdtEim, 5,57)

### Join Eimeria read counts with rest of the informaion
sdt <- plyr::join(sdt, sdtEim, by= "labels")

### Add a new variable that will contain the proportion of Eimeria reads by the total per sample
sdt %>%
  mutate(Eimeria_abundance = sdt$ReadsEim/sdt$TotalReads) -> sdt 

#write.csv(sdt, file = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/sdt.csv", row.names = FALSE)
rm(PS.eimeria, PS.genus, sdtEim)