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

if(!exists("PS.18S")){
  source("3_MA_18S.R")
}else{
  ##############Load data######################
  PS.18S<- readRDS(file="/SAN/Victors_playground/Eimeria_microbiome/PS_18S.Rds")  
}

############# General overview ##################
summarize_phyloseq(PS.18S) 
rank_names(PS.18S)
tax_table(PS.18S)

############# Filtering ####################### 
###Step 1) Eliminate samples with 0 read counts
PS.18S <- prune_samples(sample_sums(PS.18S)>0, PS.18S)

###Step 2) Eliminate reads Unassigned at Superkingdom level
###Check taxa
table(tax_table(PS.18S)[, "superkingdom"], exclude = NULL)
### Really low proportion just 1 ASV not assigned as Eukaryote
#kable(table(tax_table(PS)[, "superkingdom"], exclude = NULL), format = "html") 
PS.18S <- subset_taxa(PS.18S, !is.na(superkingdom) & !superkingdom %in% c("", "uncharacterized"))

###Step 3) Eliminate reads from host DNA 
##Check for host DNA reads
## Host read numbers
hist(rowSums(otu_table(subset_taxa(PS.18S, genus%in%"Mus"))))
###NO Mus DNA (Super good!) 

###Step 4) Taxonomic filter (eliminate spurius taxa)
##Taxonomic filtering
# Create table, number of features for each phyla
table(tax_table(PS.18S)[, "phylum"], exclude = NULL)

# Compute prevalence of each feature, store as data.frame
PrevAll <- apply(X = otu_table(PS.18S),
                 MARGIN = ifelse(taxa_are_rows(PS.18S), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
PrevAll <- data.frame(Prevalence = PrevAll,
                      TotalAbundance = taxa_sums(PS.18S),
                      tax_table(PS.18S))

plyr::ddply(PrevAll, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Define phyla to filter (present in less than 5 samples)
filterPhyla <- c("Basidiomycota", "Mucoromycota") ###For the full dataset

# Filter entries with unidentified Phylum.
PS.18S <- subset_taxa(PS.18S, !phylum %in% filterPhyla)
rm(PrevAll, filterPhyla)

###Check how many Reads have every superkingdom
###Read counts 
counts_18S <- data.frame(rowSums(otu_table(PS.18S)))
counts_18S[,2] <- rownames(counts_18S)
###Eukaryotic counts
counts_18S[,3] <- as.data.frame(rowSums(otu_table(subset_taxa(PS.18S, superkingdom%in%"Eukaryota"))))
counts_18S[,4] <- as.data.frame(rowSums(otu_table(subset_taxa(PS.18S, genus%in%"Eimeria"))))
colnames(counts_18S) <- c("Read_counts_18S", "labels", "Eukaryota_18S_reads", "ReadsEim18S")
rownames(counts_18S) <- c(1:nrow(counts_18S))
counts_18S <- data.frame(labels = counts_18S$labels, 
                           Read_counts_18S = counts_18S$Read_counts_18S,
                           Eukaryota_18S_reads= counts_18S$Eukaryota_18S_reads,
                           ReadsEim18S= counts_18S$ReadsEim18S)

hist(counts_18S$Read_counts_18S)
summary(counts_18S$Read_counts_18S)
sum(counts_18S$Read_counts_18S)
sum(counts_18S$Eukaryota_18S_reads)
sum(counts_18S$ReadsEim18S)

## Eukaryote reads proportion 
sum(otu_table(subset_taxa(PS.18S, superkingdom %in% "Eukaryota")))/sum(otu_table(PS.18S))
## Eimeria reads proportion 
sum(otu_table(subset_taxa(PS.18S, genus %in% "Eimeria")))/sum(otu_table(PS.18S))
## ASVs by superkingdom
table(tax_table(PS.18S)[, "superkingdom"], exclude = NULL)
## ASVs by genus
table(tax_table(PS.18S)[, "genus"], exclude = NULL) ## Eimeria has 13 ASVs

####Taxa detected
as.data.frame(table(tax_table(PS.18S)[, "phylum"]))
as.data.frame(table(tax_table(PS.18S)[, "genus"]))

counts_18S%>%
  dplyr::select(labels, Read_counts_18S,ReadsEim18S)%>%
  dplyr::mutate(Eimeria_abundance_18S= ReadsEim18S/Read_counts_18S)-> sdt18SEim

#write.csv(sdt18SEim, file = "/SAN/Victors_playground/Eimeria_microbiome/sdt18SEim.csv", row.names = FALSE)