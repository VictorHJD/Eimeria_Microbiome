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
PS <- subset_taxa(PS, !genus%in%c( "Mus")) ##Taxa en OTU table decreased dramatically, check and then apply 

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
counts_multi[,5] <- as.data.frame(rowSums(otu_table(subset_taxa(PS, genus%in%"Eimeria"))))
counts_multi[,6] <- as.data.frame(rowSums(otu_table(subset_taxa(PS, species%in%"Eimeria ferrisi"))))
counts_multi[,7] <- as.data.frame(rowSums(otu_table(subset_taxa(PS, species%in%"Eimeria falciformis"))))
counts_multi[,8] <- as.data.frame(rowSums(otu_table(subset_taxa(PS, species%in%"Eimeria papillata"))))
counts_multi[,9] <- as.data.frame(rowSums(otu_table(subset_taxa(PS, species%in%"Eimeria telekii"))))
counts_multi[,10] <- as.data.frame(rowSums(otu_table(subset_taxa(PS, species%in%"Eimeria vermiformis"))))
counts_multi[,11] <- as.data.frame(rowSums(otu_table(subset_taxa(PS, species%in%"Eimeria scholtysecki"))))
counts_multi[,12] <- as.data.frame(rowSums(otu_table(subset_taxa(PS, species%in%"Eimeria meleagrimitis"))))
counts_multi[,13] <- as.data.frame(rowSums(otu_table(subset_taxa(PS, species%in%"Eimeria arnyi"))))
colnames(counts_multi) <- c("Read_counts", "labels", "Bacteria_reads", "Eukaryota_reads", 
                            "ReadsEim", "ReadsEimfer", "ReadsEimfal", "ReadsEimpap", "ReadsEimtel",
                            "ReadsEimver", "ReadsEimsch", "ReadsEimmel", "ReadsEimarn")
rownames(counts_multi) <- c(1:nrow(counts_multi))
counts_multi <- data.frame(labels = counts_multi$labels, 
                             Read_counts = counts_multi$Read_counts,
                             Bacteria_reads= counts_multi$Bacteria_reads,
                             Eukaryota_reads= counts_multi$Eukaryota_reads,
                             ReadsEim= counts_multi$ReadsEim,
                             ReadsEimfer= counts_multi$ReadsEimfer,
                           ReadsEimfal= counts_multi$ReadsEimfal,
                           ReadsEimpap= counts_multi$ReadsEimpap,
                           ReadsEimtel= counts_multi$ReadsEimtel,
                           ReadsEimver= counts_multi$ReadsEimver,
                           ReadsEimsch= counts_multi$ReadsEimsch,
                           ReadsEimmel= counts_multi$ReadsEimmel,
                           ReadsEimarn= counts_multi$ReadsEimarn)

hist(counts_multi$Read_counts)
summary(counts_multi$Read_counts)
sum(counts_multi$Read_counts)
sum(counts_multi$Bacteria_reads)
sum(counts_multi$Eukaryota_reads)
sum(counts_multi$ReadsEim)
sum(counts_multi$ReadsEimfer)
sum(counts_multi$ReadsEimfal)
sum(counts_multi$ReadsEimpap)
sum(counts_multi$ReadsEimtel)
sum(counts_multi$ReadsEimver)
sum(counts_multi$ReadsEimsch)
sum(counts_multi$ReadsEimmel)
sum(counts_multi$ReadsEimarn)

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
table(tax_table(PS)[, "species"], exclude = NULL) ## Eimeria ferrisi has 116 ASVs

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
counts_multi%>%
  select(labels, ReadsEim, ReadsEimfer, ReadsEimfal, ReadsEimpap, ReadsEimtel, ReadsEimver, ReadsEimsch, ReadsEimmel, ReadsEimarn)%>%
  distinct(labels, .keep_all = TRUE)-> sdtEim

### Join Eimeria read counts with rest of the informaion
sdt <- plyr::join(sdt, sdtEim, by= "labels")

### Add a new variable that will contain the proportion of Eimeria reads by the total per sample
sdt %>%
  mutate(Eimeria_abundance = sdt$ReadsEim/sdt$TotalReads) -> sdt 

#write.csv(sdt, file = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/sdt.csv", row.names = FALSE)
rm(PS.eimeria, sdtEim)

### Extract information by primer pair
##Functions
sumSeqByTax <- function(PS.l, rank) {
  lapply(PS.l, function(x){ 
    counts <- data.frame(cbind(asvCount=colSums(otu_table(x)), tax_table(x)))
    counts$asvCount <- as.numeric(as.character(counts$asvCount))
    tapply(counts$asvCount, counts[, rank], sum)})
}

##Total reads by primer pair 
counts_multi_primers <- as.data.frame(unlist(lapply(PS.l, function(x){
  data.frame(cbind(asvCount=colSums(data.frame(rowSums(otu_table(x))))))
})))

counts_multi_primers[,2] <- rownames(counts_multi_primers)
rownames(counts_multi_primers) <- c(1:nrow(counts_multi_primers))
colnames(counts_multi_primers) <- c("TotalReads_primer","Primer_comb_ID")
counts_multi_primers <- data.frame(Primer_comb_ID = counts_multi_primers$Primer_comb_ID, TotalReads_primer = counts_multi_primers$TotalReads_primer) 
counts_multi_primers$Primer_comb_ID <- gsub(".asvCount", "\\1", counts_multi_primers$Primer_comb_ID)

##Total ASVs by primer pair
asvcounts <- as.data.frame(unlist(lapply(PS.l, function(x){
  data.frame(cbind(asvnumber=ncol((otu_table(x)))))
})))

asvcounts[,2] <- rownames(asvcounts)
rownames(asvcounts) <- c(1:nrow(asvcounts))
colnames(asvcounts) <- c("Number_ASVs","Primer_comb_ID")
asvcounts <- data.frame(Primer_comb_ID = asvcounts$Primer_comb_ID, Number_ASVs = asvcounts$Number_ASVs) 
asvcounts$Primer_comb_ID <- gsub(".asvnumber", "\\1", asvcounts$Primer_comb_ID)

## Merge both 
counts_multi_primers <- join(counts_multi_primers, asvcounts, by= "Primer_comb_ID")
rm(asvcounts)

##Create a table with total number of Reads and ASVs per primer pair
#write.csv(counts_multi_primers, file = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/Reads_per_Primer_Pair.csv")

##Extract sample counts information by primer combination 
samplecounts<- data.frame()

for (i in 1: length(PS.l)) ### Start a loop: for every element in the list ...
{ 
  tmp <- data.frame() #### make an individual data frame ...
  
  {
    tmp <- as.data.frame(sapply(PS.l[i], function(x) sample_sums(x)))  ###Make a data frame with the sample sums for each primer combination 
    tmp[,2]<-rownames(tmp) ### And use the rownames as information of the second column
    tmp[,3] <- names(PS.l[i])
    colnames(tmp)<- c("TotalReads_primer","labels", "Primer_comb_ID")
    rownames(tmp) <- c(1:nrow(tmp))
  }
  samplecounts <- rbind(samplecounts, tmp) ### Join all the "individual" data frames into the list 
}

rm(tmp, i)

##Check which primers actually amplified something assigned as Eimeria 

readNumByGenus <- sumSeqByTax(PS.l = PS.l, rank = "genus")

AbGen <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByGenus)) ### Start a loop: for every element in the list ...
{ 
  genus <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByGenus[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    genus[1,1] <- 0    ### Add a zero in the first column
    genus[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    genus <- as.data.frame((readNumByGenus[[i]]))  ###Make a data frame with the data included in each element of the list 
    genus[,2] <- rownames(genus) ### And use the rownames as information of the second column 
  }
  
  genus[,3] <- names(readNumByGenus)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(genus) <- c("Read_count", "Genus", "Primer_comb_ID") ### change the names for the columns 
  AbGen <- rbind(AbGen, genus) ### Join all the "individual" data frames into the final data frame 
  
}   ### close loop

rownames(AbGen) <- c(1:nrow(AbGen)) ### change the rownames to consecutive numbers 
AbGen <- data.frame(Primer_comb_ID = AbGen$Primer_comb_ID, Genus = AbGen$Genus, Read_count = AbGen$Read_count) ###change the order of the columns
AbGen$Read_count <- as.numeric(AbGen$Read_count)

readNumByGenus <- sumSeqByTax(PS.l, rank = "genus")

AbGen <- data.frame() ###Create the data frame 
for (i in 1: length(readNumByGenus)) ### Start a loop: fro every element in the list ...
{ 
  genus <- data.frame() #### make an individual data frame ...
  if(nrow(as.data.frame(readNumByGenus[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    genus[1,1] <- 0    ### Add a zero in the first column
    genus[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    genus <- as.data.frame((readNumByGenus[[i]]))  ###Make a data frame with the data included in each element of the list 
    genus[,2] <- rownames(genus) ### And use the rownames as information of the second column 
  }
  genus[,3] <- names(readNumByGenus)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(genus) <- c("Read_count", "Genus", "Primer_comb_ID") ### change the names for the columns 
  AbGen <- rbind(AbGen, genus) ### Join all the "individual" data frames into the final data frame 
}   ### close loop
rownames(AbGen) <- c(1:nrow(AbGen)) ### change the rownames to consecutive numbers 
AbGen <- data.frame(Primer_comb_ID = AbGen$Primer_comb_ID, Genus = AbGen$Genus, Read_count = AbGen$Read_count) ###change the order of the columns
AbGen$Read_count <- as.numeric(AbGen$Read_count)

AbGen %>%
  group_by(Genus)%>%
  filter(Genus=="Eimeria")-> AbGen_Eimeria 

rm(genus, i, readNumByGenus)
##At species level which ones are amplifing for E. ferrisi

readNumBySpecies <- sumSeqByTax(PS.l, rank = "species")

AbSp <- data.frame() ###Create the data frame 
for (i in 1: length(readNumBySpecies)) ### Start a loop: fro every element in the list ...
{ 
  species <- data.frame() #### make an individual data frame ...
  if(nrow(as.data.frame(readNumBySpecies[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    species[1,1] <- 0    ### Add a zero in the first column
    species[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    species <- as.data.frame((readNumBySpecies[[i]]))  ###Make a data frame with the data included in each element of the list 
    species[,2] <- rownames(species) ### And use the rownames as information of the second column 
  }
  
  species[,3] <- names(readNumBySpecies)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(species) <- c("Read_count", "Species", "Primer_comb_ID") ### change the names for the columns 
  AbSp <- rbind(AbSp, species) ### Join all the "individual" data frames into the final data frame 
  
}   ### close loop

rownames(AbSp) <- c(1:nrow(AbSp)) ### change the rownames to consecutive numbers 
AbSp <- data.frame(Primer_comb_ID = AbSp$Primer_comb_ID, Species = AbSp$Species, Read_count = AbSp$Read_count) ###change the order of the columns
AbSp$Read_count <- as.numeric(AbSp$Read_count)

AbSp %>%
  group_by(Species)%>%
  filter(Species=="Eimeria ferrisi")-> AbSp_Eimfer 

rm(species, i, readNumBySpecies)

### Subset information from Wang1141_13.Nem_0425_6 combination to compair with single amplicon amplification of the same

PS.18S.Multi <- phyloseq(otu_table(PS.l$wang1141_13_F.Nem_0425_6_3_R), 
                   sample_data(PS.l$wang1141_13_F.Nem_0425_6_3_R), 
                   tax_table(PS.l$wang1141_13_F.Nem_0425_6_3_R))

sum(otu_table(PS.18S.Multi)) ##Total denoised reads = 450,144 reads

counts_multi_18S <- data.frame(rowSums(otu_table(PS.18S.Multi)))
counts_multi_18S[,2] <- rownames(counts_multi_18S)
###Eukaryotic counts
counts_multi_18S[,3] <- as.data.frame(rowSums(otu_table(subset_taxa(PS.18S.Multi, superkingdom%in%"Eukaryota"))))
counts_multi_18S[,4] <- as.data.frame(rowSums(otu_table(subset_taxa(PS.18S.Multi, genus%in%"Eimeria"))))
colnames(counts_multi_18S) <- c("Read_counts_18S", "labels", "Eukaryota_18S_reads", "ReadsEim18SMulti")
rownames(counts_multi_18S) <- c(1:nrow(counts_multi_18S))
counts_multi_18S <- data.frame(labels = counts_multi_18S$labels, 
                         Read_counts_18S = counts_multi_18S$Read_counts_18S,
                         Eukaryota_18S_reads= counts_multi_18S$Eukaryota_18S_reads,
                         ReadsEim18SMulti= counts_multi_18S$ReadsEim18SMulti)

counts_multi_18S%>%
  dplyr::select(labels, Read_counts_18S,ReadsEim18SMulti)%>%
  dplyr::mutate(Eimeria_abundance_18S= ReadsEim18SMulti/Read_counts_18S)-> sdt18SEimMulti


