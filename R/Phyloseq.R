###Alpha diversity graphics

library(ggplot2)
library(phyloseq)
library(vegan)
library("microbiome")
library(data.table)
library(gridExtra)
library(grid)

if(!exists("PS")){
  source("MA_General.R")
}

##############Load data######################
#ps.l <- ps.l[!grepl("Klin|ADM", names(ps.l))] just 16S data 

PS <- readRDS(file="/SAN/Victors_playground/Eimeria_microbiome/PhyloSeqData.Rds")
PS.16SA <- readRDS(file = "/SAN/Victors_playground/Eimeria_microbiome/PS_16S_A.Rds")
PS.16SB <- readRDS(file = "/SAN/Victors_playground/Eimeria_microbiome/PS_16S_B.Rds")
PS.18S<- readRDS(file = "/SAN/Victors_playground/Eimeria_microbiome/PS_18S.Rds")

###General check of data
summarize_phyloseq(PS.18S) ##Change dataset 
rank_names(PS)
tax_table(PS)

##Taxonomic filtering
# Create table, number of features for each phyla
table(tax_table(PS.18S)[, "phylum"], exclude = NULL)

###Remove features with ambiguous phylum annotation
PS <- subset_taxa(PS, !is.na(phylum) & !phylum %in% c("", "uncharacterized"))
PS.16SA <- subset_taxa(PS.16SA, !is.na(phylum) & !phylum %in% c("", "uncharacterized"))
PS.16SB <- subset_taxa(PS.16SB, !is.na(phylum) & !phylum %in% c("", "uncharacterized"))
PS.18S <- subset_taxa(PS.18S, !is.na(phylum) & !phylum %in% c("", "uncharacterized"))

# Compute prevalence of each feature, store as data.frame
PrevAll <- apply(X = otu_table(PS),
               MARGIN = ifelse(taxa_are_rows(PS), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

Prev16SA <- apply(X = otu_table(PS.16SA),
                  MARGIN = ifelse(taxa_are_rows(PS.16SA), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})

Prev16SB <- apply(X = otu_table(PS.16SB),
                  MARGIN = ifelse(taxa_are_rows(PS.16SB), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})

Prev18S <- apply(X = otu_table(PS.18S),
                  MARGIN = ifelse(taxa_are_rows(PS.18S), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
PrevAll <- data.frame(Prevalence = PrevAll,
                    TotalAbundance = taxa_sums(PS),
                    tax_table(PS))

Prev16SA <- data.frame(Prevalence = Prev16SA,
                      TotalAbundance = taxa_sums(PS.16SA),
                      tax_table(PS.16SA))

Prev16SB <- data.frame(Prevalence = Prev16SB,
                      TotalAbundance = taxa_sums(PS.16SB),
                      tax_table(PS.16SB))

Prev18S <- data.frame(Prevalence = Prev18S,
                      TotalAbundance = taxa_sums(PS.18S),
                      tax_table(PS.18S))


plyr::ddply(PrevAll, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
plyr::ddply(Prev16SA, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
plyr::ddply(Prev16SB, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
plyr::ddply(Prev18S, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Define phyla to filter (present in less than 1% of samples)
filterPhyla <- c("Basidiomycota", "Mucoromycota") ###For the full dataset
filterPhyla <- c("Streptophyta", "Verrucomicrobia") ###For 16SA
###For 16SB nothing should be eliminated
filterPhyla <- c("Basidiomycota", "Mucoromycota") ###For 18S


# Filter entries with unidentified Phylum.
PS1 <- subset_taxa(PS, !phylum %in% filterPhyla)
PS1.16SA <- subset_taxa(PS.16SA, !phylum %in% filterPhyla)
PS1.18S <- subset_taxa(PS.18S, !phylum %in% filterPhyla)

####Prevalence filtering 
# Subset to the remaining phyla
PrevAll2 <- subset(PrevAll, phylum %in% get_taxa_unique(PS1, "phylum"))
ggplot(PrevAll2, aes(TotalAbundance, Prevalence / nsamples(PS),color=phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance (read count)") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none")

Prev16SA2 <- subset(Prev16SA, phylum %in% get_taxa_unique(PS1.16SA, "phylum"))
ggplot(Prev16SA2, aes(TotalAbundance, Prevalence / nsamples(PS.16SA),color=phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance (read count)") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none")

Prev16SB2 <- subset(Prev16SB, phylum %in% get_taxa_unique(PS.16SB, "phylum"))
ggplot(Prev16SB2, aes(TotalAbundance, Prevalence / nsamples(PS.16SA),color=phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance (read count)") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none")

Prev18S2 <- subset(Prev18S, phylum %in% get_taxa_unique(PS1.18S, "phylum"))
ggplot(Prev18S2, aes(TotalAbundance, Prevalence / nsamples(PS.18S),color=phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance (read count)") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none")


# Define prevalence threshold as 2% of total samples
prevalenceThreshold <- 0.02 * nsamples(PS)
prevalenceThreshold


# How many genera would be present after filtering?
length(get_taxa_unique(PS, taxonomic.rank = "genus"))


###Tree plots
plot_tree(physeq = PS1.16SA, method = "treeonly",
          ladderize = "left",
          title = "Before Agglomeration") 

plot_tree(PS1.16SA, color = "phylum", size= "abundance")

###Seq depths
qplot(log10(rowSums(otu_table(PS1))),binwidth=0.2) +
  xlab("Log10 counts-per-sample")+ ylab("Count")+
  theme_bw()+
  labs(tag = "A)")


qplot(log10(rowSums(otu_table(PS1.16SA))),binwidth=0.2) +
  xlab("Log10 counts-per-sample")+ ylab("Count")+
  theme_bw()+
  labs(tag = "B)")


qplot(log10(rowSums(otu_table(PS.16SB))),binwidth=0.2) +
  xlab("Log10 counts-per-sample")+ ylab("Count")+
  theme_bw()+
  labs(tag = "C)")

qplot(log10(rowSums(otu_table(PS1.18S))),binwidth=0.2) +
  xlab("Log10 counts-per-sample")+ ylab("Count")+
  theme_bw()+
  labs(tag = "D)")

####Normalizing data :S 
PSlog <- transform_sample_counts(PS, function(x) log(1 + x))
PS16SAlog <- transform_sample_counts(PS1.16SA, function(x) log(1 + x))
PS16SBlog <- transform_sample_counts(PS.16SB, function(x) log(1 + x))
PS18Slog <- transform_sample_counts(PS1.18S, function(x) log(1 + x))

##Plot abundances (alpha diversity)
sample_data(PS1)$dpi <- as.factor(sample_data(PS1)$dpi)

plot_bar(PS1, fill="phylum") + facet_wrap(~dpi, scales= "free_x", nrow=1)
plot_richness(PS1, x="dpi", color="dpi", measures=c("Chao1", "Shannon", "Simpson"))+ 
  geom_boxplot()+
  theme_bw()

plot_richness(PS1.16SB, x="dpi", color="dpi", measures=c("Chao1", "Shannon", "Simpson"))+ 
  geom_boxplot()+
  theme_bw()


alphaDiv <- estimate_richness(PS1)
alphaDiv16SA <- estimate_richness(PS1.16SA)
alphaDiv16SB <- estimate_richness(PS1.16SB)
alphaDiv18S <- estimate_richness(PS.18S)

eveAll <- evenness(PS1, c("pielou", "simpson", "evar", "bulla"))

###Difference accross DPI?
pairwise.wilcox.test(alphaDiv16SB$Shannon, sample_data(PS1.16SB)$dpi)

#Beta diversity 

##PCoA plot using the unweighted UniFrac as distance (just by marker because it requires phylogenetic tree data)
sample_data(PS1.16SA)$dpi <- as.factor(sample_data(PS1.16SA)$dpi)
PS1.16SA <- prune_samples(sample_sums(PS.16SA)>1, PS.16SA)


wunifrac_dist <- phyloseq::distance(PS1.16SA, method="unifrac", weighted=F)
ordination <- ordinate(PS1.16SA, method="PCoA", distance=wunifrac_dist)
plot_ordination(PS1.16SA, ordination, color="dpi") + 
  theme(aspect.ratio=1)+
  theme_bw()#+
  #stat_ellipse()

sample_data(PS.16SB)$dpi <- as.factor(sample_data(PS.16SB)$dpi)
PS1.16SB <- prune_samples(sample_sums(PS.16SB)>1, PS.16SB)


wunifrac_dist <- phyloseq::distance(PS1.16SB, method="unifrac", weighted=F)
ordination <- ordinate(PS1.16SB, method="PCoA", distance=wunifrac_dist)
plot_ordination(PS1.16SB, ordination, color="dpi") + theme(aspect.ratio=1)+
  theme_bw()#+
  #stat_ellipse()

sample_data(PS1.18S)$dpi <- as.factor(sample_data(PS1.18S)$dpi)
PS1.18S <- prune_samples(sample_sums(PS1.18S)>1, PS1.18S)

wunifrac_dist <- phyloseq::distance(PS1.18S, method="unifrac", weighted=F)
ordination <- ordinate(PS1.18S, method="PCoA", distance=wunifrac_dist)
plot_ordination(PS1.18S, ordination, color="dpi") + theme(aspect.ratio=1)+
  theme_bw()#+
  #stat_ellipse()


##Permanova
require("vegan")
adonis(wunifrac_dist ~ sample_data(PS1.18S)$dpi)

###Using Negative Bionomial in microbiome differential abundance testing (not working check it again :( )


require("DESeq2")
dipdds <- phyloseq_to_deseq2(PS1.16SA, ~ dpi)
DESeq(dipdds)

# calculate geometric means prior to estimate size factors
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(dipdds), 1, gm_mean)
dipdds <- estimateSizeFactors(dipdds, geoMeans)
dipdds <- DESeq(dipdds, fitType="local")

###Ordination
OrdPS16SA <- ordinate(PS16SAlog, method = "MDS", distance = "wunifrac")
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2] / evals[1]))

###Composition
require("dplyr")
PS1c <- PS1 %>%
  aggregate_taxa(level = "phylum") %>%  
  microbiome::transform(transform = "compositional")

PS1c %>%
plot_composition(sample.sort = "Eimeria", otu.sort = "abundance", group_by= "dpi") + theme(legend.position="none")

plot_composition(PS1c,
                 average_by = "dpi")+ theme(legend.position="bottom")


PS1.16SAc <- PS1.16SA %>%
  aggregate_taxa(level = "family") %>%  
  microbiome::transform(transform = "compositional")

g<- plot_composition(PS1.16SAc,average_by = "dpi")+ 
    theme(legend.position="bottom")+
    xlab("Day post infection")+
    theme_bw()+
    labs(tag = "A)")


PS1.16SBc <- PS1.16SB %>%
  aggregate_taxa(level = "family") %>%  
  microbiome::transform(transform = "compositional")

h<- plot_composition(PS1.16SBc,average_by = "dpi")+ 
    theme(legend.position="bottom")+
    xlab("Day post infection")+
    theme_bw()+
    labs(tag = "B)")

PS1.18Sc <- PS1.18S %>%
  aggregate_taxa(level = "family") %>%  
  microbiome::transform(transform = "compositional")

i<- plot_composition(PS1.18Sc, average_by = "dpi")+ 
    theme(legend.position="bottom")+
    xlab("Day post infection")+
    theme_bw()+
    labs(tag = "C)")

pdf(file = "~/AA_Microbiome/Figures/Composition_Infection.pdf", width = 15, height = 15)
grid.arrange(g,h,i, ncol= 1, nrow= 3)
dev.off()

###Summarize seq depth 
sdt <- data.table(as(sample_data(PS1), "data.frame"),
                    TotalReads= sample_sums(PS1), keep.rownames = T)

setnames(sdt, "rn", "Sample_ID")

###Summarize seq depth 18S
sdt.18S <- data.table(as(sample_data(PS1.18S), "data.frame"),
                  TotalReads= sample_sums(PS1.18S), keep.rownames = T)

setnames(sdt.18S, "rn", "Sample_ID")

saveRDS(PS1.18S, file = "/SAN/Victors_playground/Eimeria_microbiome/PS1_18S.Rds") ##Filtered 
##Check for Eimeria
PS.eimeria <- subset_taxa(PS1, genus%in%"Eimeria")

###Eimeria 
sdtEim <- data.table(as(sample_data(PS.eimeria), "data.frame"),
                     ReadsEim= sample_sums(PS.eimeria), keep.rownames = T)


sdtEim <- dplyr::select(sdtEim, 5,57)

require(dplyr)
sdt %>%
  mutate(Eimeria_abundance = sdt$ReadsEim/sdt$TotalReads) -> sdt ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

####OPG vs reads Eimeria 
require(ggpubr)
opgre <- ggplot(sdt, aes(OPG, ReadsEim))+
  geom_point(size= 3)+
  geom_smooth(method = lm)+
  xlab("Oocyst of Eimeria per gram of feces")+
  ylab("Sequence reads count (Eimeria)")+
  labs(tag= "A)")+
  theme_bw()+
  stat_cor(label.x = 3000000, label.y = 2000, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 3000000, label.y = 2300)

opgea <- ggplot(sdt, aes(OPG, Eimeria_abundance))+
  geom_point(size= 3)+
  geom_smooth(method = lm)+
  xlab("Oocyst of Eimeria per gram of feces")+
  ylab("Relative abundance of Eimeria (Sequencing reads)")+
  labs(tag= "B)")+
  theme_bw()+
  stat_cor(label.x = 3000000, label.y = .2, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 3000000, label.y = .25)

pdf(file = "~/AA_Microbiome/Figures/OPG_Eimeria_Reads.pdf", width = 15, height = 15)
grid.arrange(opgre,opgea, ncol= 1, nrow= 2)
dev.off()

###Course of infection 
compare_means(OPG ~ dpi,  data = sdt)

a<- ggplot(sdt, aes(dpi, OPG, fill= dpi))+
    geom_boxplot()+
    xlab("Day post infection")+
    ylab("Oocyst of Eimeria per gram of feces")+
    labs(tag= "A)")+
    theme_bw()

b<- ggplot(sdt, aes(dpi, ReadsEim, fill= dpi))+
    geom_boxplot()+
    xlab("Day post infection")+
    ylab("Count of Eimeria assigned reads")+
    labs(tag= "B)")+
    theme_bw()

c<- ggplot(sdt, aes(dpi, Eimeria_abundance, fill= dpi))+
    geom_boxplot()+
    xlab("Day post infection")+
    ylab("Relative proportion of Eimeria assigned reads")+
    labs(tag= "C)")+
    theme_bw()

require("grid")
require("gridExtra")
pdf(file = "~/AA_Microbiome/Figures/Course_of_Eimeria_Infection.pdf", width = 15, height = 15)
grid.arrange(a,b,c, ncol= 1, nrow= 3)
dev.off()

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
