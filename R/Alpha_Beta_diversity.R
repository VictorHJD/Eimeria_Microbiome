library("ggplot2")
library("microbiome")
library("phyloseq")
require(tidyr)
data(dietswap)
pseq <- dietswap ### Use here Julia's dataset 


###Alpha diversity 
richness <- as.data.frame(estimate_richness(pseq, measures=c("Observed", "InvSimpson", "Shannon", "Chao1")))
plot_richness(pseq)
plot_richness(pseq, x="timepoint", measures=c("Chao1", "Shannon")) ###plot alpha diversity by DPI (Julia's dataset)

####Add richness alpha diversity data to phyloseq 
sample_data(pseq)$Shannon <- richness$Shannon
sample_data(pseq)$InvSimpson <- richness$InvSimpson

pseqdata <- meta(pseq)


###Compair differences in Shannon idex between groups (e.g. time points)
# create a list of pairwise comaprisons
#bmi <- levels(pseqdata$bmi_group) # get the variables

# make a pairwise list that we want to compare.
#bmi.pairs <- combn(seq_along(bmi), 2, simplify = FALSE, FUN = function(i)bmi[i])
#stat_compare_means(comparisons = bmi.pairs) 

###with ggplot
ggplot(pseqdata, aes(x= as.factor(timepoint), y= Shannon, fill= as.factor(timepoint)))+
  geom_boxplot()+
  geom_jitter()+
  facet_grid(bmi_group~.)+ 
  theme_classic()



###Visualize the data
plot_regression(Shannon ~ timepoint, meta(pseq)) + xlab("Time point") + ylab("Diversity (Shannon)")


###
plot_landscape(pseq, "MDS", "bray", col = "timepoint")

###
dominance(pseq)
rarity(pseq)
evenness(pseq)

###Transform to compositional data 
pseqcomp <- microbiome::transform(pseq, "compositional")


# Extract abundance matrix from the phyloseq object, transpose taxa to colum and sample to rows, and Coerce to data.frame
OTUdf <- as.data.frame(t(as(otu_table(pseq), "matrix")))
rownames(OTUdf) <- gsub("-", ".", rownames(OTUdf))

# Extract sample data
sample <- as(sample_data(pseq), "data.frame")
rownames(sample) <- gsub("-", ".", rownames(sample))


# Merge richeness, OTU and sample 
richness <- merge(richness, sample, by= "row.names")
row.names(richness)<- richness$Row.names
richness[,1]<- NULL

richness <- richness[order(richness$timepoint),]

ggplot(richness, aes(x=timepoint, y=Shannon)) +
  #geom_line(aes(color=sex))+
  geom_point(aes(color=bmi_group))+
  geom_smooth(aes(color=bmi_group))
#facet_grid(sex~.)+ theme(legend.position="none")


####Mean, 95%CI
df_tidy_mean <- df_tidy %>%
  filter(!is.na(Ratio)) %>%
  group_by(Time, id) %>%
  summarise(n = n(),
            mean = mean(Ratio),
            median = median(Ratio),
            sd = sd(Ratio)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)


###Beta diversity
ordi <- ordinate(pseq, method="MDS", distance="bray")

plot_ordination(pseq, ordi, color= "timepoint", shape = "bmi_group")+
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=bmi_group))


#sample_data(PSFox)$readsum <-  sample_sums(PSFox)
#plot_ordination(PSFox, ordi, color="readsum", label="readsum")
