library(ggplot2)
library(MultiAmplicon)
library(reshape)
library(phyloseq)
library(data.table)
library(taxonomizr)
library(taxize)
library(parallel)

if(!exists("MA")){
  source("MA_General.R")
}

#######Load data
MA <- readRDS("/SAN/Victors_playground/Eimeria_microbiome/MASample.Rds")

##Extract sequences from the 3 primers 
SeqV1V2 <- getSequences(MA@sequenceTableNoChime$`27M_F_98_F.Klin0341_CR_18_R`)
names(SeqV1V2) <- SeqV1V2

SeqV3V4 <- getSequences(MA@sequenceTableNoChime$`515F_Y_118_F.806R_118_R`)
names(SeqV3V4) <- SeqV3V4

Seq18S <- getSequences(MA@sequenceTableNoChime$wang1141_13_F.Nem_0425_6_3_R)
names(Seq18S) <- Seq18S

##Alignments
require(DECIPHER)

AlignV1V2<- AlignSeqs(DNAStringSet(SeqV1V2), anchor= NA, verbose= FALSE)
AlignV3V4<- AlignSeqs(DNAStringSet(SeqV3V4), anchor=NA, verbose = FALSE)
Align18S<- AlignSeqs(DNAStringSet(Seq18S), anchor=NA, verbose = FALSE)

##Trees
require(phangorn)

phangAlignV1V2 <- phyDat(as(AlignV1V2, "matrix"), type="DNA")
dmV1V2 <- dist.ml(phangAlignV1V2) ## Distance matrix
treeNJV1V2 <- NJ(dmV1V2) # Note, tip order != sequence order
plot(treeNJV1V2, type= "unrooted", use.edge.length= FALSE, no.margin= TRUE, show.tip.label= FALSE)
fit1 <- pml(treeNJV1V2, data=phangAlignV1V2)
fitGTRV1V2 <- update(fit1, k=4, inv=0.2)
fitGTRV1V2 <- optim.pml(fitGTRV1V2, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
plot(fitGTRV1V2, type= "unrooted", use.edge.length= FALSE, no.margin= TRUE, show.tip.label= FALSE)


phangAlignV3V4 <- phyDat(as(AlignV3V4, "matrix"), type="DNA")
dmV3V4 <- dist.ml(phangAlignV3V4)
treeNJV3V4 <- NJ(dmV3V4) # Note, tip order != sequence order
plot(treeNJV3V4, type= "unrooted", use.edge.length= FALSE, no.margin= TRUE, show.tip.label= FALSE)
fit2 <- pml(treeNJV3V4, data=phangAlignV3V4)
fitGTRV3V4<- update(fit2, k=4, inv=0.2)
fitGTRV3V4 <- optim.pml(fitGTRV3V4, model="GTR", optInv=TRUE, optGamma=TRUE,
                        rearrangement = "stochastic", control = pml.control(trace = 0))

plot(fitGTRV3V4, type= "unrooted", use.edge.length= FALSE, no.margin= TRUE, show.tip.label= FALSE)

phangAlign18S <- phyDat(as(Align18S, "matrix"), type="DNA")
dm18S <- dist.ml(phangAlign18S)
treeNJ18S <- NJ(dm18S) # Note, tip order != sequence order
plot(treeNJ18S, type= "unrooted", use.edge.length= FALSE, no.margin= TRUE, show.tip.label= FALSE)
fit3 <- pml(treeNJ18S, data=phangAlign18S)
fitGTR18S<- update(fit3, k=4, inv=0.2)
fitGTR18S <- optim.pml(fitGTR18S, model="GTR", optInv=TRUE, optGamma=TRUE,
                        rearrangement = "stochastic", control = pml.control(trace = 0))

plot(fitGTR18S, type= "unrooted", use.edge.length= FALSE, no.margin= TRUE, show.tip.label= FALSE)
#detach("package:phangorn", unload=TRUE)

###Create individual phyloseq objects by primer pair 

PS.16SA <- phyloseq(otu_table(PS.l$`27M_F_98_F.Klin0341_CR_18_R`),
                    sample_data(PS.l$`27M_F_98_F.Klin0341_CR_18_R`), 
                    tax_table(PS.l$`27M_F_98_F.Klin0341_CR_18_R`), 
                    phy_tree(fitGTRV1V2$tree))

#saveRDS(PS.16SA, file="/SAN/Victors_playground/Eimeria_microbiome/PS_16S_A.Rds") ###Information from 16S V1-V2


PS.16SB <- phyloseq(otu_table(PS.l$`515F_Y_118_F.806R_118_R`), 
                    sample_data(PS.l$`515F_Y_118_F.806R_118_R`), 
                    tax_table(PS.l$`515F_Y_118_F.806R_118_R`), 
                    phy_tree(fitGTRV3V4$tree))

#saveRDS(PS.16SB, file="/SAN/Victors_playground/Eimeria_microbiome/PS_16S_B.Rds") ###Information from 16S V3-V4

PS.18S <- phyloseq(otu_table(PS.l$wang1141_13_F.Nem_0425_6_3_R), 
                   sample_data(PS.l$wang1141_13_F.Nem_0425_6_3_R), 
                   tax_table(PS.l$wang1141_13_F.Nem_0425_6_3_R), 
                   phy_tree(fitGTR18S$tree))

#saveRDS(PS.18S, file="/SAN/Victors_playground/Eimeria_microbiome/PS_18S.Rds") ###Information from 18S
