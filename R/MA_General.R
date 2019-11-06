## Please uncomment the first time you run this and re-install packages

## require(devtools)
## devtools::install_github("derele/MultiAmplicon", force= T)
## devtools::install_github("derele/dada2", force= T)

library(ggplot2)
library(MultiAmplicon)
library(reshape)
library(phyloseq)
library(data.table)
library(taxonomizr)
library(taxize)
library(parallel)

## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doFilter <- TRUE

doMultiAmp <- TRUE

doTax <- TRUE
## But remember: if you change the MultiAmplicon Analysis, the
## taxonomic annotation might be out of sync...

###################Full run Microbiome#######################
#Preparation of files

##These are the same steps that are followed by the DADA2 pipeline

path <- "/SAN/Victors_playground/Eimeria_microbiome/2018_22_Mmb_1/" ## change according to where you downloaded
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 


samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
samples<- gsub("S\\d+-", "\\1", basename(samples))

#Extra step in the pipeline: quality plots of the reads 
## plotQualityProfile(fastqF[[1]])
## plotQualityProfile(fastqF[[2]])
## plotQualityProfile(fastqR[[1]])
## plotQualityProfile(fastqR[[2]])

#Creation of a folder for filtrated reads 

filt_path <- "/SAN/Victors_playground/Eimeria_microbiome/filtered/"

#Pipeline filtration 
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

## some files will be filtered out completely, therefore allowing 50
## files less present and still don't redo filtering
if(doFilter){
  lapply(seq_along(fastqF),  function (i) {
    filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                  truncLen=c(250,250), minLen=c(250,250), 
                  maxN=0, maxEE=2, truncQ=2, rm.phix = TRUE,
                  compress=TRUE, verbose=TRUE)
  })
}


names(filtFs) <- names(filtRs) <- samples
files <- PairedReadFileSet(filtFs, filtRs)

#Preparation of primer file ### Here stats the Multiamplicon pipeline from Emanuel

#Primers used in the arrays 
ptable <- read.csv(file = "/SAN/Victors_playground/Eimeria_microbiome/primer.file.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "TS.SequenceF"]
primerR <- ptable[, "TS.SequenceR"]
names(primerF) <- as.character(ptable[, "corrected.NameF"])
names(primerR) <- as.character(ptable[, "corrected.NameR"])

primer <- PrimerPairsSet(primerF, primerR)

##Multi amplicon pipeline
if(doMultiAmp){
  MA <- MultiAmplicon(primer, files)
  filedir <- "/SAN/Victors_playground/Eimeria_microbiome/Stratified_files"
  if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
  MA <- sortAmplicons(MA, n=1e+05, filedir=filedir) ## This step sort the reads into amplicons based on the number of primer pairs
  
  errF <-  learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                       verbose=0, multithread = 12)
  errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                      verbose=0, multithread = 12)
  
  MA <- derepMulti(MA, mc.cores=12) 
  MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                   verbose=0, mc.cores=12)
  MA <- mergeMulti(MA, mc.cores=12) 
  
  propMerged <- MultiAmplicon::calcPropMerged(MA)
  
  MA <- mergeMulti(MA, justConcatenate=propMerged<0.8, mc.cores=12) 
  
  MA <- makeSequenceTableMulti(MA, mc.cores=12) 
  
  MA <- removeChimeraMulti(MA, mc.cores=12)
  
  saveRDS(MA, "/SAN/Victors_playground/Eimeria_microbiome/MA.RDS")
} else{
  MA <- readRDS("/SAN/Victors_playground/Eimeria_microbiome/MA.RDS")
}

trackingF <- getPipelineSummary(MA) 
## plotPipelineSummary(trackingF) 
## plotPipelineSummary(trackingF) + scale_y_log10()


plotAmpliconNumbers(MA)

###New taxonomic assignment 

MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    infasta = "/SAN/Victors_playground/Eimeria_microbiome/in.fasta",
                    outblast = "/SAN/Victors_playground/Eimeria_microbiome/out.fasta",
                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                    num_threads = 20)

saveRDS(MA, file="/SAN/Victors_playground/Eimeria_microbiome/MATax.Rds")


### Add sample information
require(phyloseq)
sample.data <- read.csv("/SAN/Victors_playground/Eimeria_microbiome/sample_data_infb_Exp005.csv")
exp.des <- read.csv("/SAN/Victors_playground/Eimeria_microbiome/Inf1b_Exp005.DESIGN.csv")
exp.des$InfectionStrain <- NULL
exp.des$Departure <- NULL
exp.des$NOTE <- NULL

###Use function from Alice to calculate OPG
calculateOPG <- function(sample.data){
  sample.data$mean_Neubauer <- 
    (sample.data$oocyst_sq1 + sample.data$oocyst_sq2 + sample.data$oocyst_sq3 + sample.data$oocyst_sq4) / 4
  # NB! Limit of detection = 1 oocysts
  sample.data$mean_Neubauer[sample.data$oocyst_sq1 + sample.data$oocyst_sq2 + sample.data$oocyst_sq3 + sample.data$oocyst_sq4 == 1] <- 0
  sample.data$oocysts.per.tube <- sample.data$mean_Neubauer * 10000 * sample.data$dilution
  sample.data$OPG <- sample.data$oocysts.per.tube / sample.data$fecweight_flot
  ## If we don't have the fecal weight BUT we counted in Neubauer chamber 0, then OPG = 0
  sample.data$oocysts.per.tube[sample.data$fecweight_flot == 0 & sample.data$mean_Neubauer == 0] <- 0
  sample.data$OPG[sample.data$fecweight_flot == 0 & sample.data$mean_Neubauer == 0] <- 0
  return(sample.data)
}

sample.data<- calculateOPG(sample.data = sample.data)

sample.data <- merge(sample.data, exp.des, by= "EH_ID", all= TRUE) ##merge sample data with genotype of mice 
sample.data <- sample.data[-243,] ##eliminate empty row
rownames(sample.data) <- sample.data$labels

###Estimate microbiota density (MD) from Contijoch et al. 2019
### MD = total DNA per sample (µg)/ mg of fresh feces

require(dplyr)

##considering 30µL of elution volume 
sample.data %>%
  mutate(Total_DNA = (sample.data$Conc_DNA*30)*0.001) -> sample.data ### Add a new variable that will contain total DNA extracted per sample in µg

sample.data %>%
  mutate(Microbial_density = sample.data$Total_DNA/(sample.data$fecweight_DNA*1000)) -> sample.data ### Total DNA extracted per sample in µg by feces weight in mg


sample.data$labels<- as.vector(sample.data$labels)
rownames(sample.data) <- make.unique(sample.data$labels)
MA <- addSampleData(MA, sample.data)

saveRDS(MA, file="/SAN/Victors_playground/Eimeria_microbiome/MASample.Rds")  ###START from here now! 

####Raw counts
rawcounts <- data.frame(colSums(getRawCounts(MA)))
rawcounts[,2] <- rownames(rawcounts)
colnames(rawcounts) <- c("Raw_counts", "labels")
rownames(rawcounts) <- c(1:nrow(rawcounts))
rawcounts <- data.frame(labels = rawcounts$labels, Raw_counts = rawcounts$Raw_counts) 

hist(rawcounts$Raw_counts)
summary(rawcounts$Raw_counts)
boxplot(rawcounts$Raw_counts)

sum(rawcounts$Raw_counts) ###Total raw counts reads = 2789639

##To phyloseq
##Sample data
PS <- toPhyloseq(MA, colnames(MA))
sum(otu_table(PS)) ##Total denoised reads = 2015804
##Primer data
PS.l <- toPhyloseq(MA, colnames(MA),  multi2Single=FALSE)

saveRDS(PS.l, file="/SAN/Victors_playground/Eimeria_microbiome/PhyloSeqList.Rds") ###For primer analysis (Victor)
saveRDS(PS, file="/SAN/Victors_playground/Eimeria_microbiome/PhyloSeqData.Rds") ###For Sample analysis (Julia and Victor)

####Old taxonomic assignment
###Extract sequences to do taxonomic assignment 

#STNC <- getSequenceTableNoChime(MA)

#sequences <- unlist(lapply(STNC, colnames))
#names(sequences) <- paste0("asv_", 1:length(sequences))


###Taxonomic assignment BLAST
#if(doTax){
#  library(taxonomizr)
#  library(taxize)
  
#  Biostrings::writeXStringSet(DNAStringSet(unlist(sequences)),
#                              "/SAN/Victors_playground/Eimeria_microbiome/Microbiome_seq.fasta")
  
 # clusters <- plotAmpliconNumbers(MA) 
  
  ###BLAST
  ## blastn -negative_gilist /SAN/db/blastdb/uncultured.gi -query /SAN/Victors_playground/Eimeria_microbiome/Microbiome_seq.fasta -db /SAN/db/blastdb/nt/nt -outfmt 11 -evalue 1e-5 -num_threads 10 -out /SAN/Victors_playground/Eimeria_microbiome/asv_vs_nt_microbiome.asn
  
  ## blast_formatter -archive /SAN/Victors_playground/Eimeria_microbiome/asv_vs_nt_microbiome.asn -outfmt "10 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid" > /SAN/Victors_playground/Eimeria_microbiome/asv_vs_nt_microbiome.blttax
  
  ###Read blast result 
  ## we read that ouput into R blast <-
#  blast <- read.csv("/SAN/Victors_playground/Eimeria_microbiome/asv_vs_nt_microbiome.blttax", header=FALSE)
  
#  names(blast) <- c("query", "subject", "pident", "length", "mismatch",
 #                   "gapopen", "qstart", "qend", "sstart", "send", "evalue",
 #                   "bitscore", "staxid")
#  blast <- as.data.table(blast)
#  blast$staxid <- as.character(blast$staxid)
  
#  read.nodes.sql("/SAN/db/taxonomy/nodes.dmp",
#                 "/SAN/db/taxonomy/taxonomizr.sql")
#  read.names.sql("/SAN/db/taxonomy/names.dmp",
#                 "/SAN/db/taxonomy/taxonomizr.sql")
  
#  blast.tax <- getTaxonomy(unique(blast$staxid),
#                           "/SAN/db/taxonomy/taxonomizr.sql")
  
#  blast.tax <- as.data.table(blast.tax, keep.rownames="staxid")
#  blast.tax$staxid <- gsub("\\s*", "", blast.tax$staxid)
  
#  blt <- merge(blast, blast.tax, by="staxid", all=TRUE)
  
  ## ## ## We need to be more clever if we want to use multiple
  ## ## ## hsps, this does not work for whole genome subjects eg.
  ### blt <- blt[,.(bitsum=sum(bitscore),
  ###              superkingdom, phylum, class, order, family, genus, species),
  ###           by=c("query", "subject")]
  
  ###    blt <- unique(blt)
  
#  blt <- blt[,.(bitdiff= bitscore - max(bitscore),
#                superkingdom, phylum, class, order, family, genus, species),
#             by=c("query")]
  
 # get.unique.or.na <- function (x){
    ## unique taxa at that level excluding potential NA's 
#    ux <- unique(as.character(x[!is.na(x)]))
    ## but return NA if they are not unique
#    if(length(ux)==1){return(ux)} else {as.character(NA)}
#  }
  
#  genus <- blt[bitdiff>-2, .(genus=get.unique.or.na(genus)),
#               by=query]
  
#  family <- blt[bitdiff>-7, .(family=get.unique.or.na(family)),
#                by=query]
  
#  order <- blt[bitdiff>-12, .(order=get.unique.or.na(order)),
#               by=query]
  
#  class <- blt[bitdiff>-20, .(class=get.unique.or.na(class)),
#               by=query]
  
#  phylum <- blt[bitdiff>-30, .(phylum=get.unique.or.na(phylum)),
#                by=query]
  
#  superkingdom <- blt[bitdiff>-50, .(superkingdom=get.unique.or.na(superkingdom)),
#                      by=query]
  
#  annot <- cbind(superkingdom[,c("query", "superkingdom")],
#                 phylum[,"phylum"],
#                 class[,"class"],
#                 order[,"order"],
#                 family[,"family"],
#                 genus[,"genus"])
#  
#  seqnametab <- as.data.table(cbind(query=names(sequences), sequences))
#  seqnametab <- merge(seqnametab, annot)
  
#  dupseq <- seqnametab$sequences[duplicated(seqnametab$sequences)]
  
#  seqnametab <- seqnametab[!duplicated(seqnametab$sequences),]
  
#  annot.list <- lapply(STNC, function (x) {
#    setkey(seqnametab, sequences)
#    seqnametab[colnames(x),
#               c("superkingdom", "phylum", "class", "order", "family", "genus")]
#  })
  
#  saveRDS(annot.list, file="/SAN/Victors_playground/Eimeria_microbiome/Microbiome_blast_tax.Rds")
#} else{
#  annot.list <- readRDS(file="/SAN/Victors_playground/Eimeria_microbiome/Microbiome_blast_tax.Rds")
#}


## ## Not needed anymore
## keep <- unlist(lapply(annot.list, nrow))>0
## annot.list <- annot.list[keep]
## STNC <- STNC[keep]


## name the annotation lists to have the names of the taxa 
#annot.list <- lapply(seq_along(annot.list), function (i){
#  an <- as.matrix(annot.list[[i]])
 # rownames(an) <- colnames(STNC[[i]])
#  an
#})

#names(annot.list) <- names(STNC)

#phylalist <- lapply(annot.list, function (x) {
 # if(nrow(x)>0){
 #   table(x[, "phylum"])
#  }
#})


#tabulate.taxa <- function(taxtab, taxon, phylumsubset){
 # if(nrow(taxtab)>0){
#    t <- taxtab[taxtab[, "phylum"]%in%phylumsubset, ]
 #   if(!is.null(ncol(t))){
#      table(t[, taxon])
#    } else {NULL} 
#  }else {NULL} 
#}


## Tabulate by specific phylum
#lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Cestoda"))
#lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Nematoda"))
#lapply(annot.list, function (x) tabulate.taxa(x,  "family", "Apicomplexa"))
#lapply(annot.list, function (x) tabulate.taxa(x, "genus",  "Platyhelminthes"))
#lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Streptophyta"))
#lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Ascomycota"))
#lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Chordata"))
#lapply(annot.list, function (x) tabulate.taxa(x, "phylum", "Ascomycota"))
#lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Bacteria"))



###Merge sample information
#PS.l <- lapply(seq_along(STNC)[keep], function(i){
#  phyloseq(otu_table(STNC[[i]], taxa_are_rows=FALSE),
#           sample_data(sample.data[rownames(STNC[[i]]), ]),
#           tax_table(annot.list[[i]]))
#})

#sumSeqByTax <- function (Phy, tax) {
#  counts <- data.frame(cbind(asvCount=colSums(otu_table(Phy)), tax_table(Phy)))
#  counts$asvCount <- as.numeric(as.character(counts$asvCount))
#  tapply(counts$asvCount, counts[, tax], sum)
#}

#readNumByPhylum <- lapply(PS.l, sumSeqByTax, "phylum")
#names(readNumByPhylum) <- names(STNC)[keep]


#readNumByGenus <- lapply(PS.l, sumSeqByTax, "genus") ## Change "text" in order to get counts per a different taxonomic level
#names(readNumByGenus) <- names(STNC)[keep]


#readNumByfamily <- lapply(PS.l, sumSeqByTax, "family") ## Change "text" in order to get counts per a different taxonomic level
#names(readNumByfamily) <- names(STNC)[keep]


####
#fill <- fillSampleTables(MA)
#MA@sequenceTableFilled <- fill@sequenceTableFilled


## Analyse all at once for now
#ALL <- Reduce(cbind, fill@sequenceTableFilled[keep])

## Problem: over all amplicons some ASVs are identical...
#table(duplicated(colnames(ALL)))

## sum up same reads over amplicons
#ALL.u <- do.call(rbind, by(t(ALL), rownames(t(ALL)), colSums))

## same for tax
#all.tax <- Reduce(rbind, annot.list[rownames(MA)[keep]])
#all.tax <- all.tax[rownames(ALL.u), ]


#PS <- phyloseq(otu_table(ALL.u, taxa_are_rows=TRUE),
#               sample_data(sample.data[rownames(ALL), ]),
#               tax_table(all.tax))

#prune_both_zero <- function (ps) {
#  p <- prune_samples(sample_sums(ps) > 0 , ps)
#  prune_taxa(taxa_sums(p) > 0 , p)
#}

#PS <- prune_both_zero(PS)
#PS.l <- lapply(PS.l, prune_both_zero)

################# ## HOW TO GO ON FROM HERE ## ######################
#### PS is now a single Phyloseq object over all amplicons. 

## For Phyloseq see: https://joey711.github.io/phyloseq/tutorials-index.html


