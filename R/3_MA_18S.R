## Please uncomment the first time you run this and re-install packages
##Sequence cleaning and Multiamplicon pipeline for Access Array of Eimeria infection experiment
## require(devtools)
## devtools::install_github("derele/MultiAmplicon", force= T)
## devtools::install_github("derele/dada2", force= T)
library("lifecycle", lib.loc="/usr/local/lib/R/site-library") 
library("ggplot2")
library("MultiAmplicon")
library("reshape")
library("phyloseq")
library("data.table")
library("taxonomizr")
library("taxize")
library("parallel")

## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doFilter <- FALSE

doMultiAmp <- FALSE

doTax <- FALSE

doPhyloseq<- FALSE

doSample<- FALSE

doMicrobiome<- FALSE
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

#trackingF <- getPipelineSummary(MA) 
## plotPipelineSummary(trackingF) 
## plotPipelineSummary(trackingF) + scale_y_log10()

###New taxonomic assignment 
if(doTax){
MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    infasta = "/SAN/Victors_playground/Eimeria_microbiome/in.fasta",
                    outblast = "/SAN/Victors_playground/Eimeria_microbiome/out.fasta",
                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                    num_threads = 20)

saveRDS(MA, file="/SAN/Victors_playground/Eimeria_microbiome/MATax.Rds")
}else{
  MA <- readRDS("/SAN/Victors_playground/Eimeria_microbiome/MATax.Rds")
}

### Add sample information
if(doSample){
if(!exists("sample.data")){
  source("~/GitProjects/Eimeria_Microbiome/R/1_Data_preparation.R")
}

MA <- addSampleData(MA, sample.data)

saveRDS(MA, file="/SAN/Victors_playground/Eimeria_microbiome/MASample.Rds")  ###START from here now! 
}else{
  MA<- readRDS("/SAN/Victors_playground/Eimeria_microbiome/MASample.Rds")
}
####Raw counts 18S
rawcounts18S <- as.data.frame(MA@rawCounts)
rawcounts18S <- as.data.frame(t(rawcounts18S))
rawcounts18S[,1]<- NULL ###16S V1_V1
rawcounts18S[,2]<-NULL ###16S V3_V4
rawcounts18S[,2] <- rownames(rawcounts18S)
colnames(rawcounts18S) <- c("Raw_counts_18S", "labels")
rownames(rawcounts18S) <- c(1:nrow(rawcounts18S))
rawcounts18S <- data.frame(labels = rawcounts18S$labels, Raw_counts_18S = rawcounts18S$Raw_counts_18S) 

hist(rawcounts18S$Raw_counts_18S)
summary(rawcounts18S$Raw_counts_18S)
boxplot(rawcounts18S$Raw_counts_18S)

sum(rawcounts18S$Raw_counts_18S) ###Total raw counts reads 18S= 1,116,513

##To phyloseq
##Sample data
if(doMicrobiome){
PS <- toPhyloseq(MA, colnames(MA))
sum(otu_table(PS)) ##Total denoised reads = 2015804
saveRDS(PS, file="/SAN/Victors_playground/Eimeria_microbiome/PhyloSeqData.Rds") ###For Sample analysis (Susana and Victor)
}

##Primer data
PS.l <- toPhyloseq(MA, colnames(MA),  multi2Single=FALSE)

#saveRDS(PS.l, file="/SAN/Victors_playground/Eimeria_microbiome/PhyloSeqList.Rds") ###For Microbiome analysis (Victor)

PS.18S <- phyloseq(otu_table(PS.l$wang1141_13_F.Nem_0425_6_3_R), 
                   sample_data(PS.l$wang1141_13_F.Nem_0425_6_3_R), 
                   tax_table(PS.l$wang1141_13_F.Nem_0425_6_3_R))

sum(otu_table(PS.18S)) ##Total denoised reads = 853,134
#saveRDS(PS.18S, file="/SAN/Victors_playground/Eimeria_microbiome/PS_18S.Rds") ###Information from 18S

rm(MA, PS.l)