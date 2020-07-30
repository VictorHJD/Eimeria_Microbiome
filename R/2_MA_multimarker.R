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
## But remember: if you change the MultiAmplicon Analysis, the
## taxonomic annotation might be out of sync...

###################Full run Microbiome#######################
#Preparation of files

##These are the same steps that are followed by the DADA2 pipeline

#path <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/2018_22_Eie_TestRun/" ## Test run 24.06.2020
path <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/2018_22_Eie_FullRun_1/" ## Full run 29.06.2020
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 


samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
samples<- gsub("S\\d+-", "\\1", basename(samples))
samples<-gsub("-", "_", basename(samples))

#Extra step in the pipeline: quality plots of the reads 
## plotQualityProfile(fastqF[[1]])
## plotQualityProfile(fastqF[[2]])
## plotQualityProfile(fastqR[[1]])
## plotQualityProfile(fastqR[[2]])

#Creation of a folder for filtrated reads 

#filt_path <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/filtered_TestRun/"
filt_path <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/filtered_FullRun_1/"

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
ptable <- read.csv(file = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/primer.file.multi.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "TS.SequenceF"]
primerR <- ptable[, "TS.SequenceR"]
names(primerF) <- as.character(ptable[, "corrected.NameF"])
names(primerR) <- as.character(ptable[, "corrected.NameR"])

primer <- PrimerPairsSet(primerF, primerR)

##Multi amplicon pipeline
if(doMultiAmp){
  MA <- MultiAmplicon(primer, files)
  #filedir <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/Stratified_files_TestRun"
  filedir <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/Stratified_files_FullRun_1"
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
  
  #saveRDS(MA, "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MA_Multi_TestRun.RDS")
  saveRDS(MA, "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MA_Multi_FullRun_1.RDS")
} else{
  #MA <- readRDS("/SAN/Victors_playground/Eimeria_microbiome/MA_Multi_TestRun.RDS")
  MA <- readRDS("/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MA_Multi_FullRun_1.RDS")
}

trackingF <- getPipelineSummary(MA) 
## plotPipelineSummary(trackingF) 
## plotPipelineSummary(trackingF) + scale_y_log10()


plotAmpliconNumbers(MA)

###New taxonomic assignment 
if(doTax){
MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    #infasta = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/in_TestRun.fasta",
                    infasta = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/in_FullRun_1.fasta",
                    #outblast = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/out_TestRun.fasta",
                    outblast = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/out_FullRun_1.fasta",
                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                    num_threads = 20) ##Change for use more power!!

#saveRDS(MA, file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MATax_TestRun.Rds")
saveRDS(MA, file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MATax_FullRun_1.Rds")
}else{
  MA1 <- readRDS("/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MATax_FullRun_1.Rds")
  MA2 <- readRDS("/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MATax_TestRun.Rds")  
}
### Add sample information
if(!exists("sample.data")){
  source("~/GitProjects/Eimeria_Microbiome/R/1_Data_preparation.R")
}

MA1 <- addSampleData(MA1, sample.data)
MA2 <- addSampleData(MA2, sample.data)

#saveRDS(MA2, file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MASample_TestRun.Rds")  ###START from here now! 
#saveRDS(MA1, file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MASample_FullRun_1.Rds")

####Raw counts
temp1 <- data.frame(colSums(getRawCounts(MA1)))
temp1[,2] <- rownames(temp1)
colnames(temp1) <- c("Raw_counts", "labels")
rownames(temp1) <- c(1:nrow(temp1))
temp1 <- data.frame(labels = temp1$labels, Raw_counts_1 = temp1$Raw_counts)

temp2 <- data.frame(colSums(getRawCounts(MA2)))
temp2[,2] <- rownames(temp2)
colnames(temp2) <- c("Raw_counts", "labels")
rownames(temp2) <- c(1:nrow(temp2))
temp2 <- data.frame(labels = temp2$labels, Raw_counts_2 = temp2$Raw_counts)

rawcounts<- join(temp1, temp2, by= "labels")
rm(temp1, temp2)

rawcounts%>%
  dplyr::mutate(Raw_counts= Raw_counts_1 + Raw_counts_2)->rawcounts

hist(rawcounts$Raw_counts)
summary(rawcounts$Raw_counts)
boxplot(rawcounts$Raw_counts)

sum(rawcounts$Raw_counts) ###Total raw counts reads = 5,360,843

##To phyloseq
if(doPhyloseq){
##Sample data
PS1 <- toPhyloseq(MA1, colnames(MA1))
sum(otu_table(PS1)) ##Total denoised reads = 3,610,637 for FullRun

PS2<- toPhyloseq(MA2, colnames(MA2))
sum(otu_table(PS2)) ##Total denoised reads = 267,293 for TestRun

##Primer data
## just sorting out primers whithout any taxannot
MA1 <- MA1[which( !unlist(lapply(MA1@taxonTable, is.null))), ] ##Make the next function work 
PS1.l <- toPhyloseq(MA1, colnames(MA1),  multi2Single=FALSE) 
#saveRDS(PS1.l, file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/PhyloSeqList_FullRun_1.Rds")

MA2 <- MA2[which( !unlist(lapply(MA2@taxonTable, is.null))), ] ##Make the next function work 
PS2.l <- toPhyloseq(MA2, colnames(MA2),  multi2Single=FALSE) 
#saveRDS(PS1.l, file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/PhyloSeqList_TestRun.Rds")

along<- names(PS2.l) ## Run with less primers working
PS.l <- lapply(along, function(i) merge_phyloseq(PS1.l[[i]], PS2.l[[i]])) ##Merge all the information from both experiments
names(PS.l) <- names(PS2.l) ###Use the names from test list 

#saveRDS(PS.l, file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/PhyloSeqList_All.Rds") ###For primer analysis
#saveRDS(PS2, file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/PhyloSeqData_TestRun.Rds") ###For Sample analysis (Susana and Victor)
#saveRDS(PS1, file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/PhyloSeqData_FullRun_1.Rds") ###For Sample analysis (Susana and Victor)
rm(along,PS1.l, PS2.l)
}else{
####Merge phyloseq objects 
  
  PS2 <-  readRDS(file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/PhyloSeqData_TestRun.Rds")
  PS1<-  readRDS(file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/PhyloSeqData_FullRun_1.Rds")
  
  PS <- merge_phyloseq(PS1, PS2) ###Works! 
  #saveRDS(PS, file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/PhyloSeqData_All.Rds") ###Results from full + test run 
}

rm(MA1, MA2, PS1, PS2)