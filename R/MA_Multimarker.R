## Please uncomment the first time you run this and re-install packages
##Sequence cleaning and Multiamplicon pipeline for Access Array of Eimeria infection experiment
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

path <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/2018_22_Eie_TestRun/" ## Test run 24.06.2020
#path <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/2018_22_Eie_FullRun_1/" ## Full run 
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

filt_path <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/filtered_TestRun/"

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
  filedir <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/Stratified_files_TestRun"
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
  
  saveRDS(MA, "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MA_Multi_TestRun.RDS")
} else{
  MA <- readRDS("/SAN/Victors_playground/Eimeria_microbiome/MA_Multi_TestRun.RDS")
}

trackingF <- getPipelineSummary(MA) 
## plotPipelineSummary(trackingF) 
## plotPipelineSummary(trackingF) + scale_y_log10()


plotAmpliconNumbers(MA)

###New taxonomic assignment 

MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    infasta = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/in_TestRun.fasta",
                    outblast = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/out_TetsRun.fasta",
                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                    num_threads = 12) ##Change for use more power!!

saveRDS(MA, file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MATax_TestRun.Rds")


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

library("lifecycle", lib.loc="/usr/local/lib/R/site-library")
require(dplyr)

##considering 30µL of elution volume 
sample.data %>%
  mutate(Total_DNA = (sample.data$Conc_DNA*30)*0.001) -> sample.data ### Add a new variable that will contain total DNA extracted per sample in µg

sample.data %>%
  mutate(Microbial_density = sample.data$Total_DNA/(sample.data$fecweight_DNA*1000)) -> sample.data ### Total DNA extracted per sample in µg by feces weight in mg

sample.data$labels<- as.vector(sample.data$labels)
rownames(sample.data) <- make.unique(sample.data$labels)
MA <- addSampleData(MA, sample.data)

saveRDS(MA, file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MASample_TestRun.Rds")  ###START from here now! 

####Raw counts
rawcounts <- data.frame(colSums(getRawCounts(MA)))
rawcounts[,2] <- rownames(rawcounts)
colnames(rawcounts) <- c("Raw_counts", "labels")
rownames(rawcounts) <- c(1:nrow(rawcounts))
rawcounts <- data.frame(labels = rawcounts$labels, Raw_counts = rawcounts$Raw_counts) 

hist(rawcounts$Raw_counts)
summary(rawcounts$Raw_counts)
boxplot(rawcounts$Raw_counts)

sum(rawcounts$Raw_counts) ###Total raw counts reads = 389165 for TestRun

##To phyloseq
##Sample data
PS <- toPhyloseq(MA, colnames(MA))
sum(otu_table(PS)) ##Total denoised reads = 267293 for TestRun
##Primer data
#PS.l <- toPhyloseq(MA, colnames(MA),  multi2Single=FALSE) Not working, check later!

#saveRDS(PS.l, file="/SAN/Victors_playground/Eimeria_microbiome/PhyloSeqList.Rds") ###For primer analysis
saveRDS(PS, file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/PhyloSeqData_TestRun.Rds") ###For Sample analysis (Susana and Victor)
