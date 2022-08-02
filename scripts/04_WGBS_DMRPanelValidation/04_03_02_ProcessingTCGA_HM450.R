#Script for reading in raw HM450 data (idat files), followed by processing and basic QC

#Generates processed dataset 

#load libraries
library(minfi)
library(wateRmelon)

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#User will need to download all TCGA PRAD level 1 idats listed in targetsboth.csv from the NCI GDC Data portal into the local data directory ('WGBS/data/TCGA_Level1/'). This can be done manually, or using a function such as readIDATDNAmethylation in the TCGAbiolinks package - https://rdrr.io/bioc/TCGAbiolinks

#read in sample info table generated in '04_03_01_SelectingSamplesTCGA_HM450'
targetsboth<-read.csv("WGBS/data/targetsboth.csv")

#read in raw idats
RGsetExt<-read.450k.exp(targets=targetsboth,extended=TRUE)

#perform basic QC
MSet<-preprocessRaw(RGsetExt)
qc <- getQC(MSet)

#Filter out bad samples 
gset.quantile <- preprocessQuantile(RGsetExt,removeBadSamples = TRUE, badSampleCutoff = 10.5)
dim(gset.quantile)
#[1] 485512 437
#note, 4 samples removed due to poor QC

RGsetExt2<-RGsetExt[,sampleNames(gset.quantile)]

#Perform functional normalisation
RGsetExt_fun<-preprocessFunnorm(RGsetExt2, nPCs=2, sex = NULL, bgCorr = TRUE,
                                dyeCorr = TRUE, verbose = TRUE)

#identify poorly performing probes
pn<-detectionP(RGsetExt2, type = "m+u")

newpfilter<-function(x){length(which(x>0.01))}
failedprobes<-apply(pn,1,newpfilter)

#create vector for probes to be removed that have failed in at least 10% of patients
tenperc<-(dim(RGsetExt_fun)[2]/100)*10
toremove<-which(failedprobes>tenperc)

length(toremove)
# [1] 2636

#identify probes with a low beadcount (likely to give poor data)
bc<-beadcount(RGsetExt2)
newbc<-function(x){length(which(is.na(x)==TRUE))}
failedbeads<-apply(bc,1,newbc)

#create vector for probes to be removed with low beadcount in at least 10% of patients
toremove2<-(which(failedbeads>tenperc))
length(toremove2)
#[1] 95

#save processed RGset
save(toremove,toremove2,RGsetExt_fun,file="WGBS/data/RGsetExt_fun.RData")
