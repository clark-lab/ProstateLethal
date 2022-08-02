#Script for preparing downsampled CACNA2D4 MBPS data for visualisation and survival analysis

## Load processed downsampled multiplex data (generated from processing multiplex sequencing reads through methpanel https://github.com/thinhong/MethPanel)

# Processed data for each downsampled level is exported from methpanel as a .tsv file called a "bigtable". There were 2 bigtables for each level for this project as the multiplex data was run in 2 sequencing batches.

# Bigtables contain methylation percent values (".ratio") and coverage (".cov") for each CpG site and each sample.

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

options(scipen=999)

#Load data
R1_list<-R2_list<-list()

R1_list[[1]] <- read.table("MBPS/data/bigTable1_fullRead_1000.tsv",sep="\t",header=T,na.strings = c("NA","."))
R1_list[[2]] <- read.table("MBPS/data/bigTable1_fullRead_10000.tsv",sep="\t",header=T,na.strings = c("NA","."))
R1_list[[3]] <- read.table("MBPS/data/bigTable1_fullRead_50000.tsv",sep="\t",header=T,na.strings = c("NA","."))
R1_list[[4]] <- read.table("MBPS/data/bigTable1_fullRead_100000.tsv",sep="\t",header=T,na.strings = c("NA","."))
names(R1_list)<-c("cov_1000","cov_10000","cov_50000","cov_100000")

R2_list[[1]] <- read.table("MBPS/data/bigTable2_fullRead_1000.tsv",sep="\t",header=T,na.strings = c("NA","."))
R2_list[[2]] <- read.table("MBPS/data/bigTable2_fullRead_10000.tsv",sep="\t",header=T,na.strings = c("NA","."))
R2_list[[3]] <- read.table("MBPS/data/bigTable2_fullRead_50000.tsv",sep="\t",header=T,na.strings = c("NA","."))
R2_list[[4]] <- read.table("MBPS/data/bigTable2_fullRead_100000.tsv",sep="\t",header=T,na.strings = c("NA","."))
names(R1_list)<-c("cov_1000","cov_10000","cov_50000","cov_100000")

#minor name change to all tables
namefunc<-function(k){names(k)<-gsub(x = names(k), pattern = "X", replacement = "S")
return(k)} 

R1_listN<-lapply(R1_list,namefunc)
R2_listN<-lapply(R2_list,namefunc)

#create function for extracting ratio data and setting to NA if coverage <100
extractfunc<-function(k){
  R_cov<-k[,grep("cov", colnames(k))]
  R_ratio<-k[,grep("ratio", colnames(k))]
  a_name<-sapply(colnames(R_cov),function(x){strsplit(x,"[.]")[[1]][1]})
  b_name<-sapply(colnames(R_ratio),function(x){strsplit(x,"[.]")[[1]][1]})
  R_meth<-R_ratio
  R_cov_ord<-R_cov[,match(b_name,a_name)] 
  R_meth[(R_cov_ord<100)]<-"NA" 
  return(list(R_meth,R_cov_ord))
}

Meth1_list<-lapply(R1_listN,extractfunc)
Meth2_list<-lapply(R2_listN,extractfunc)

#Combine the two runs together
if(identical(names(Meth1_list),names(Meth2_list))){
  M<-Cov<-list()
  for(i in 1:length(Meth1_list)){
    M[[i]]<-cbind(Meth1_list[[i]][[1]],Meth2_list[[i]][[1]])
    Cov[[i]]<-cbind(Meth1_list[[i]][[2]],Meth2_list[[i]][[2]])
  }
}

names(M)<-names(Cov)<-names(Meth1_list)

#save out methylation data averaged by amplicon and ordered by sample 
save(M2, Cov2, M_tab,Cov_tab, file="MBPS/ProcessedData/ProcessedDownsampledMultiplex.RData")

