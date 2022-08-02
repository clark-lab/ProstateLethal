#Script for preparing MBPS data from MethPanel ready for visualisation and survival analysis ######

#set user defined working directory and create new folders to store data and output results
BASEDIR<-*user_defined_file_path*
dir.create("MBPS")
setwd(paste(BASEDIR,"/MBPS",sep=""))
dir.create("data")
dir.create("results")
setwd(BASEDIR)

#Multiplex data was generated from processing multiplex sequencing reads through methpanel https://github.com/thinhong/MethPanel using the config file in "05_01_ProcessMBPSdata.sh"

#Processed data is in the format of a .tsv file called a "bigtable" and has been copied from ".../MBPS/bigtable/")to a new folder in this github repository (".../MBPS/data/"). There were 2 bigtables for this project as the multiplex data was run in 2 sequencing batches.

#Each bigtable contains methylation percent values (".ratio") and coverage (".cov") for each CpG site and each sample.

#Load processed multiplex data
R1<-read.table("MBPS/data/bigTable1.tsv",sep="\t",header=T,na.strings = c("NA","."))
R2<-read.table("MBPS/data/bigTable2.tsv",sep="\t",header=T,na.strings = c("NA","."))

#minor name change
names(R1) <- gsub(x = names(R1), pattern = "X", replacement = "S") 
names(R2) <- gsub(x = names(R2), pattern = "X", replacement = "S")

#create vector specifying which columns contain methylation percent (".ratio") or coverage (".cov") for dataset 1
cov<-grep("cov", colnames(R1))
ratio<-grep("ratio", colnames(R1))

#extract gene names at each CpG site
genes<-sapply(R1$Amplicon,function(x){strsplit(as.character(x),"::")[[1]][1]})
genes<- gsub(x = genes, pattern = "lethal_", replacement="")

#extract unique gene names 
namelist<-unique(genes)

#extract genomic position of each CpG site
R_pos<-vector() 
for(i in 1:nrow(R1)){
  R_pos[i]<-strsplit(as.character(R1$Amplicon[i]),":")[[1]][4]
  R_pos[i]<-as.numeric(strsplit(as.character(R_pos[i]),"-")[[1]][1])+R1$Cpos[i]
}

#create vector specifying which columns contain methylation percent (".ratio") or coverage (".cov") for dataset 2
cov1<-grep("cov", colnames(R2))
ratio1<-grep("ratio", colnames(R2))

#extract gene names at each CpG site
genes1<-sapply(R2$Amplicon,function(x){strsplit(as.character(x),"::")[[1]][1]})
genes1<- gsub(x = genes1, pattern = "lethal_", replacement="")

#extract unique gene names 
namelist1<-unique(genes1)

#genes1 and namelist1 should be exactly the same as genes and namelist from run 1 above 
genes1 == genes #TRUE
namelist1 == namelist #TRUE

#extract genomic position of each CpG site
R2_pos<-vector() 
for(i in 1:nrow(R2)){
  R2_pos[i]<-strsplit(as.character(R2$Amplicon[i]),":")[[1]][4]
  R2_pos[i]<-as.numeric(strsplit(as.character(R2_pos[i]),"-")[[1]][1])+R2$Cpos[i]
}

R2_pos == R_pos #TRUE

#Extract methylation data for Run 1

#For individual CpGs, exclude if coverage <100
R1_cov<-R1[,cov]
R1_ratio<-R1[,ratio]

a<-colnames(R1_cov)
b<-colnames(R1_ratio)

a_name<-sapply(a,function(x){strsplit(x,"[.]")[[1]][1]})
b_name<-sapply(b,function(x){strsplit(x,"[.]")[[1]][1]})

R1_cov_ord<-R1_cov[,match(b_name,a_name)] 
R1_meth<-R1_ratio
R1_meth[(R1_cov_ord<100)]<-"NA" 

#Extract methylation data for Run 2

#For individual CpGs, exclude if coverage <100
R2_cov<-R2[,cov]
R2_ratio<-R2[,ratio]

a<-colnames(R2_cov)
b<-colnames(R2_ratio)

a_name<-sapply(a,function(x){strsplit(x,"[.]")[[1]][1]})
b_name<-sapply(b,function(x){strsplit(x,"[.]")[[1]][1]})

R2_cov_ord<-R2_cov[,match(b_name,a_name)] 
R2_meth<-R2_ratio
R2_meth[(R2_cov_ord<100)]<-"NA" 

#Combine the two runs together
M<-cbind(R1_meth,R2_meth)

##extract coverage data 
# run 1
meancovtab<-matrix(data=NA,nrow=length(namelist),ncol=length(cov))
for(i in 1:length(namelist)){
  
  temp<-namelist[i]
  gen<-which(genes==temp)
  len<-length(gen) 
  
  if(len>1){
    meancovtab[i,]<-as.numeric(colMeans(R1[gen,cov]))
  } else {
    meancovtab[i,]<-as.numeric(R1[gen,cov])
  }
}
#Adding column and row names
colnames(meancovtab)<-colnames(R1[,cov])
rownames(meancovtab)<-namelist

#run 2
meancovtab1<-matrix(data=NA,nrow=length(namelist),ncol=length(cov1))
for(i in 1:length(namelist1)){
  
  temp<-namelist1[i]
  gen<-which(genes1==temp)
  len<-length(gen) 
  
  if(len>1){
    meancovtab1[i,]<-as.numeric(colMeans(R2[gen,cov1]))
  } else {
    meancovtab1[i,]<-as.numeric(R2[gen,cov1])
  }
}
#Adding column and row names
colnames(meancovtab1)<-colnames(R2[,cov1])
rownames(meancovtab1)<-namelist1

#save out prepared methylation data and gene annotation
save(M, genes, namelist, R_pos, meancovtab, meancovtab1, file="MBPS/results/ProcessedMultiplex.RData")

#remove the word "ratio" from extract columns in M that contain methylation percent values (".ratio")
colnames(M)<-strsplit(colnames(M),".ratio")

#read in validation cohort clinical information (Available upon request. Download into metadata folder)
clin<-read.csv("metadata/Validation_Clinical.csv",header=T)

#re-order methyation matrix M to match order of samples in sample table, and re-name to M_samp 
M_samp<-M[,match(clin$Sample_ID,colnames(M))]

#convert M_samp into a dataframe for analysis
M_samp<-as.data.frame(sapply(M_samp,as.numeric))

#calculate the mean methylation across each amplicon for each sample
avM_samp<-matrix(data=NA,nrow=length(namelist),ncol=length(colnames(M_samp)))
for(i in 1:length(namelist)) {
  temp<-namelist[i]
  gen<-which(genes==temp)
  len<-length(gen)
  
  if(len>1) {
    avM_samp[i,]<-as.numeric(colMeans((M_samp[gen,])))
  }
  else {
    avM_samp[i,]<-as.numeric(M_samp[gen,])
  }
}
#Adding column and row names
colnames(avM_samp)<-colnames(M_samp)
rownames(avM_samp)<-namelist

#save out methylation data averaged by amplicon and ordered by sample 
save(avM_samp, file="MBPS/results/ProcessedMultiplexMeans.RData")


