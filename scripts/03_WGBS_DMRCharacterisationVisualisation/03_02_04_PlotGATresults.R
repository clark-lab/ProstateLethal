#Script for plotting overlap between DMRs and LNCaP and PrEC ChromHMM segmentation data
#Generates Figure 1f and Figure S3

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#read in data of overlaps between DMRs and PrEC 15 states generated in '03_02_03_RunGAT_DMRChromHMM.py'
HyperDMR_GAT_PreC15<-read.table("WGBS/results/GAT_HyperDMR_PrEC15.tsv",skip=9,header=T)

HypoDMR_GAT_PreC15<-read.table("WGBS/results/GAT_HypoDMR_PrEC15.tsv",skip=9,header=T)

#Make into list
listall_PreC15<-list(HyperDMR_GAT_PreC15,HypoDMR_GAT_PreC15)
nameslist<-c("HyperDMR_GAT_PreC15","HypoDMR_GAT_PreC15")

#Extract number of overlaps between DMRs and different PrEC chromHMM states - put data into matrix
matchnames<-HyperDMR_GAT_PreC15$annotation

listall_ObsPrEC<-matrix(data=NA,ncol=length(listall_PreC15),nrow=15)
for(i in 1:length(listall_PreC15)){
  x<-listall_PreC15[[i]]
  x2<-x[match(matchnames,x$annotation),]
  listall_ObsPrEC[,i]<-x2$overlap_nsegments
  
  labels<-x2$annotation
  labels<-gsub("_"," ",labels)
}

rownames(listall_ObsPrEC)<-labels
colnames(listall_ObsPrEC)<-c(nameslist)


#read in data of overlaps between DMRs and LNCaP 15 states generated in '03_02_03_RunGAT_DMRChromHMM.py'
HyperDMR_GAT_LNCaP15<-read.table("WGBS/results/GAT_HyperDMR_LNCaP15.tsv",skip=9,header=T)

HypoDMR_GAT_LNCaP15<-read.table("WGBS/results/GAT_HypoDMR_LNCaP15.tsv",skip=9,header=T)

#Make into list
listall_LNCaP15<-list(HyperDMR_GAT_LNCaP15,HypoDMR_GAT_LNCaP15)
nameslist<-c("HyperDMR_GAT_LNCaP15","HypoDMR_GAT_LNCaP15")

matchnames<-HypoDMR_GAT_LNCaP15$annotation

#Extract number of overlaps between DMRs and different LNCaP chromHMM states - put data into matrix
listall_ObsLNCaP<-matrix(data=NA,ncol=length(listall_LNCaP15),nrow=15)
for(i in 1:length(listall_LNCaP15)){
  x<-listall_LNCaP15[[i]]
  x2<-x[match(matchnames,x$annotation),]
  listall_ObsLNCaP[,i]<-x2$overlap_nsegments
  
  labels<-x2$annotation
  labels<-gsub("_"," ",labels)
}

rownames(listall_ObsLNCaP)<-labels
colnames(listall_ObsLNCaP)<-c(nameslist)

#reorder and rename ChromHMM overlap results based on the 15 states
E1<-"Active TSS"
E2<-"Flanking active TSS"
E3<-"Transcription at gene5and3"
E4<-"Strong transcription"
E5<-"Weak transcription"
E6<-"Genic Enhancer"
E7<-"Enhancer"
E8<-"ZNF genes and repeats"
E9<-"Heterochromatin"
E10<-"Bivalent/Poised TSS"
E11<-"Flanking Bivalent TSS/Enh"
E12<-"Bivalent enhancer"
E13<-"Repressed polycomb"
E14<-"Weak repressed polycomb"
E15<-"Quiescent"

Evec<-c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12","E13","E14","E15")
Enames<-c(E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13,E14,E15)

listall_ObsLNCaPORD<-listall_ObsLNCaP[match(Evec,rownames(listall_ObsLNCaP)),]
listall_ObsLNCaPORDname<-data.frame(listall_ObsLNCaPORD,Enames)

listall_ObsPrECORD<-listall_ObsPrEC[match(Evec,rownames(listall_ObsPrEC)),]
listall_ObsPrECORDname<-data.frame(listall_ObsPrECORD,Enames)


#calculate percentage of observed overlap
percHyperOL_L<-(listall_ObsLNCaPORDname[,1]/sum(listall_ObsLNCaPORDname[,1])*100)
percHypoOL_L<-(listall_ObsLNCaPORDname[,2]/sum(listall_ObsLNCaPORDname[,2])*100)

percHyperOL_P<-(listall_ObsPrECORDname[,1]/sum(listall_ObsPrECORDname[,1])*100)
percHypoOL_P<-(listall_ObsPrECORDname[,2]/sum(listall_ObsPrECORDname[,2])*100)

pdf("WGBS/plots/HyperDMR_overlap_ChromHMM_PrEC_LNCaP_difference.pdf",width=10)
par(mar=c(12,5,2,2))
bp<-barplot((percHyperOL_P-percHyperOL_L),main="Hyper DMRs - % overlapping segments",las=2,ylim=c(-30,30),names.arg=listall_ObsLNCaPORDname[,3],las=2,ylab="% overlap PrEC minus % overlap LNCaP ")
dev.off()

pdf("WGBS/plots/HypoDMR_overlap_ChromHMM_PrEC_LNCaP_difference.pdf",width=10)
par(mar=c(12,5,2,2))
bp<-barplot((percHypoOL_P-percHypoOL_L),main="Hypo DMRs - % overlapping segments",las=2,ylim=c(-20,20),names.arg=listall_ObsLNCaPORDname[,3],las=2,ylab="% overlap PrEC minus % overlap LNCaP ")
dev.off()
