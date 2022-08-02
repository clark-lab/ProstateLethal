#Script for extracting DNA methylation at selected DMRs in blood using public data
#Generates data for use in Figure S4
#Generates data for use in Table S4
#Generates bedGraph files used to create Figure S1

#load library
library(GenomicRanges)

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#Download public blood methylation data from GEO GSE40279 into ~/WGBS/data/

#read corresponding sample data table
samples <- read.table("WGBS/data/GSE40279.txt"), header=TRUE, stringsAsFactors=FALSE, sep="\t")
load("WGBS/data/GSE40279_methyl.Rdata"))

#reduce beta value table to only male samples (as this is a prostate cancer study, therefore males only)
B <- B[, samples$id[samples$gender=="M"]]

ncol(B) #318


#Download Illumina HM450 hg19 annotation data from Illumina website into folder "WGBS/annotationdata/" (File: "HumanMethylation450 v1.2 Manifest File (CSV Format)" from https://sapac.support.illumina.com/downloads/infinium_humanmethylation450_product_files.html)

#read in above annotation data for 450K array
fTab<-read.csv("WGBS/annotationdata/humanmethylation450_15017482_v1-2.csv",header=T,skip=7)

#turn 450K annotation data into Genomic Ranges object
fTabGR<-GRanges(seqnames=paste("chr",as.character(fTab$CHR[1:485512]),sep=""), ranges=IRanges(start=fTab$MAPINFO[1:485512], width=1),data=fTab$IlmnID[1:485512])

#read in DMRs
DMRs <- read.table("metadata/SelectedRegions.bed"))

#convert DMRs into genomic ranges object and re-order in ascending chromosomal position
DMRGR<-GRanges(seqnames=DMRs[,1],ranges=IRanges(start=DMRs[,2],end=DMRs[,3]))
seqlevels(DMRGR)<-seqlevels(DMRGR)[c(1,8,10:14,2:7,9,15)]
DMRGR_ord<-sort(DMRGR)

#find overlaps between DMRs and 450K annotation data
ol<-findOverlaps(DMRGR_ord,fTabGR)

#order beta values to be the same as 
mat<-match(fTab$IlmnID,rownames(B))

B_ord<-B[mat,]

#for each sample average methylation across probes in a DMR
uniqdmr<-unique(ol@from)

#i<-1
av_dmr<-matrix(data=NA,ncol=ncol(B_ord),nrow=length(uniqdmr))
for(i in 1:length(uniqdmr)){
  temp<-which(ol@from==uniqdmr[i])
  if(length(temp)>1){
    av_dmr[i,]<-colMeans(B_ord[ol@to[temp],])
  } else {
    av_dmr[i,]<-B_ord[ol@to[temp],]
  }
}
templist<-list()
for(i in 1:length(uniqdmr)){
    templist[[i]]<-which(ol@from==uniqdmr[i])
}
#DMR18 not present in the data for some reason - likely because they have removed chromosome X from their analaysis

#insert an empty row where CACNA2D4 would be (as no probes cover this region on the 450K array)
av_dmr2<-rbind(av_dmr[1:9,],rep(NA,ncol(av_dmr)),av_dmr[10:17,])

rownames(av_dmr2)<-paste("DMR",1:18,sep="")

#write out average methylation at each DMR for each sample, data to be used for Figure S4
write.csv(av_dmr2,"WGBS/results/AvMeth_DMRs_maleblood.csv")

#summarise results and write out
mean<-rowMeans(av_dmr2)
sd<-apply(av_dmr2,1,sd)

summary<-data.frame(mean,sd)
rownames(summary)<-rownames(av_dmr2)  

#write out average methylation at each DMR, averaged across samples for Table S4
write.csv(summary,"WGBS/results/Summary_DMRs_maleblood.csv")

#####writing out bed file with 10,000bp from center
res<-resize(DMRGR,width=10000,fix="center",ignore.strand=FALSE)

#find overlaps between DMRs and 450K annotation data
ol_new<-findOverlaps(res,fTabGR)

B_CpG<-B_ord[unique(ol_new@to),]

avprobes<-rowMeans(B_CpG,na.rm=T)

bloodbed<-data.frame(seqnames=paste("chr",fTab$CHR[unique(ol_new@to)],sep=""),start=fTab$MAPINFO[unique(ol_new@to)]-1, end=fTab$MAPINFO[unique(ol_new@to)],scores=avprobes)

#write out bedgraph files for IGV plots in Figure S1
write.table(bloodbed,
             "WGBS/results/blood_bed.bedGraph",
             quote=F, sep="\t", row.names=F, col.names=F,na="0")
