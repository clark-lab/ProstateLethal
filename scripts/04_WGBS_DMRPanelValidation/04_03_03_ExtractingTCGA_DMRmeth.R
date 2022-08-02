#Script for removing noisy probes from TCGA PRAD HM450 data and extracting methylation data at DMRs, and comparing DMR methylation and expression of nearest protein-coding gene

#Generates data for Table S4
#Generates plot for Figure S4b
#Generates bedgraph files for IGV plots in Figure S1
#Generates plots for Figure S5
#Generates plot for Figure S11b

#load libraries
library(minfi)
library(GenomicRanges)

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#load processed HM450 TCGA PRAD data, generated in '04_03_02_ProcessingTCGA_HM450.R'
load("WGBS/data/RGsetExt_fun.RData")

#create vector of probes to be removed for high detection p-value or low beadcount
remove<-unique(c(names(toremove),names(toremove2)))

#read in csv file containing details of "noisy probes" that should be excluded from the analysis - Additional File 2 from http://www.biomedcentral.com/1471-2164/15/51
forfiltering<-read.csv("metadata/NoisyProbes.csv",row.names=1)

#create vector of probes that cross-hybridse to multiple locations or map to repeat sequences
discard<-which(forfiltering$MultiMap>0 | forfiltering$Repeat>0)
probediscard<-rownames(forfiltering)[discard]

discardall<-unique(c(remove,probediscard))

#extract methylation beta values
beta1<-getBeta(RGsetExt_fun)

#remove all suspicious probes
matchdiscard<-match(discardall,rownames(beta1))

beta2<-beta1[-matchdiscard,]

#dimensions of final beta value table used for analysis
dim(beta2)
#[1] 414133    437
rm(beta1)

#read in annotation data for 450K array (instructions for download in '04_02_ExtractingBlood_DMRmeth.R')
fTab<-read.csv("WGBS/annotationdata/humanmethylation450_15017482_v1-2.csv",header=T,skip=7)

#match rows in annotation data to those in the beta value table 
fTab3<-fTab[match(rownames(beta2),fTab$IlmnID),]

#read in genomic coordinates for DMRs selected for validation
OrigDMR<-read.table("metadata/SelectedRegions.bed")

#manually adding gene names
#note, some of these changed later as reviewer suggested to use alternative gene names
OrigDMR_names<-c("CALM1","EPHB3","CDH17","PAH","BHLHB9","CACNA2D4","PRDM8","LECT1","RP11-927P21.4","CD34","AC074091.13","TBX1","CCDC67","PARP6","ZNF655","EFCAB4B","CDO1","Y-RNA")

#Convert DMR coordinate tables to bed files
OrigDMR_GR<-GRanges(seqnames=OrigDMR[,1],IRanges(start=OrigDMR[,2],end=OrigDMR[,3]))

#convert TCGA annotation data to GRanges
fTab_ranges<-GRanges(seqnames=paste("chr",fTab3$CHR,sep=""),IRanges(start=fTab3$MAPINFO,end=fTab3$MAPINFO))


#read in sample info table generated in '04_03_01_SelectingSamplesTCGA_HM450'
targetsboth3<-read.csv("WGBS/data/targetsboth.csv",header=T)

##Create vectors of sample ID, and sample type (tumour or normal)
id<-targetsboth3$Sample_Name[match(colnames(beta2),as.character(targetsboth3$Array.Data.File))]

#function for tissue type info
sample_type<-function(x){
  substr(as.character(x),14,15)
}
type<-sample_type(id)
status<-statusN<-vector()
status[which(type=="11")]<-"N"
status[which(type=="01")]<-"T"
status<-as.factor(status)
statusN[which(type=="11")]<-0
statusN[which(type=="01")]<-1
statusN<-as.factor(statusN)

sample_name<-function(x){
  substr(as.character(x),1,12)
}
name<-sample_name(id)


#extract mean methylation in Tumour and Normal TCGA samples for Table S4

#find overlap between DMRs selected for validation and probes in TCGA 450K 
OLs<-findOverlaps(OrigDMR_GR,fTab_ranges)
length(unique(OLs@from))#17/18 regions covered by TCGA 450K 
length(OLs@from)#comprising 96 probes in total
idx_OLs<-OLs@to

#extract methylation beta values at the 96 probes overlapping selected DMRs
ALLmeth_cpg_DMP<-beta2[idx_OLs,]

#extract annotation data for the 96 probes overlapping selected DMRs
ALLannot_meth_cpg_DMP<-data.frame(OrigDMR_names[OLs@from],OrigDMR_GR[OLs@from],fTab_ranges[OLs@to])

#calculate mean methylation across each selected DMR in the TCGA data for each sample
uniqOLs<-unique(OLs@from)

methlist_ALL<-list()
for(i in 1:length(uniqOLs)){
	temp<-which(OLs@from==uniqOLs[i])
	if(length(temp)>1){
		temp2<-colMeans(ALLmeth_cpg_DMP[temp,])
	}else if(length(temp)==1){
		temp2<-ALLmeth_cpg_DMP[temp,]
	}
methlist_ALL[[i]]<-temp2	
}

meth_cpg_DMR_ALL<-do.call(rbind,methlist_ALL)
rownames(meth_cpg_DMR_ALL)<-OrigDMR_names[uniqOLs]

#write out mean methylation of tumour and normal samples at each DMR in TCGA as csv file, with t-test statistics for Table S4
TvsNtableDMR<-matrix(data=NA,ncol=4,nrow=nrow(meth_cpg_DMR_ALL))
for(i in 1:nrow(meth_cpg_DMR_ALL)){
	temp<-t.test(meth_cpg_DMR_ALL[i,]~statusN)
	TvsNtableDMR[i,]<-c(temp$estimate,temp$statistic,temp$p.)
}
colnames(TvsNtableDMR)<-c("MethNorm","MethTum","T_est","T_pval")
rownames(TvsNtableDMR)<-OrigDMR_names[unique(OLs@from)]	

#export data for Table S4
write.csv(TvsNtableDMR,"WGBS/results/TCGA_TvsN_AverageDMR.csv")



#re-order DMRs in meth_cpg_DMR_ALL by genomic position

genes<-unique(ALLannot_meth_cpg_DMP$OrigDMR_names.OLs.from.)
dmr_coord<-ALLannot_meth_cpg_DMP[match(genes,ALLannot_meth_cpg_DMP$OrigDMR_names.OLs.from.),2:4]

numonly<-sapply(dmr_coord$seqnames,function(x){substr(as.character(x),4,nchar(as.character(x)))})
reord_idx<-order(as.numeric(numonly),dmr_coord$start)

Av_DMRs_ord<-meth_cpg_DMR_ALL[reord_idx,]

#for plotting add an NA where CACNA2D4 would have been if there were probes on the array

Av_DMRs_ord2<-rbind(Av_DMRs_ord[1:9,],rep(NA,ncol(Av_DMRs_ord)),Av_DMRs_ord[10:17,])

rownames(Av_DMRs_ord2)<-paste("DMR",1:18,sep="")

#prepare data for plotting
norm<-1:45
tumour<-46:437

# Load in blood data, generated in '04_02_ExtractingBlood_DMRmeth.R'
Blood<-as.matrix(read.csv("WGBS/results/AvMeth_DMRs_maleblood.csv", row.names=1))

#to match up with WGBS boxplot:
seq_tumour<-seq(from=1,by=7,length.out=18)
seq_norm<-seq(from=2,by=7,length.out=18)
seq_blood<-seq(from=3,by=7,length.out=18)

#export boxplot for Figure S4b
pdf("WGBS/plots/MeanMeth_bygroup_all_withBlood.pdf",
width=15,height=5,useDingbats=FALSE)
plot(seq_tumour,Av_DMRs_ord2[,1]*100,ylim=c(0,100),xlim=c(0,max(seq_blood)),xlab="",ylab="DNA methylation (%)",xaxt="n",pch=16, col="white",cex=0.8)

for(i in 1:nrow(Av_DMRs_ord2)){
  boxplot(Av_DMRs_ord2[i,tumour]*100,boxwex=1.5,at=seq_tumour[i],ylim=c(0,100),xaxt="n",xlab="",ylab="DNA Methylation (%)",add=T,outline=F,yaxt="n",col="red")
}

for(i in 1:length(seq_norm)){
  boxplot(Av_DMRs_ord2[i,norm]*100,boxwex=1.5,at=seq_norm[i],ylim=c(0,100),xaxt="n",xlab="",ylab="DNA Methylation (%)",add=T,outline=F,yaxt="n",col="blue")
}

for(i in 1:18){
  boxplot(Blood[i,]*100, boxwex=1.5, at=seq_blood[i],ylim=c(0,100),xaxt="n",xlab="",ylab="DNA Methylation (%)",add=T,outline=F,yaxt="n",yaxt="n",col="orange")
}

a<-rep("DMR",18)
b<-1:18
lab<-paste(a,b,sep="")

axis(side=1, at=seq_tumour+1,labels=lab,las=2)
dev.off()


###export bed files
#take mean of each probe for normal and tumour separately

Nmean<-rowMeans(beta2[,statusN==0],na.rm=T)
Tmean<-rowMeans(beta2[,statusN==1],na.rm=T)

#export as 2 bedGraph files
dfN<-data.frame(seqnames=paste("chr",fTab3$CHR,sep=""),start=fTab3$MAPINFO-1,end=fTab3$MAPINFO,scores=Nmean)

dfT<-data.frame(seqnames=paste("chr",fTab3$CHR,sep=""),start=fTab3$MAPINFO-1,end=fTab3$MAPINFO,scores=Tmean)

options(scipen=10)

#write out bedgraph files for IGV plots in Figure S1
write.table(dfN, file="WGBS/results/Normal_TCGA.bedGraph", quote=F, sep="\t", row.names=F, col.names=F)
write.table(dfT, file="WGBS/results/Tumour_TCGA.bedGraph", quote=F, sep="\t", row.names=F, col.names=F)




###Comparing methylation and expression data from TCGA

#RNA-seq level 3 data for TCGA v2.1.14.0 downloaded from https://tcga-data.nci.nih.gov and made into an RData file
#see command.log, run_samples.R and run_expr_values.R

#load in processed RNA-seq data for TCGA
load("WGBS/data/unc.edu_PRAD.IlluminaHiSeq_RNASeqV2.1.14.0.genes.Rdata")

#create vector of sample IDs of 450K data analysed in this study
id<-targetsboth3$Sample_Name[match(colnames(beta2),as.character(targetsboth3$Array.Data.File))]

#function for identifying sample
sample_identifier<-function(x){
  substr(as.character(x),1,16)
  
}

targs<-sample_identifier(id)

#Extract expression data matching 450K data 
RNA<-X[,match(targs,colnames(X))]
#note, same dimensions at beta2. 437 columns, but some of these have no data

#remove "|" from RNA rownames
newRNA<-vector()
for(i in 1:nrow(RNA)){
  newRNA[i]<-strsplit(rownames(RNA)[i],"[|]")[[1]][1]
}

DMR_names<-c("CD34","AC074091.13","EPHB3","PRDM8","MARCH6","CDO1","ZNF655","CDH17","CCDC67","EFCAB4B","PAH","LECT1","CALM1","PARP6","LRRC37A3","TBX1","BHLHB9")

#note CCDC67 also known as DEUP1
#note EFCAB4B also known as CRACR2A
#note LECT1 also known as CNMD

#Identify rows in expression data matching selected DMRs
vec<-list()
for(i in 1:length(DMR_names)){
  vec[[i]]<-match(DMR_names[i],newRNA)
}
vec2<-unlist(vec)

RNAmatchmeth<-RNA[vec2,]

#export plots of correlation between methylation and expression at nearest protein-coding gene for Figure S5
pdf("WGBS/plots/PlotAvDMR_Expr.pdf",width=12)
par(mfrow=c(2,4))
for(i in 1:nrow(RNAmatchmeth)){
 checkR<-length(which(is.na(RNAmatchmeth[i,])==TRUE))
    checkM<-length(which(is.na(Av_DMRs_ord[i,])==TRUE))
    if(checkR<437 & checkM <437){
        rna<-log2(RNAmatchmeth[i,])
        rna2<-which(rna=="-Inf")
        rna[rna2]<-NA
        temp<-cor.test(rna,Av_DMRs_ord[i,],na.rm=T)#pearsons
    plot(Av_DMRs_ord[i,]*100,rna,pch=16,col="blue",ylab="Gene experession (log2)",xlab="DNA methylation (%)",xlim=c(0,100),
    main=paste(rownames(RNAmatchmeth)[i]," - cor=",round(temp$est,2)," p=",signif(temp$p.,2),sep=""))
    points(Av_DMRs_ord[i,statusN==1]*100,rna[statusN==1],pch=16,col="red" )
    abline(lm(as.numeric(rna)~as.numeric(Av_DMRs_ord[i,]*100)))
    }
}
dev.off()

#write out correlation statistics for Table S4
MethvsExprtab<-matrix(data=NA,ncol=2,nrow=nrow(RNAmatchmeth))
for(i in 1:nrow(RNAmatchmeth)){
    checkR<-length(which(is.na(RNAmatchmeth[i,])==TRUE))
       checkM<-length(which(is.na(Av_DMRs_ord[i,])==TRUE))
       if(checkR<437 & checkM <437){
           rna<-log2(RNAmatchmeth[i,])
           rna2<-which(rna=="-Inf")
           rna[rna2]<-NA
           temp<-cor.test(rna,Av_DMRs_ord[i,],na.rm=T)#pearsons
       }
    MethvsExprtab[i,]<-c(temp$est,temp$p.)
}
rownames(MethvsExprtab)<-DMR_names
colnames(MethvsExprtab)<-c("Pearson_cor","P")
write.csv(MethvsExprtab,"WGBS/results/CorAvDMR_Expr.csv")
#Boxplot for figure 10 - gene expression of CACNA2D4 in Tumour vs normal

#adding CACNA2D4 back in
DMR_names<-c("CD34","AC074091.13","EPHB3","PRDM8","MARCH6","CDO1","ZNF655","CDH17","CCDC67","CACNA2D4","EFCAB4B","PAH","LECT1","CALM1","PARP6","LRRC37A3","TBX1","BHLHB9")

#note CCDC67 also known as DEUP1
#note EFCAB4B also known as CRACR2A
#note LECT1 also known as CNMD

#Identify rows in expression data matching selected DMRs
vec<-list()
for(i in 1:length(DMR_names)){
  vec[[i]]<-match(DMR_names[i],newRNA)
}
vec2<-unlist(vec)

RNAmatchmeth<-RNA[vec2,]

#export plot of tumour versus normal expression of CACNA2D4 for Figure S11b
pdf(file="WGBS/plots/TCGA_TvsN.pdf",width=12)
par(mfrow=c(2,4))
for(i in 1:nrow(RNAmatchmeth)){
	checkR<-length(which(is.na(RNAmatchmeth[i,])==TRUE))
	if(checkR<437){
      plot(statusN,log2(RNAmatchmeth[i,]),main=rownames(RNAmatchmeth)[i],ylab="Gene expression (log2)",xaxt="n",outline=F)
      stripchart(log2(RNAmatchmeth[i,])~statusN,add=T,vertical=T,method="jitter",pch=1,col=c("dark blue","red"))
  axis(1,at=c(1,2),labels=c("Normal","Tumour"))
    }}
dev.off()

#Statistics for Figure S11b - gene expression of CACNA2D4 in Tumour vs normal
checkRvec<-vector()
for(i in 1:nrow(RNAmatchmeth)){
	checkR<-length(which(is.na(RNAmatchmeth[i,])==TRUE))
	checkRvec[i]<-checkR
}
namesR<-rownames(RNAmatchmeth)[which(checkRvec==22)]

TvsNtableRNA<-matrix(data=NA,ncol=4,nrow=nrow(RNAmatchmeth))
for(i in 1:nrow(RNAmatchmeth)){
	checkR<-length(which(is.na(RNAmatchmeth[i,])==TRUE))
	checkRvec[i]<-checkR
	if(checkR<437){	
	temp<-log2(RNAmatchmeth[i,])
	temp2<-which(temp=="-Inf")
	temp[temp2]<-NA
	temp3<-t.test(temp~statusN)
	TvsNtableRNA[i,]<-c(temp3$estimate,temp3$statistic,temp3$p.)
}}
colnames(TvsNtableRNA)<-c("ExpNorm","ExpTum","T_est","T_pval")
rownames(TvsNtableRNA)<-rownames(RNAmatchmeth)

write.csv(TvsNtableRNA,"WGBS/results/TCGA_TvsN.csv")
