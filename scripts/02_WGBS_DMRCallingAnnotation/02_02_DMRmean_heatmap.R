#Script for extracting data from bigTable to calculate mean methylation for lethal, non-lethal and normal samples for each of the 1420 DMRs
#Generates data for Table 2 and Tables S2 & S4
#Generates heatmap for Figure 1c

#note, this script was run on an HPC rather than a personal computer as the BigTable data is large

#load libraries
library(data.table)
library(bsseq)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(DSS)
library(gplots)
library(RColorBrewer)

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#####extracting methylation data at DMRs

#read in WGBS methylated C and coverage data (bigTable)
tab <- fread("WGBS/data/GSE158927_BigTable.tsv.gz")

#turn WGBS data table into a Genomic Ranges object
tab <- tab[-grep("_", tab$chr),]
CpGs <- GRanges(tab$chr, IRanges(tab$position, width=1))
values(CpGs) <- as.data.frame(tab)[,-c(1:3)]
names(values(CpGs)) <- paste("X", names(values(CpGs)), sep='')
rm(tab)

#read in sample table
samples <- read.csv("metadata/WGBS_samplegroups.csv")

#reformat sample table so names match those in CpGs GRanges object, ready for analysis
samples <- data.frame(Sample=samples$Sample,
                      C=paste("X", samples$Sample, ".C", sep=""),
                      cov=paste("X", as.character(samples$Sample), ".cov", sep=""),
                      Tissue=samples$Group,
                      row.names=samples$Sample,
                      stringsAsFactors=FALSE)
                      
#read in DMR results table generated from '02_01_DMRcalling_plot.R'
analysisdmrs <- read.csv("WGBS/results/Analysis_DMRs.csv")
analysisdmrs <- GRanges(analysisdmrs$seqnames, IRanges(analysisdmrs$start, analysisdmrs$end), mcols=analysisdmrs[,7:12])
colnames(mcols(analysisdmrs)) <- gsub("mcols.", "", colnames(mcols(analysisdmrs)))


#Extract per-DMR methylation for all samples for analysis
ols <- findOverlaps(CpGs, analysisdmrs)
an_counts <- GRanges(seqnames=seqnames(analysisdmrs),
                         ranges=ranges(analysisdmrs),
                         mcols=t(sapply(unique(subjectHits(ols)), function (x) colSums(as.matrix(values(CpGs)[queryHits(ols)[subjectHits(ols)==x],])))))

#summing all CpGs (meth and coverage) within a DMR
colnames(values(an_counts)) <- gsub("mcols.", "", colnames(values(an3_counts)))
coverage <- mcols(an_counts[,samples$cov])
meth <- mcols(an_counts[,samples$C])

temp<-match(1:length(analysisdmrs),unique(subjectHits(ols)))

meth_m<-(as.matrix(meth))[temp,]
coverage_m<-(as.matrix(coverage))[temp,]
an_DMR_ratios <- meth_m/coverage_m

an_DMRs <- GRanges(seqnames(analysisdmrs), ranges(analysisdmrs), mcols=an_DMR_ratios)

#write out tables and RData files of methylation at each DMR
#Count of Cs at CpG sites within each DMR
write.csv(meth_m,"WGBS/results/C_analysisDMRs.csv")

#Total coverage at CpG sites within each DMR
write.csv(coverage_m,"WGBS/results/cov_analysisDMRs.csv")

#Proportion methylated Cs at CpG sites within each DMR
#Data used in Tables S2 and S4
write.csv(an_DMR_ratios,"WGBS/results/meth_analysisDMRs.csv")

#Proportion methylated Cs at CpG sites within each DMR in GRanges format
save(an_DMRs,file="WGBS/results/methGR_analysisDMRs.RData")

##Extract CpG-level methylation for top 50 DMRs, extended to 10,000bp

#create GRanges object of top 50 DMRs only
topdmrs<-analysisdmrs[1:50]

#maximum dmr size around 5000bp, so resizing all regions to 10,000bp from the center of the DMR
res<-resize(topdmrs,width=10000,fix="center",ignore.strand=FALSE)

#extracting methylation for CpGs in resized DMRs
an3_res_CpGs<-list()
for(i in 1:50){
    ols <- findOverlaps(CpGs, res[i])
    CpG_dmr1 <- CpGs[queryHits(ols)]
    coverage <- mcols(CpG_dmr1[,samples$cov])
    meth <- mcols(CpG_dmr1[,samples$C])
    an_res_ratios <- as.matrix(meth)/as.matrix(coverage)
    #added line to filter out any individual CpGs with coverage < 5
    an_res_ratios[which(as.matrix(coverage)<5)]<-NA
    an_res_CpGs[[i]] <- GRanges(seqnames(CpG_dmr1), ranges(CpG_dmr1), mcols=an_res_ratios)
}
save(an_res_CpGs,file="WGBS/results/methGR_analysisCpGsInDMRs_10000bp.RData")

####Calculating mean methylation at DMRs
DMR <- read.csv("WGBS/results/meth_analysisDMRs.csv")

#get rid of first column
DMR <- DMR[,-1]

colnames(DMR)<-sapply(colnames(DMR),function(x){strsplit(as.character(x),"\\.")[[1]][1]})


#define samples for ordering (including normals)
tumnorm <- paste("X",as.character(samples$Sample),sep="")

DMR<-DMR[,match(tumnorm,colnames(DMR))]

#make vector of which samples are in which group
lethal_idx<-which(samples$Tissue=="Cancer_Dead")
nonlethal_idx<-which(samples$Tissue=="Cancer_Alive")
normal_idx<-which(samples$Tissue=="Normal")

#define lethal,non-lethal and normal samples
lethal<-tumnorm[lethal_idx]
nonlethal<-tumnorm[nonlethal_idx]
normals<-tumnorm[normal_idx]

#average based on sampletype: lethal, nonlethal, normals

avM <- matrix(data=NA,nrow=nrow(DMR),ncol=3)
for(i in 1:nrow(DMR)){
  
  avM[i,1]<-mean(as.numeric(DMR[i,lethal]),na.rm=T)*100
  
  avM[i,2]<-mean(as.numeric(DMR[i,nonlethal]),na.rm=T)*100
  
  avM[i,3]<-mean(as.numeric(DMR[i,normals]),na.rm=T)*100
  
}
colnames(avM)<-c("Lethal mean methylation","Non-lethal mean methylation","Normals mean methylation")


#### Difference in methylation between the groups ####

# Lethal vs non-lethal - and indicate which were hypermethylated or hypomethylated DMRs

diff<-matrix(data=NA,nrow=nrow(DMR),ncol=2)
for (i in 1:nrow(DMR)){
  
  diff[i,1]<-avM[i,1]-avM[i,2]
  
  if (diff[i,1] > 0) {
    diff[i,2] <- "Hypermethylated"
  } else {
    diff[i,2]<- "Hypomethylated"
  }
  
}
colnames(diff)<-c("Mean difference (%)","Direction")

# now bind the two matrices, save out as .csv file. This contains proportion methylated CpGs in each DMR per patient group, used in Table S2 and S4
M <- cbind(avM,diff)
write.csv(M,"WGBS/results/methbygroup_analysisDMRs.csv")

#combine data of proprotion methylated CpGs in each DMR per sample for plotting in heatmap
DMR_mat<-as.matrix(an_DMRs@elementMetadata)

idx_leth<-vector()
for(i in 1:length(lethal)){
  idx_leth[i]<-grep(lethal[i],colnames(DMR_mat))
}

idx_nonleth<-vector()
for(i in 1:length(nonlethal)){
  idx_nonleth[i]<-grep(nonlethal[i],colnames(DMR_mat))
}

idx_norm<-vector()
for(i in 1:length(normals)){
  idx_norm[i]<-grep(normals[i],colnames(DMR_mat))
}

DMR_ord<-DMR[,c(idx_leth,idx_nonleth,idx_norm)]

a<-rep("red",length(idx_leth))
b<-rep("blue",length(idx_nonleth))
c<-rep("green",length(idx_norm))
vectortemp<-c(a,b,c)

results<-read.csv("WGBS/results/Analysis_DMRs.csv")
col1 <- c("purple","dark green")
hyp<-vector(length=nrow(results))
hyp[which(results$meanbetafc>0)]<-1

#save out Heatmap plot for Figure 1c
pdf("WGBS/plots/DMRheatmap.pdf")
A<-heatmap.2(DMR_ord[1:1420,c(normal_idx,nonlethal_idx,lethal_idx)],col=rev(brewer.pal(11,"RdBu")), trace="none",key.title=NA,ColSideColors=vectortemp[c(normal_idx,nonlethal_idx,lethal_idx)],Colv=FALSE,
          RowSideColors = col1[as.factor(hyp)])
dev.off()

#export row order in heatmap for plotting LNCaP/PrEC side by side:
ind<-A$rowInd

#Write out table specifying order of DMRs in heatmap
write.csv(ind,file="/WGBS/plots/DMRheatmap_rowIndex.csv")



