#Script contains details of WGBS data processing
#Generates WGBS PCA plot for Figure 1a
#Generates plots of methylation according to CpG context and at repetitive elements for Figure 1b and Figure S2

#note, this script was run on an HPC rather than a personal computer as the BigTable data is large

#note, script includes functions originally written by Aaron Statham as part of aaRon package: https://github.com/astatham/aaRon

#load libraries and prepare environment
library(data.table) 
library(GenomicRanges)
library(RColorBrewer)
library(ggplot2)
library(ggbiplot) 
library(vegan) 
library (BiodiversityR)
library(REMP)
library(dunn.test)

options(scipen=100)

#set user defined working directory and create new folders to output results
BASEDIR<-*user_defined_file_path*
setwd(paste(BASEDIR,"/WGBS",sep=""))
dir.create("data")
dir.create("plots")
setwd(BASEDIR)

#to run this script user should download WGBS data: GSE158927_BigTable.tsv.gz from GEO GSE158927 into newly created "~/WGBS/data/"

#read in WGBS methylated C and coverage data (bigTable)
tab <- fread("WGBS/data/GSE158927_BigTable.tsv.gz")

#note, to create the BigTable bisulfite reads in fastq files were aligned to the human genome (hg19) in 2017 using an internally developed pipeline specialised for Garvan hardware. This pipeline has since been re-written for general use on high performance computing clusters that use e.g. SGE/PBS and is publicly available at: https://github.com/luuloi/Meth10X

#convert WGBS data table into a Genomic Ranges object
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
                      
                      
#define tumour samples
tum<-as.character(samples$Sample[which(samples$Tissue== "Cancer_Dead" | samples$Tissue== "Cancer_Alive")])

#define lethal and non-lethal samples
lethal<-as.character(samples$Sample[which(samples$Tissue== "Cancer_Dead" )])
nonlethal<-as.character(samples$Sample[which(samples$Tissue== "Cancer_Alive")])

#define normal samples
normals<-as.character(samples$Sample[which(samples$Tissue== "Normal")])

vectorstatus<-c(tum,normals)

#function for calculating methylation percentage from bigTable 
methRatios <- function(x, samples, minCov=5) {
  stopifnot(all(c("Sample", "C", "cov") %in% colnames(samples)))
  tmp <- as.matrix(values(x)[samples$cov])
  tmp[tmp<minCov] <- NA
  values(x) <- as.matrix(values(x)[samples$C])/tmp
  names(values(x)) <- samples$Sample
  x
}

#calculate methylation percentage from bigTable using methRatios function
CpGs.ratios <- methRatios(CpGs,samples,minCov=5)
ratios <- as.data.frame(values(CpGs.ratios))

#from 'ratios' remove rows conatining NAs from low coverage CpG sites
idxs <- apply(ratios, 1, function (x) sum(is.na(x))==0)
ratios2 <- ratios[idxs,]

###Create PCA plot for Figure 1a###

#define colour palette
pal<-brewer.pal(4,"Set1")

#PCA with vegan package
prin_comp<-rda(t(ratios2), scale=TRUE)
pca_scores<-scores(prin_comp)

#extract variances explained by PCA axis 1 and 2 using function from BiodiversityR package
sig <- PCAsignificance (prin_comp,axes=5)
k1 <- round(sig[2,1], 1)
k2 <- round(sig[2,2], 1)

#save out PCA plot
pdf("WGBS/plots/PCAplot.pdf", useDingbats=FALSE)

plot(pca_scores$sites[,1],
     pca_scores$sites[,2],
    xlab="", ylab="",
     pch=16,
    col=pal[samples$Tissue],
    xlim=c(-60,40),
    ylim=c(-50,50))
text(pca_scores$sites[,1],pca_scores$sites[,2]+3, labels = colnames(ratios2), col=pal[samples$Tissue])

ordiellipse(prin_comp,samples$Tissue,conf=0.95)
title(main=paste("p=", nrow(ratios2), ",", k1+k2, "% of variance explained"), xlab=paste(k1, "% of variance explained"), ylab=paste(k2, "% of variance explained"))

dev.off()


###Create functions for Figure 1b and Figure S2###

#functions for calculating methylation percentages from count data from each supplied window from https://github.com/astatham/aaRon
methWindowRatios <- function(x, windows, samples, minCov=5, mc.cores=1) {
  stopifnot(all(c("Sample", "C", "cov") %in% colnames(samples)))
  tmp <- simplify2array(mclapply(1:nrow(samples),
                                 function(i) overlapRatios(x, windows, samples$C[i], samples$cov[i]),
                                 mc.cores = mc.cores))
  if (length(windows)==1) tmp <- matrix(tmp, nrow=1)
  values(windows) <- tmp
  names(values(windows)) <- samples$Sample
  windows
}

overlapRatios <- function(x, y, C, cov, minCov=5, na.rm=FALSE) {
  C.sum <- overlapSums(x, y, values(x)[[C]], na.rm=na.rm)
  cov.sum <- overlapSums(x, y, values(x)[[cov]], na.rm=na.rm)
  rat <- C.sum/cov.sum
  rat[cov.sum<minCov] <- NA
  rat
}

overlapSums <- function(x, y, z, na.rm=FALSE) {
  stopifnot(class(x)=="GRanges")
  stopifnot(class(y)=="GRanges")
  stopifnot(class(z)=="numeric" | class(z)=="integer")
  stopifnot(length(x)==length(z))
  ov <- as.matrix(findOverlaps(y, x))
  ovSums <- rep(as.numeric(NA), length(y))
  ovSums[unique(ov[,1])] <- viewSums(Views(z[ov[,2]], ranges(Rle(ov[,1]))), na.rm=na.rm)
  ovSums
}

#function for methylation density plots (from https://github.com/astatham/aaRon)
methDensityPlot <- function(x, samples, title, cols=NULL) {
  y <- Group <- Sample <- NULL #to shut up R CMD check
  stopifnot(all(c("Sample", "C", "cov") %in% colnames(samples)))
  if (is.null(samples$Group)) {
    samples$Group <- samples$Sample
  }
  if (!is.null(cols)) {
    stopifnot(all(names(cols) %in% samples$Group))
    stopifnot(all(samples$Group %in% names(cols)))
  }
  if (class(x)=="GRanges") x <- as.matrix(values(x))
  stopifnot(ncol(x)==nrow(samples))
  dists <- do.call(rbind, lapply(1:ncol(x), function(i) {
    x <- density(x[,i], na.rm=TRUE)
    data.frame(Sample=samples$Sample[i], Group=samples$Group[i], x=x$x, y=x$y, stringsAsFactors=FALSE)
  }))
  
  p <- ggplot(dists, aes(x=x, y=y, color=Group, group=Sample)) +
    ggtitle(title) + 
    xlab("Methylation ratio") + ylab("Relative Frequency") +
    geom_line() + xlim(-0.05, 1.05)
  if (!is.null(cols)) p <- p + scale_color_manual(values=cols)
  p
}


###Create CpG Island and CpG Shore plots for Figure 1b and Figure S2a###

#load in CpG island genomic positions generated from '01_01_DownloadCpGislandData.R'
load("WGBS/annotationdata/CpGislands.Rdata")

#define sample colours
groups <- c("Lethal","Nonlethal","Normal")
cols <- brewer.pal(length(groups), "Dark2")[1:length(groups)]
names(cols) <- groups
makeGroup<-c(rep("Lethal",length(lethal)),rep("Nonlethal",length(nonlethal)),rep("Normal",length(normals)))
samples$cols <- cols[makeGroup]
samples$Group <- makeGroup

#create factor of sample groups
fac<-factor(samples$Group,levels=c("Normal","Nonlethal","Lethal"))

#Figure 1b
#extract methylation data at each CpG island for each sample
CpGislands <- as(CpGislands, "GRanges")
islands<-methWindowRatios(CpGs, CpGislands, samples)

#calculate median methylation at CpG islands for each sample
cm_island2<-colMedians(as.matrix(values(islands)),na.rm=T)

#save out boxplots of median methylation at CpG islands for each sample
pdf("WGBS/plots/Boxplots_islands_MedianByGroup.pdf",useDingbats=FALSE)
boxplot(cm_island2*100~fac,ylab="DNA methylation (%)",xlab="",col=cols[c(3,2,1)],ylim=c(0,20))
stripchart(cm_island2*100~fac,add=T,vert=T,pch=16,method="jitter")
dev.off()

#test for difference in median CpG island methylation between groups
dunn.test(cm_island2*100,fac)

#extract methylation data at each CpG shore for each sample
CpGshores <- as(CpGshores, "GRanges")
shores<-methWindowRatios(CpGs, CpGshores, samples)

#calculate median methylation at CpG shores for each sample
cm_shores2<-colMedians(as.matrix(values(shores)),na.rm=T)

#save out boxplots of median methylation at CpG shores for each sample
pdf("WGBS/plots/Boxplots_shores_MedianByGroup.pdf",useDingbats=FALSE)
boxplot(cm_shores2*100~fac,ylab="DNA methylation (%)",xlab="",col=cols[c(3,2,1)],ylim=c(50,70))
stripchart(cm_shores2*100~fac,add=T,vert=T,pch=16,method="jitter")
dev.off()

#test for difference in median CpG island methylation between groups
dunn.test(cm_shores2*100,fac)

#Figure S2a

#make vector of which samples are in which group
lethal_idx<-which(samples$Tissue=="Cancer_Dead")
nonlethal_idx<-which(samples$Tissue=="Cancer_Alive")
normal_idx<-which(samples$Tissue=="Normal")

##save out CpG island density plots - averaged within group
pdf("WGBS/plots/CpG_Islands_averagedgroups.pdf")
a<-rowMeans(as.matrix(values(islands[,lethal_idx])),na.rm=T)
b<-rowMeans(as.matrix(values(islands[,nonlethal_idx])),na.rm=T)
c<-rowMeans(as.matrix(values(islands[,normal_idx])),na.rm=T)
d<-cbind(a,b,c)
samples3<-samples[c(lethal_idx[1],nonlethal_idx[1],normal_idx[1]),]
samples3$Sample<-samples3$Group
print(methDensityPlot(d,
samples3, paste0("Methylation distribution - CpG Islands"), cols))
dev.off()

#save out CpG shore density plots - averaged within group
pdf("WGBS/plots/CpG_Shores_averagedgroups.pdf")
a<-rowMeans(as.matrix(values(shores[,lethal_idx])),na.rm=T)
b<-rowMeans(as.matrix(values(shores[,nonlethal_idx])),na.rm=T)
c<-rowMeans(as.matrix(values(shores[,normal_idx])),na.rm=T)
d<-cbind(a,b,c)
print(methDensityPlot(d,
samples3, paste0("Methylation distribution - CpG Shores"), cols))
dev.off()


###Create plots of methylation at repetitive elements for Figure 1b and Figure S2a&b###

#extract genomic regions corresponding to each type of repetitive element using REMP package

Alus<-initREMP(arrayType = c("Sequencing"),REtype = c("Alu"),annotation.source = c("UCSC"),genome = c("hg19"),Seq.GR=CpGs)

Line<-initREMP(arrayType = c("Sequencing"),REtype = c("L1"),annotation.source = c("UCSC"),genome = c("hg19"),Seq.GR=CpGs)

LTR<-initREMP(arrayType = c("Sequencing"),REtype = c("LTR"),annotation.source = c("UCSC"),genome = c("hg19"),Seq.GR=CpGs)

#extract methylation data at each type of repetitive element for each sample
Alu_meth_2<-methWindowRatios(CpGs, Alus@RECpG, samples)

Line_meth_2<-methWindowRatios(CpGs, Line@RECpG, samples)

LTR_meth_2<-methWindowRatios(CpGs, LTR@RECpG, samples)

#calculate median methylation at each type of repetitive element for each sample
cm_Alu_2<-colMedians(as.matrix(values(Alu_meth_2)),na.rm=T)

cm_Line_2<-colMedians(as.matrix(values(Line_meth_2)),na.rm=T)

cm_LTR_2<-colMedians(as.matrix(values(LTR_meth_2)),na.rm=T)


#save out boxplots of median methylation at each type of repetitive element for each sample - Figure 1b and Figure S2b
pdf("WGBS/plots/Boxplots_Alu_MedianByGroup.pdf")
boxplot(cm_Alu_2*100~fac,ylab="DNA methylation (%)",xlab="",col=cols[c(3,2,1)],ylim=c(90,100))
stripchart(cm_Alu_2*100~fac,add=T,vert=T,pch=16,method="jitter")
dev.off()

#test for difference in median Alu methylation between groups
dunn.test(cm_Alu_2*100,fac)

pdf("WGBS/plots/Boxplots_Line_MedianByGroup.pdf")
boxplot(cm_Line_2*100~fac,ylab="DNA methylation (%)",xlab="",col=cols[c(3,2,1)],ylim=c(70,100))
stripchart(cm_Line_2*100~fac,add=T,vert=T,pch=16,method="jitter")
dev.off()

#test for difference in median LINE methylation between groups
dunn.test(cm_Line_2*100,fac)

pdf("WGBS/plots/Boxplots_LTR_MedianByGroup.pdf")
boxplot(cm_LTR_2*100~fac,ylab="DNA methylation (%)",xlab="",col=cols[c(3,2,1)],ylim=c(70,100))
stripchart(cm_LTR_2*100~fac,add=T,vert=T,pch=16,method="jitter")
dev.off()

#test for difference in median LINE methylation between groups
dunn.test(cm_LTR_2*100,fac)

#save out density plots of methylation at each type of repetive element -  averaged within group - Figure S2a
#Alu density plots - averaged within group
pdf("WGBS/plots/Line_averagedgroups.pdf")
a<-rowMeans(as.matrix(values(Alu_meth_2[,lethal_idx])),na.rm=T)
b<-rowMeans(as.matrix(values(Alu_meth_2[,nonlethal_idx])),na.rm=T)
c<-rowMeans(as.matrix(values(Alu_meth_2[,normal_idx])),na.rm=T)

d<-cbind(a,b,c)
print(methDensityPlot(d,
samples3, paste0("Methylation distribution - Alu"), cols))
dev.off()

#LTR density plots - averaged within group
pdf("WGBS/plots/LTR_averagedgroups.pdf")
a<-rowMeans(as.matrix(values(LTR_meth_2[,lethal_idx])),na.rm=T)
b<-rowMeans(as.matrix(values(LTR_meth_2[,nonlethal_idx])),na.rm=T)
c<-rowMeans(as.matrix(values(LTR_meth_2[,normal_idx])),na.rm=T)

d<-cbind(a,b,c)
print(methDensityPlot(d,
samples3, paste0("Methylation distribution - LTR"), cols))
dev.off()

#Line density plots - averaged within group
pdf("WGBS/plots/Line_averagedgroups.pdf")
a<-rowMeans(as.matrix(values(Line_meth_2[,lethal_idx])),na.rm=T)
b<-rowMeans(as.matrix(values(Line_meth_2[,nonlethal_idx])),na.rm=T)
c<-rowMeans(as.matrix(values(Line_meth_2[,normal_idx])),na.rm=T)

d<-cbind(a,b,c)
print(methDensityPlot(d,
samples3, paste0("Methylation distribution - Line1"), cols))
dev.off()
