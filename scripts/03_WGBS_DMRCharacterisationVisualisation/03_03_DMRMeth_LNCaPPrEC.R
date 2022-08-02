#Script to extract methylation at all DMRs from LNCaP and PrEC WGBS data
#Generates data used in Table S4
#Generates 2 columns of heatmap for Figure 1c

#load libraries
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#Read in BigWig files corresponding to data in GEO GSE86832
#In house data, large file so available at users request
#Note, LNCaP and PrEC bigwig files loaded into IGV as tracks to create Figure S1
bw.file <- "metadata/LNCaP.bw"
bw.fileP <- "metadata/PrEC.bw"

# Read in table of DMRs generated from '02_01_DMRcalling_plot.R' and convert to a Genomic Ranges object
DMRs<-read.csv("WGBS/results/Analysis_DMRs.csv")

selection<-GRanges(seqnames=DMRs$seqnames,IRanges(start=DMRs$start,end=DMRs$end))

# Extract BigWig methylation data from LNCaP and PrEC at each DMR
list_L<-list()
for(i in 1:nrow(DMRs)){
  list_L[[i]]<- import.bw(con = BigWigFile(bw.file),
                          selection = BigWigSelection(selection[i,]))
}

list_P<-list()
for(i in 1:nrow(DMRs)){
  list_P[[i]]<- import.bw(con = BigWigFile(bw.fileP),
                          selection = BigWigSelection(selection[i,]))
}

#Calculate mean methylation at each DMR in LNCaP and PrEC
meanL<-vector()
for(i in 1:length(list_L)){
  meanL[i]<-mean(list_L[[i]]$score)
}

meanP<-vector()
for(i in 1:length(list_P)){
  meanP[i]<-mean(list_P[[i]]$score)
}

MeanDMR<-cbind(DMRs,meanL*100,meanP*100)

#Write out table of mean methylation at each DMR in LNCaP and PrEC for use in Table S4
write.csv(MeanDMR,file="/WGBS/results/LNCaP_PreC_DMRmeth_1420.csv")

#read in order of DMRs used in patient heatmap - index generated in 02_02_DMRmean_heatmap.R
ind<-read.csv("WGBS/plots/DMRheatmap_rowIndex.csv")

#load libraries for plotting heatmap
library(gplots)
library(RColorBrewer)

#re-order LNCaP-PrEC methylation data in same order as DMRs in patient heatmap
LP<-cbind(MeanDMR[,13],MeanDMR[,14])
LP2<-LP[ind$x,]
colnames(LP2)<-c("LNCaP","PrEC")

#create heatmap for LNCaP and PrEC methylation at DMRs to append to patient heatmap in Figure 1c

pdf("WGBS/plots/DMRheatmap_LNCaPPrEC.pdf")
heatmap.2(LP2[1420:1,], trace = 'n', dendrogram = "none", labRow=F, cexCol = 1.5, srtCol=0, col=rev(brewer.pal(11,"RdBu")), density.info="none", key.title="", na.color="grey", Rowv=F, Colv=F)
dev.off()


