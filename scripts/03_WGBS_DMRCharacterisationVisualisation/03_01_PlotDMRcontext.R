#Script for assessing overlap of DMRs with CpG islands and genic regions
#Generates Figure 1e

#load libraries
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#define function required for annotating DMRs (from aaRon) (https://github.com/astatham/aaRon/)
strip <- function(temp) {
    if (class(temp)=="GRangesList") return(endoapply(temp, strip))
    temp <- temp[seqnames(temp) %in% names(Hsapiens)[1:24]]
    seqlevels(temp) <- names(Hsapiens)[1:25]
    values(temp) <- NULL
    strand(temp) <- "*"
    reduce(temp)
}

#Read in table of DMR annotation generated from '02_04_AnnotateDMRs.R'
DMRs<-read.csv("WGBS/results/Annotated_Analysis_DMRs.csv")

#split DMRs into whether hyper or hypomethylated
hyperDMRs<-DMRs[which(DMRs$meanbetafc>0),]
hypoDMRs<-DMRs[which(DMRs$meanbetafc<0),]

#Quantify the overlap between DMRs and CpG island regions
#CpG islands
proportioncpg<-function(x){colSums(x[,c("CpGisland","CpGshores","nonCpG")])/nrow(x)}

a<-proportioncpg(hyperDMRs)
b<-proportioncpg(hypoDMRs)

c<-rbind(a,b)

#Quantify the overlap between DMRs and Protein-coding genic regions
proportioncpg<-function(x){colSums(x[,c("promoter_prot","genebody_prot","intergenic_prot")])/nrow(x)}

aaa<-proportioncpg(hyperDMRs)
bbb<-proportioncpg(hypoDMRs)

ccc<-rbind(aaa,bbb)


#read in background genomic regions generated from '01_01_DownloadCpGislandData.R' and '02_03_DownloadTranscriptData.R'
load("WGBS/annotationdata/CpGislands.Rdata")
load("WGBS/annotationdata/hg19.tx.Rdata")

#comparing CpG island data to whole genome data to establish regions outside of CpG island (i.e. define non-CpG island regions)
temp<-resize(CpGislands, width(CpGislands)+4000, fix="center")
temp2<-reduce(temp)
temp3<-setdiff(temp2,CpGislands)

genomeGR <- GRanges(seqlevels(tx), IRanges(1, seqlengths(tx)))
diff1<-setdiff(genomeGR, CpGislands)
diff2<-setdiff(diff1,CpGshores)

#calculate percent of genome covered by CpG islands, shores and non-CpG islands
summaryCpG<-data.frame(c(length(CpGislands),length(CpGshores),length(diff2)),
                       c(sum(width(CpGislands)),sum(width(CpGshores)),sum(as.numeric(width(diff2)))),
                       c(mean(width(CpGislands)),mean(width(CpGshores)),mean(width(diff2))))
rownames(summaryCpG)<-c("CpGislands","CpGshores","Non_CpG")
colnames(summaryCpG)<-c("NumberRegions","LengthRegions (bp)","MeanWidth")

PercentageGenome<-round((summaryCpG[,2]/sum(summaryCpG[,2]))*100,2)
summaryCpG2<-cbind(summaryCpG,PercentageGenome)

w_isl<-(sum(width(CpGislands))/sum(as.numeric(width(genomeGR))))*100
w_shore<-(sum(width(CpGshores))/sum(as.numeric(width(genomeGR))))*100
w_non<-(sum(as.numeric(width(diff2)))/sum(as.numeric(width(genomeGR))))*100

df2<-data.frame(c(w_non,w_shore,w_isl),factor(c("Non_CpGisland","CpGshores","CpGislands"),levels=c("Non_CpGisland","CpGshores","CpGislands")))
colnames(df2)<-c("Meth","Region")

#calculate percent of genome covered by Protein coding gene promoters, gene bodies and non-genic regions

txp<-tx[tx$tx_type=="protein_coding"]
#promoters
promotersGR_prot <- suppressWarnings(reduce(strip(resize(resize(txp, 1, fix="start"), 4000, fix="center"))))
genebodyGR_prot <- suppressWarnings(setdiff(reduce(strip(txp)), promotersGR_prot))
intergenicGR_prot <- suppressWarnings(setdiff(setdiff(genomeGR, genebodyGR_prot), promotersGR_prot))

w_prom_prot<-(sum(width(promotersGR_prot))/sum(as.numeric(width(genomeGR))))*100
w_gb_prot<-(sum(width(genebodyGR_prot))/sum(as.numeric(width(genomeGR))))*100
w_int_prot<-(sum(width(intergenicGR_prot))/sum(as.numeric(width(genomeGR))))*100

summaryCpG4<-data.frame(c(length(promotersGR_prot),length(genebodyGR_prot),length(intergenicGR_prot)),
                        c(sum(width(promotersGR_prot)),sum(width(genebodyGR_prot)),sum(as.numeric(width(intergenicGR_prot)))),
                        c(mean(width(promotersGR_prot)),mean(width(genebodyGR_prot)),mean(width(intergenicGR_prot))),
                        c(w_prom_prot,w_gb_prot,w_int_prot))
rownames(summaryCpG4)<-c("Promoters","Genebody","Intergenic")
colnames(summaryCpG4)<-c("NumberRegions","LengthRegions (bp)","MeanWidth","PercentageGenome")

df_prot<-data.frame(c(w_int_prot,w_gb_prot,w_prom_prot),factor(c("Intergenic","Gene Body","Promoter"),levels=c("Intergenic","Gene Body","Promoter")))
colnames(df_prot)<-c("Meth","Region")

#make bar charts of overlap between DMRs and genomic features for Figure 1e
pdf("WGBS/plots/DMROverlap.pdf",width=10)
par(mfrow=c(1,3))
c2<-rbind(c*100,rev(df2[,1]))
barplot(t(c2),names.arg=c("Hyper","Hypo","Genomic bkgrd"),main="Overlap with CpG islands",legend.text=colnames(c2),args.legend="top left", ylab="Percent of DMRs")

cc3<-rbind(ccc*100,rev(df_prot[,1]))
barplot(t(cc3),names.arg=c("Hyper","Hypo","Genomic bkgrd"),main="Overlap with Genic Regions (prot only)",legend.text=colnames(cc3),args.legend="top left", ylab="Percent of DMRs")
dev.off()

#Calculate median distance to TSS for manuscript
medTSS_hyper<-median(hyperDMRs$distanceTSS)
# 745
medTSS_hypo<-median(hypoDMRs$distanceTSS)
#12799
wilcox.test(hyperDMRs$distanceTSS,hypoDMRs$distanceTSS)$p.
#4.643222e-21



