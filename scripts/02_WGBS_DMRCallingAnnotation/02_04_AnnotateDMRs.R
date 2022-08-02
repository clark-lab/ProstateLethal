#Script for annotating DMRs with CpG and genic context and PrEC/LNCaP ChromHMM state
#Generates data for Table 2 and Tables S2 & S4, and in later scripts to generate Figure 1e

#load libraries
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#read R code for annotating DMRs, adapted from a function in the aaRon package https://github.com/astatham/aaRon/
source("metadata/AnnotateRegions_updown_update.R")

#define functions required for annotating DMRs (from aaRon)
strip <- function(temp) {
    if (class(temp)=="GRangesList") return(endoapply(temp, strip))
    temp <- temp[seqnames(temp) %in% names(Hsapiens)[1:24]]
    seqlevels(temp) <- names(Hsapiens)[1:25]
    values(temp) <- NULL
    strand(temp) <- "*"
    reduce(temp)
}

unvalue <- function(x) {
       values(x) <- NULL
       x
}

coverageRatio <- function(query, subject, ratio=TRUE) {
  if (length(unique(strand(subject)))>1) {
    warning("Supplied 'subject' contains strand information - this is ignored by coverageRatio")
    strand(subject) <- "*"
  }
  if (nrow(as.matrix(findOverlaps(subject, drop.self=TRUE, drop.redundant=TRUE)))>0) {
    warning("Supplied 'subject' contains overlapping ranges, reducing...")
    subject <- reduce(subject)
  }
  # Views uses  RangesLists which do not preserve order
  oo <- GenomicRanges::order(query)
  chr.oo <- unique(as.character(seqnames(query[oo])))
  seqlevels(query) <- chr.oo
  query.cov <- viewSums(Views(coverage(subject)[chr.oo], as(query[oo], "IntegerRangesList")))
  query$result[oo] <- if (ratio) unlist(query.cov)/width(query[oo]) else unlist(query.cov)
  as.numeric(query$result)
}

#load in CpG island data and hg19 transcript annotation data generated from '01_01_DownloadCpGislandData.R' and '02_03_DownloadTranscriptData.R'
load("WGBS/annotationdata/CpGislands.Rdata")
load("WGBS/annotationdata/hg19.tx.Rdata")

#define chromosomes to be used in the transcript tx object
seqlevels(tx, force=TRUE) <- seqlevels(CpGislands, force=TRUE) <- seqlevels(Hsapiens)[1:24]

#read in DMR file generated from '02_01_DMRcalling_plot.R' for annotation and convert into GenomicRanges object
DMRs<-read.csv("WGBS/results/Analysis_DMRs.csv")

DMRsGR<-GRanges(seqnames=DMRs$seqnames,IRanges(DMRs$start,DMRs$end))

#annotate regions for proximity to CpG islands and genic regions using annotatRegions_updown_update function (modified from aaRon)
DMRs2 <- annotateRegions_updown_update(DMRsGR, tx, CpGislands)

# write annotated DMRs out as a csv files for Table 2 and Tables S2 & S4 and later use in Figure 1e
DMRs2 <- as.data.frame(DMRs2)
colnames(DMRs2)[1] <- "chr"
DMRs3<-cbind(DMRs,DMRs2)
write.csv(DMRs3,"WGBS/results/Annotated_Analysis_DMRs.csv")

#read in ChromHMM files for LNCaP and PrEC
#ChromHMM data generated as part of Giles 2019 Epigenetics and Chromatin: doi.org/10.1186/s13072-019-0258-9
#Note, these bed files also uploaded into IGV as tracks to create Figure S1
LNCaP<-read.table("metadata/LNCaP_15_segments.bed",row.names=NULL)
PrEC<-read.table("metadata/PrEC_15_segments.bed",row.names=NULL)

#convert ChromHMM data into GenomicRanges object
LNCaP_GR<-GRanges(seqnames=LNCaP[,1],IRanges(LNCaP[,2],LNCaP[,3]),data=LNCaP[,4])
PrEC_GR<-GRanges(seqnames=PrEC[,1],IRanges(PrEC[,2],PrEC[,3]),data=PrEC[,4])

#Overlap DMRs with LNCaP ChromHMM
OL<-findOverlaps(DMRsGR,LNCaP_GR)
OLfromUniq<-(unique(OL@from))
subset_temp<-list()
    for(j in 1:length(OLfromUniq)){
    subset<-LNCaP_GR[OL@to[which(OL@from==OLfromUniq[j])],]
    cat<-unique(subset@elementMetadata$data)
    cat_temp<-vector()
    for(k in 1:length(cat)){
        cat_temp[k]<-sum(width(subset[which(subset@elementMetadata$data==cat[k]),]))/sum(width(subset))
    }
    names(cat_temp)<-cat
    subset_temp[[j]]<-cat_temp
}
lev<-levels(LNCaP_GR@elementMetadata$data)
levord<-lev[c(1,8:15,2:7)]
mat<-matrix(data=NA,ncol=length(levord),nrow=length(subset_temp))
for(l in 1:length(subset_temp)){
    mat[l,match(names(subset_temp[[l]]),levord)]<-subset_temp[[l]]
}
colnames(mat)<-levord
temp_mat<-cbind(DMRs,mat)

#write out LNCaP ChromHMM overlap for use in Table S4
write.csv(temp_mat,file="WGBS/results/annotatedChromHMMLNCaP_DMRs.csv")

#Overlap DMRs with PrEC ChromHMM
OL<-findOverlaps(DMRsGR,PrEC_GR)
OLfromUniq<-(unique(OL@from))
subset_temp<-list()
for(j in 1:length(OLfromUniq)){
    subset<-PrEC_GR[OL@to[which(OL@from==OLfromUniq[j])],]
    cat<-unique(subset@elementMetadata$data)
    cat_temp<-vector()
    for(k in 1:length(cat)){
      cat_temp[k]<-sum(width(subset[which(subset@elementMetadata$data==cat[k]),]))/sum(width(subset))
    }
    names(cat_temp)<-cat
    subset_temp[[j]]<-cat_temp
  }
  lev<-levels(PrEC_GR@elementMetadata$data)
  levord<-lev[c(1,8:15,2:7)]
  mat<-matrix(data=NA,ncol=length(levord),nrow=length(subset_temp))
  for(l in 1:length(subset_temp)){
    mat[l,match(names(subset_temp[[l]]),levord)]<-subset_temp[[l]]
  }
  colnames(mat)<-levord
  temp_mat<-cbind(DMRs,mat)

#write out PrEC ChromHMM overlap for use in Table S4
  write.csv(temp_mat,file="WGBS/results/annotatedChromHMMPrEC_DMRs.csv")

#Note states are the following
#E1<-"Active TSS"
#E2<-"Flanking active TSS"
#E3<-"Transcription at gene5and3"
#E4<-"Strong transcription"
#E5<-"Weak transcription"
#E6<-"Genic Enhancer"
#E7<-"Enhancer"
#E8<-"ZNF genes and repeats"
#E9<-"Heterochromatin"
#E10<-"Bivalent/Poised TSS"
#E11<-"Flanking Bivalent TSS/Enh"
#E12<-"Bivalent enhancer"
#E13<-"Repressed polycomb"
#E14<-"Weak repressed polycomb"
#E15<-"Quiescent"

