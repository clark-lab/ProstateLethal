#Script for performing gene ontology analysis on DMR associated genes
#Generating data for Table S3 and Figure 1g

#load libraries
library(GenomicRanges)
library(RITAN)
library(RITANdata)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(liftOver)
library(ggplot2)
library(DOSE)

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#read in Genehancer Interactions Double Elite regions (originally downloaded from UCSC Table Browser)
GH_interactions_elite <- read.table("metadata/GH_Interactions_Double_Elite_hg38.txt", sep="\t")

#Liftover Genehancer Regions from hg38 to hg19
interactionsgr.hg38 <- GRanges(GH_interactions_elite$V1, IRanges(GH_interactions_elite$V2, GH_interactions_elite$V3), mcols=GH_interactions_elite[,4:18])
genome(interactionsgr.hg38) <- "hg38"

#chain for lift over originally downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz
ch <- import.chain("metadata/hg38ToHg19.over.chain")
seqlevelsStyle(interactionsgr.hg38) <-  "UCSC"  

interactionsgr.hg19 <- unlist(liftOver(interactionsgr.hg38, ch))

#define genhance background gene set
background <- as.character(unique(interactionsgr.hg19$mcols.V17)) 

#read in DMR table generated in '02_01_DMRcalling_plot.R'
analysis.dmrs <- read.csv("WGBS/results/Analysis_DMRs.csv")
analysis.ranges <- GRanges(analysis.dmrs$seqnames, 
                            IRanges(analysis.dmrs$start, analysis.dmrs$end),
                            mcols=analysis.dmrs[,5:12])
ol <- findOverlaps(analysis.ranges, interactionsgr.hg19)
getcandgenes <- function(dmr){
  cands <- as.character(unique(interactionsgr.hg19$mcols.V17[subjectHits(ol)[queryHits(ol)==dmr]]))
  paste(cands, collapse=", ")
}
analysis.ranges$candgenes <- sapply(1:length(analysis.ranges), getcandgenes)

analysis.ranges_up <- analysis.ranges[analysis.ranges$mcols.meanbetafc > 0]
analysis_up <- unique(unlist(strsplit(paste(analysis.ranges_up$candgenes, collapse=", "), ", ")))
analysis_up <- analysis_up[grep("^ENSG", analysis_up, invert = T)]
analysis_up <- analysis_up[nchar(analysis_up) > 0]

analysis_up_res <- term_enrichment(analysis_up, resources=names(geneset_list), all_symbols = background,report_resources_separately = T)
analysis_up_res <- lapply(analysis_up_res, function (x) x[x$q < 0.05,])
analysis_up_res <- plyr::rbind.fill(analysis_up_res)
analysis_up_res <- analysis_up_res[order(analysis_up_res$q),]
write.csv(analysis_up_res, file="WGBS/results/RITAN_GSE_analysis_GH_up_withbkg.csv")

analysis.ranges_down <- analysis.ranges[analysis.ranges$mcols.maxbetafc < 0]
analysis_down <- unique(unlist(strsplit(paste(analysis.ranges_down$candgenes, collapse=", "), ", ")))
analysis_down <- analysis_down[grep("^ENSG", analysis_down, invert = T)]
analysis_down <- analysis_down[nchar(analysis_down) > 0]

analysis_down_res <- term_enrichment(analysis_down, resources=names(geneset_list), all_symbols = background,report_resources_separately = T)
analysis_down_res <- lapply(analysis_down_res, function (x) x[x$q < 0.05,])
analysis_down_res <- plyr::rbind.fill(analysis_down_res)
analysis_down_res <- analysis_down_res[order(analysis_down_res$q),]
write.csv(analysis_down_res, file="WGBS/results/RITAN_GSE_analysis_GH_down_withbkg.csv")

up2<-analysis_up_res[which(analysis_up_res$n.set>15 & analysis_up_res$n.set<500 ),]
write.csv(up2,"WGBS/results/HyperDMR_withbkg/HyperFilteredbySize.csv")


#separate RITAN results into different categories
listcat<-c("MSigDB_Hallmarks", "MSigDB_C2", "MSigDB_C3.TFtargets", "MSigDB_C3.miRNA", "MSigDB_C5", "MSigDB_C7", "GO[.]", "GO_slim_PIR", "GO_slim_generic", "PathwayCommonsPathways", "ReactomePathways", "NetPath_Gene_regulation", "Chaussabel_Modules", "Blood_Translation_Modules", "KEGG_filtered_canonical_pathways", "DisGeNet")


#write out up2 with categories separated
cat<-list()
for(i in 1:length(listcat)){
  cat[[i]]<-grep(listcat[i],up2$name )
}
catvec<-vector()
for(i in 1:length(cat)){
  catvec[cat[[i]]]<-listcat[i]
}

up_cat<-cbind(catvec,up2)

namevec<-vector()
for(j in 1:length(listcat)){
  temp<-grep(listcat[j],up_cat$catvec)
  for(i in 1:length(temp)){
    namevec[temp[i]]<-(gsub(listcat[j],"",up_cat$name[temp[i]]))
  }
}

a<-grep("GO.",up_cat$catvec)
b<-grep("GO_slim",up_cat$catvec)
a2<-a[-match(b,a)]
for(i in 1:length(a2)){
  namevec[a2[i]]<-(gsub(listcat[7],"",up_cat$name[a2[i]]))
}

up3<-cbind(up_cat,namevec)

#write out RITAN results with categories as separate column - data for Table S3
write.csv(up3,"WGBS/results/HyperFilteredbySize_separated.csv")

#work out how many significant terms in each category and create tables for each category
subset<-list()
for(i in 1:length(listcat)){
  subset[[i]]<-up2[grep(listcat[i],up2$name ),]
}
names(subset)<-listcat

sum(unlist(lapply(subset,nrow))) #298 (have all been divided appropriately)

unlist(lapply(subset,nrow))

#create bar plot of top 10 results from MSigDB - for Figure 1g
pdf("WGBS/plots/EnrichmentBarplot_MSigDB_C2_top10.pdf",width=12)
i<-2
sub<-subset[[i]][order(subset[[i]]$q),]
sub$name<-gsub("MSigDB_C2.","",sub$name)
sub$name<-factor(sub$name,levels=rev(sub$name))
sub<-sub[1:10,]
p <- ggplot(sub, aes(x=name, y = n, fill = q)) +
  theme_dose(font.size=8) +
  scale_fill_continuous(low="red", high="blue", name = "q", guide=guide_colorbar(reverse=TRUE))
print(p + geom_bar(stat = "identity") + coord_flip() +
        ggtitle(names(subset)[i]) + xlab(NULL) + ylab(NULL))
print(p + geom_bar(stat = "identity") +
        ggtitle(names(subset)[i]) + xlab(NULL) + ylab(NULL))+
  theme(axis.text.x = element_text(angle = 90))
dev.off()

