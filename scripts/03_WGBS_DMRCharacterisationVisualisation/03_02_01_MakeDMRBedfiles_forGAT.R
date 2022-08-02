#Script for making bed files of DMR genomic coordinates for use in GAT analysis

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#read in table of DMRs and coordinates from csv file generated in '02_01_DMRcalling_plot.R'
analysis.dmrs <- read.csv("WGBS/results/Analysis_DMRs.csv")

#subset DMRs as to whether they are hyper or hypomethylated in the samples from patients with lethal prostate cancer
hyperDMR<-analysis.dmrs[which(analysis.dmrs$meanbetafc>0),]
hypoDMR<-analysis.dmrs[which(analysis.dmrs$meanbetafc<0),]

hyperDMR2<-hyperDMR[,c("seqnames","start","end")]
hypoDMR2<-hypoDMR[,c("seqnames","start","end")]

#write out coordinates of hyper and hypo DMRs as bed files
write.table(hyperDMR2,file="WGBS/results/hyperDMR_analysis.bed",quote=F,sep="\t",row.names=F,col.names=F)
write.table(hypoDMR2,file="WGBS/results/hypoDMR_analysis.bed",quote=F,sep="\t",row.names=F,col.names=F)

