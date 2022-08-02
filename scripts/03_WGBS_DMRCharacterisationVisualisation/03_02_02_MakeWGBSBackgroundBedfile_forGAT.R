#Script for making bed file of WGBS CpG site genomic coordinates for use as background dataset in GAT analysis

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#read in full results from DSS analysis, generated in '02_01_DMRcalling_plot.R'
load("WGBS/results/DSSresults.RData")

#convert CpG sites from DSS into a data frame and bed file
ra<-results_dml[,1:2]

df <- data.frame(seqnames=ra[,1], start=ra[,2]-1, end=ra[,2], names=c(rep(".", nrow(ra))), scores=c(rep(".", nrow(ra))), strands=c(rep("*", nrow(ra)))) 


write.table(df, file="WGBS/results/AllCpGs_analysis.bed", quote=F, sep="\t", row.names=F, col.names=F)
