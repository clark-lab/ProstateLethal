#Script for extracting QC data for MBPS data
#Generates summary statistics for Table S5

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#MBPS QC data generated from processing multiplex sequencing reads through methpanel https://github.com/thinhong/MethPanel

#read in Metrics table that were exported by methpanel
mettab1<-read.table("MBPS/data/metrics_table1.tsv",header=T)
met1<-read.table("MBPS/data/metrics1.tsv",header=T)
mettab2<-read.table("MBPS/data/metrics_table2.tsv",header=T)
met2<-read.table("MBPS/data/metrics2.tsv",header=T)

#Total number of reads - run 1
sum(met1$Number_Mapped_Read)
#201402566

#Total number of reads - run 2
sum(met2$Number_Mapped_Read)
#205189303

#Total number of reads after QC - run 1
sum(met1$Number_Mapped_Converted_40mapq_Read)
#175163574

#Total number of reads after QC - run 2
sum(met2$Number_Mapped_Converted_40mapq_Read)
#179587048

#Average Coverage of CpGs - run 1
mean(mettab1$Depth.of.CpG.with.full.read..X.)
#100161

#Average Coverage of CpGs - run 2
mean(mettab2$Depth.of.CpG.with.full.read..X.)
#102722.6


