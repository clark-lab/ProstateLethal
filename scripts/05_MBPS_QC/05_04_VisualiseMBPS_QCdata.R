#Script for visualising MBPS QC data
#Generates plots for Figures S6, S7 & S8

#set user defined working directory and create new folder to output results
BASEDIR<-*user_defined_file_path*
setwd(paste(BASEDIR,"/MBPS",sep=""))
dir.create("plots")
setwd(BASEDIR)

#load libraries
library(ggplot2) #ggplot2_3.3.5
library(corrplot) #corrplot_0.84 
library(limma) #limma_3.42.0 
library(tidyverse) #tidyverse_1.3.1
library(gridExtra) #gridExtra_2.3

#load data

#load in mean multiplex methylation data, prepared in '05_02_PrepareMBPSdata.R'
load("MBPS/results/ProcessedMultiplexMeans.RData")

#avM_samp contains methylation percent values for each amplicon and sample

#read in validation cohort clinical information (Available upon request. Downloaded into metadata folder in '05_02_PrepareMBPSdata.R')
clin<-read.csv("metadata/Validation_Clinical.csv",header=T)

##create MDS plot comparing 2 sequencing runs (Figure S6)
plotMDS1<-plotMDS(avM_samp, top = nrow(avM_samp), cex = 0.5, main = "All samples", labels = colnames(avM_samp))
toplot <- data.frame(distance.matrix = plotMDS1$distance.matrix, x = plotMDS1$x, y = plotMDS1$y, PCR = factor(clin$Batch))
MDSplot <- ggplot(toplot$distance.matrix, aes(x = toplot$x, y = toplot$y, colour=as.factor(clin$Batch))) + ggtitle("MDS All Samples (All amplicons)") +
  geom_point(size =4, alpha=0.8) +  theme_classic() + ylab("PC2") + xlab("PC1") + theme(plot.title = element_text(hjust = 0.5))
  
#export MDS plot (Figure S6)
pdf("MBPS/plots/MDS_batch_QC.pdf", width = 9)
print(MDSplot)
dev.off()


#load in multiplex methylation data, prepared in '05_02_PrepareMBPSdata.R'
load("MBPS/results/ProcessedMultiplex.RData")

##create correlation plot (Figure S8)
#transform avM_samp so that the correlation is between amplicons not samples ...
tavM_samp <- t(avM_samp)
tavM_samp <- as.data.frame(tavM_samp)

#remove amplicon that failed completely from namelist (list of amplicons) and order namelist to how amplicons ordered in manuscript
namelisttrim<-namelist[-7]

namelist_ordered<-c("CD34","AC074091.13","EPHB3","PRDM8","Y_RNA","CDO1","ZNF655",
                    "CCDC67","CACNA2D4","EFCAB4B","PAH","LECT1","CALM1",
                    "PARP6_CELF6","RP11-927P21.4","TBX1","BHLHB9")

a<-match(namelist_ordered,namelisttrim) 

#remove amplicon that failed completely from methylation data
tavM_samp1 <- tavM_samp[,-7]

#reformat methylation data frame so that amplicons are in same order and numbered as in other parts of manuscript
tavM_samp2<-tavM_samp1[,a]

colnames(tavM_samp2)<- factor(c(1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18))

#Perform correlation between amplicons for all samples using spearman correlation
#spearman rho
cor_samp_spearman2 <- cor(tavM_samp2,use = "complete.obs", method = "spearman")

#p-values and confidence intervals
res_spearman2 <- cor.mtest(tavM_samp2, method= "spearman", conf.level = .95)

##Wxport correlation plot (Figure S8)
pdf("MBPS/plots/correlation_amplicons_QC.pdf",width=7,useDingbats = FALSE)
corrplot(cor_samp_spearman2, method="number",type="upper", order="hclust",diag=FALSE, p.mat = res_spearman2$p, sig.level=0.05, title = "Correlation between lethal amplicons",insig="blank", mar=c(0,0,1,0))
dev.off()

#create boxplots of sequencing coverage (Figure S7)

#remove amplicon that failed completely from coverage data
meancovtab<-meancovtab[-7,]
meancovtab1<-meancovtab1[-7,]

#reformat coverage data frame so that amplicons are in same order and numbered as in other parts of manuscript

meancovtab<-meancovtab[a,]
meancovtab1<-meancovtab1[a,]

#combine coverage dataframes from different sequencing runs
meancovtab_merged<-cbind(meancovtab,meancovtab1)

meancov <- as.data.frame(meancovtab_merged)

#reformat coverage data frame to long format with amplicons numbered as in other parts of manuscript
meancov$num <- factor(c(1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18))
cov_long<-gather(meancov,"sample","coverage",-num)

plot5 <-ggplot(data=cov_long, aes(x=num, y=coverage, group = num)) +
  stat_boxplot(geom = "errorbar", width = 0.2, lwd =0.4) +
  geom_boxplot(fill='cyan3', lwd=0.25, outlier.shape = 21, outlier.size = 1) +
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=0.5))+
  xlab("DMR")+
  ylab("Coverage")+
  scale_y_continuous(labels = scales::comma)+ #change from scientific notation
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

#change the y-axis to amplify the lower coverage boxplots
plot6 <-ggplot(data=cov_long, aes(x=num, y=coverage, group = num)) +
  stat_boxplot(geom = "errorbar", width = 0.2, lwd =0.4) +
  geom_boxplot(fill='cyan3', lwd=0.25, outlier.shape = 21, outlier.size = 1) +
  ylim(0,100000)+
  geom_hline(yintercept=100, linetype="dashed", color = "gray")+
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=0.5))+
  xlab("DMR")+
  ylab("Coverage")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))


#Export Figure S7
pdf("MBPS/plots/coverageboxplots_QC.pdf",width=15,useDingbats=FALSE)
grid.arrange(plot5, plot6, nrow=2, ncol=1)
print(plot5)
print(plot6)
dev.off()







