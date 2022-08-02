#Script for creating boxplot of LNCaP and PrEC WGBS data at 18 selected DMRs
#Generates Figure S4c

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#read in LNCaP PreC WGBS methylation data for all DMRs, generated in '03_03_DMRMeth_LNCaPPrEC.R'
Meth<-read.csv("/WGBS/results/LNCaP_PreC_DMRmeth_1420.csv")

#read in DMRs that were assayed in multiplex
selecDMR<-read.table("metadata/SelectedRegions.bed")

#extract LNCaP PrEC data for 18 DMRs only
M3<-Meth[match(selecDMR$V2,Meth$start),]

# rearrange to DMR order according to genomic positions:
M3 <- M3[c(10,11,2,7,18,17,15,3,13,6,16,4,8,1,14,9,12,5),]

# subset  methylation data to just LNCaP and PrEC samples
M_data <- as.matrix(M3[,c("meanL...100", "meanP...100")])

tM_data <- t(M_data)

#export boxplot of LNCaP and PrEC WGBS data at 18 DMRs - Figure S4c
a<-1:18

pdf("WGBS/plots/LNCaP_PreC_DMRmeth.pdf",width=15,height=5,useDingbats=FALSE)
barplot(tM_data, beside=TRUE,
        xlab = "DMR",ylab = "DNA methylation (%)", 
        ylim=c(0,101),las =1,
        names.arg=a,width=5,
        col = c("purple","lightblue"), axis.lty=1,
        )
box()
dev.off()


