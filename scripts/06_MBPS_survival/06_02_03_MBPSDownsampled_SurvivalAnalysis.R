#Script for repeating log rank, univariable Cox and multivariable Cox analysis of CACNA2D4 with downsampled MBPS data
#Generates data for Tables S8 and S9

#load libraries
library(survival) #survival_3.1-8 
library(survminer) #survminer_0.4.6 
library(ggplot2) #ggplot2_3.3.5
library(dplyr) #dplyr_1.0.7
library(Hmisc) #Hmisc_4.3-0 
library(MXM) #MXM_1.4.5

#set user defined working directory and create new folders to store data and output results
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#load processed data generated in '06_02_02_MBPSDownsampled_Processing.R'
load("MBPS/ProcessedData/ProcessedDownsampledMultiplex.Rdata")

#read in validation cohort clinical information (Available upon request. Downloaded into metadata folder in '05_02_PrepareMBPSdata.R')
clin<-read.csv("metadata/Validation_Clinical.csv",header=T)

#prepare clinicopathological data for analysis

#extract clinicopathological data from clinical data table to use in survival analysis

#extracting PSA information
time_PSA<-vector()
psa<-which(clin$PSA_relapse_days!="")
time_PSA[psa]<-clin$PSA_relapse_days[psa]
no_psa<-which(is.na(clin$PSA_relapse_days)==TRUE)
time_PSA[no_psa]<-clin$Time_last_FU[no_psa]

#version of time_PSA in months for KM plots
time_PSAm<-time_PSA/(365/12)

#create PSA_status where 0 = censored (no PSA relapse, lost to follow up), 1 = PSA relapse
PSA_status<-vector()
PSA_status[psa]<-1
PSA_status[no_psa]<-0

#then add new vectors as additional columns in clinical dataframe, clin
clin<-cbind(clin,time_PSA,time_PSAm,PSA_status)

#extracting clinical relapse information - some patients died of PCa but had no clinical relapse - should be included in the group
time_clinical<-vector()
cr<-which(clin$Clinical_relapse_days!="")
time_clinical[cr]<-clin$Clinical_relapse_days[cr]
cr_clin<-which(is.na(clin$Clinical_relapse_days)==TRUE&clin$PPa1.Cause_of_death=="Prostate Cancer") #67,68
time_clinical[cr_clin]<-clin$Time_last_FU[cr_clin]
no_cr<-which(is.na(clin$Clinical_relapse_days)==TRUE&clin$PPa1.Cause_of_death!="Prostate Cancer")
time_clinical[no_cr]<-clin$Time_last_FU[no_cr]

#version of time_clinical in months for KM plots
time_clinicalm<-time_clinical/(365/12)

clinical_status<-vector()
clinical_status[c(cr,cr_clin)]<-1
clinical_status[no_cr]<-0

clin<-cbind(clin,time_clinical,time_clinicalm,clinical_status)

#PCa death only
clin$Death_PCa_status<-ifelse(clin$PPa1.Cause_of_death=="Prostate Cancer",1,0)

#version of time last FU in months for KM plot
time_FUm<-clin$Time_last_FU/(365/12)

clin<-cbind(clin,time_FUm)
# remove sample in row 37 - no PSA level information so uninformative for this analysis

clin<-clin[-37,]

#Dischotamise PSA value
clin$PSA_factor<-ifelse(clin$Pro_PreOpPSAvalue<10,0,1)
clin$PSA_factor<-factor(clin$PSA_factor)

#Dischotamise margin status variable
clin$margin_group<-ifelse(clin$PPr1.Margin_Status=="Not reported" | clin$PPr1.Margin_Status == "Negative",0,1)

#1. Univariable methylation analysis
## Kaplan-Meier plots and log-rank p-vals for univariable methylation analysis (Table S8)####

# Dichotomise CACNA2D4 methylation variable using 75% (third-quartile) cut-off

m<-match(clin$Sample_ID,rownames(M_tab))
M_samp<-M_tab[m,]
M_factor75<-matrix(data=NA,ncol=4,nrow=nrow(M_samp))
for(i in 1:ncol(M_samp)) {
    M_factor75[,i]<-ifelse(M_samp[,i]<=quantile(M_samp[,i],0.75,na.rm=T),0,1)
}

colnames(M_factor75)<-colnames(M_tab)
rownames(M_factor75)<-rownames(M_tab)
clin_75<-cbind(M_factor75,clin)

pval_mat_surv<-matrix(data=NA,ncol=4,nrow=3)
colnames(pval_mat_surv)<-colnames(M_samp)
rownames(pval_mat_surv)<-c("PSA","Clinical","PCSM")

for(i in 1:4){
  ####BCR
  fit1<-survfit(Surv(time_PSAm,PSA_status)~clin_75[,i],data=clin_75)
  
  ####MR
  fit2<-survfit(Surv(time_clinicalm,clinical_status)~clin_75[,i],data=clin_75)
  
  ####PCSM
  fit3<-survfit(Surv(time_FUm,Death_PCa_status)~clin_75[,i],data=clin_75)
  
  pval_mat_surv[1,i]<-surv_pvalue(fit1)$pval.txt
  pval_mat_surv[2,i]<-surv_pvalue(fit2)$pval.txt
  pval_mat_surv[3,i]<-surv_pvalue(fit3)$pval.txt
  
  plot_list=list()
  
  plot_list[[1]]<-ggsurvplot(fit1,data=clin_75,palette=c("mediumblue","firebrick1"),
                             pval=TRUE,xlab="Time to BCR (months)",title=paste("Biochemical recurrence - ",colnames(M_tab)[i],sep=""),
                             break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                             risk.table.y.text.col=T,
                             legend.labs = 
                               c("< 75th percentile", ">= 75th percentile"))
  
  
  plot_list[[2]]<-ggsurvplot(fit2,data=clin_75,palette=c("mediumblue","firebrick1"),
                             pval=TRUE,xlab="Time to clinical relapse (months)",title=paste("Clinical relapse - ",colnames(M_tab)[i],sep=""),
                             break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                             risk.table.y.text.col=T,
                             legend.labs = 
                               c("< 75th percentile", ">= 75th percentile"))
  
  
  plot_list[[3]]<-ggsurvplot(fit3,data=clin_75,palette=c("mediumblue","firebrick1"),
                             pval=TRUE,xlab="Survival time (months)",title=paste("Prostate Cancer mortality - ",colnames(M_tab)[i],sep=""),
                             break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                             risk.table.y.text.col=T,
                             legend.labs = 
                               c("< 75th percentile", ">= 75th percentile"))
  
#export new KM curves with downsample data
  pdf(paste("MBPS/results/colnames(M_tab)[i],"_survival_KM.pdf",sep=""),useDingbats = FALSE)
  for (j in 1:3){
    print(plot_list[[j]])
  }
  dev.off()
}

write.csv(pval_mat_surv,"MBPS/results/survival_KM_pval.csv")

##Cox analysis and c-index for univariable methylation analysis
output_cindexM<-matrix(data=NA,ncol=18,nrow=4)
rownames(output_cindexM) <- colnames(M_tab)
colnames(output_cindexM) <- c("Biochemical relapse: HR", "Biochemical relapse: HR 95%CI", 
                              "Biochemical relapse: HR pval","Biochemical relapse: coxph C-index",
                              "Biochemical relapse: rcorrcens C-index", "Biochemical relapse: C-index pval",
                              "Clinical relapse: HR", "Clinical relapse: HR 95%CI", 
                              "Clinical relapse: HR pval","Clinical relapse: coxph C-index",
                              "Clinical relapse: rcorrcens C-index", "Clinical relapse: C-index pval",
                              "Prostate cancer mortality: HR", "Prostate cancer mortality: HR 95%CI", 
                              "Prostate cancer mortality: HR pval","Prostate cancer mortality: coxph C-index",
                              "Prostate cancer mortality: rcorrcens C-index", "Prostate cancer mortality: C-index pval"
)

for(i in 1:4){
  temp<-i
  surv<-Surv(clin_75$time_PSA,clin_75$PSA_status)
  output_cindexM[i,1]<-round(summary(coxph(surv~clin_75[,temp]))$conf.int[,1],digits=2)
  output_cindexM[i,2]<-paste( round(summary(coxph(surv~clin_75[,temp]))$conf.int[,3],digits=2),"-",
                              round(summary(coxph(surv~clin_75[,temp]))$conf.int[,4], digits=2),sep="")
  output_cindexM[i,3]<-coef(summary(coxph(surv~clin_75[,temp])))[,5]
  output_cindexM[i,4]<-summary(coxph(surv~clin_75[,temp]))$concordance[1]
  output_cindexM[i,5]<-rcorrcens(surv~I(-1*clin_75[,temp]))[1]
  output_cindexM[i,6]<-rcorrcens(surv~I(-1*clin_75[,temp]))[6]
  
  surv2<-Surv(clin_75$time_clinical,clin_75$clinical_status)
  output_cindexM[i,7]<-round(summary(coxph(surv2~clin_75[,temp]))$conf.int[,1],digits=2)
  output_cindexM[i,8]<-paste( round(summary(coxph(surv2~clin_75[,temp]))$conf.int[,3],digits=2),"-",
                              round(summary(coxph(surv2~clin_75[,temp]))$conf.int[,4], digits=2),sep="")
  output_cindexM[i,9]<-coef(summary(coxph(surv2~clin_75[,temp])))[,5]
  output_cindexM[i,10]<-summary(coxph(surv2~clin_75[,temp]))$concordance[1]
  output_cindexM[i,11]<-rcorrcens(surv2~I(-1*clin_75[,temp]))[1]
  output_cindexM[i,12]<-rcorrcens(surv2~I(-1*clin_75[,temp]))[6]
  
  surv3<-Surv(clin_75$Time_last_FU,clin_75$Death_PCa_status)
  output_cindexM[i,13]<-round(summary(coxph(surv3~clin_75[,temp]))$conf.int[,1],digits=2)
  output_cindexM[i,14]<-paste( round(summary(coxph(surv3~clin_75[,temp]))$conf.int[,3],digits=2),"-",
                               round(summary(coxph(surv3~clin_75[,temp]))$conf.int[,4], digits=2),sep="")
  output_cindexM[i,15]<-coef(summary(coxph(surv3~clin_75[,temp])))[,5]
  output_cindexM[i,16]<-summary(coxph(surv3~clin_75[,temp]))$concordance[1]
  output_cindexM[i,17]<-rcorrcens(surv3~I(-1*clin_75[,temp]))[1]
  output_cindexM[i,18]<-rcorrcens(surv3~I(-1*clin_75[,temp]))[6]
}

#write out output_cindex for second part of Table S8
write.csv(output_cindexM, file ="MBPS/results/Methylation_univariate_cox_cindex.csv")


# Multivariable analysis with downsampled CACNA2D4

####Multivariable clinopathological and methylation analysis
#Using SES for model selection and then extracting results for manuscript (using CACNA2D4 only as this is the only downsampled methylation data)
covnames<-colnames(clin_75)[1:4]
for(i in 1:length(covnames)){
  varInt<-c(covnames[i],"ISUP_groups","Stage_groups","margin_group","PSA_factor")
  data_75split<-clin_75[,varInt]
  ses_psa_75<- SES(surv,data_75split,test="censIndCR")
  ses_psa_75.model<-ses.model(surv,data_75split,test="censIndCR",sesObject=ses_psa_75)
  lp.ses_psa_75.model<-predict(ses_psa_75.model$mod,predict="lp")
  
  ses_cr_75<- SES(surv2,data_75split,test="censIndCR")
  ses_cr_75.model<-ses.model(surv2,data_75split,test="censIndCR",sesObject=ses_cr_75)
  lp.ses_cr_75.model<-predict(ses_cr_75.model$mod,predict="lp")

  ses_death_75<- SES(surv3,data_75split,test="censIndCR")
  ses_death_75.model<-ses.model(surv3,data_75split,test="censIndCR",sesObject=ses_death_75)
  lp.ses_death_75.model<-predict(ses_death_75.model$mod,predict="lp")

  print(names((ses_psa_75.model[[1]])$assign))
  print(names((ses_cr_75.model[[1]])$assign))
  print(names((ses_death_75.model[[1]])$assign))
  
  output_cindex_SES_75split<-matrix(data=NA,nrow=3,ncol=6)
  rownames(output_cindex_SES_75split) <- c("CACNA2D4 + Grade + PSA (Biochemical recurrence)",
                                         "Grade + Margin Status (Clinical relapse)",
                                         "CACNA2D4 + Grade (Prostate Cancer Mortality)")
  colnames(output_cindex_SES_75split) <- c("HR", "HR 95%CI", 
                                         "HR pval","coxph C-index",
                                         "rcorrcens C-index", "C-index pval")

  output_cindex_SES_75split[1,1]<-paste(round(summary(ses_psa_75.model$mod)$conf.int[1,1],digits=2),";",round(summary(ses_psa_75.model$mod)$conf.int[2,1],digits=2),";",round(summary(ses_psa_75.model$mod)$conf.int[3,1],digits=2),sep="")
  output_cindex_SES_75split[1,2]<-paste(round(summary(ses_psa_75.model$mod)$conf.int[1,3],digits=2),"-",
                                      round(summary(ses_psa_75.model$mod)$conf.int[1,4], digits=2),";",
                                      round(summary(ses_psa_75.model$mod)$conf.int[2,3],digits=2),"-",
                                      round(summary(ses_psa_75.model$mod)$conf.int[2,4], digits=2),";",
                                      round(summary(ses_psa_75.model$mod)$conf.int[3,3],digits=2),"-",
                                      round(summary(ses_psa_75.model$mod)$conf.int[3,4], digits=2),sep="")
  output_cindex_SES_75split[1,3]<-paste(round(coef(summary(ses_psa_75.model$mod))[1,5],digits=5),";",round(coef(summary(ses_psa_75.model$mod))[2,5],digits=5),";",round(coef(summary(ses_psa_75.model$mod))[3,5],digits=5),sep="")
  output_cindex_SES_75split[1,4]<-summary(ses_psa_75.model$mod)$concordance[1]
  output_cindex_SES_75split[1,5]<-rcorrcens(surv~I(-1*lp.ses_psa_75.model))[1]
  output_cindex_SES_75split[1,6]<-rcorrcens(surv~I(-1*lp.ses_psa_75.model))[6]

  output_cindex_SES_75split[2,1]<-paste(round(summary(ses_cr_75.model$mod)$conf.int[1,1],digits=2),";",round(summary(ses_cr_75.model$mod)$conf.int[2,1],digits=2),sep="")
  output_cindex_SES_75split[2,2]<-paste(round(summary(ses_cr_75.model$mod)$conf.int[1,3],digits=2),"-",
                                      round(summary(ses_cr_75.model$mod)$conf.int[1,4], digits=2),";",
                                      round(summary(ses_cr_75.model$mod)$conf.int[2,3],digits=2),"-",
                                      round(summary(ses_cr_75.model$mod)$conf.int[2,4], digits=2),sep="")
  output_cindex_SES_75split[2,3]<-paste(round(coef(summary(ses_cr_75.model$mod))[1,5],digits=5),";",round(coef(summary(ses_cr_75.model$mod))[2,5],digits=5),sep="")
  output_cindex_SES_75split[2,4]<-summary(ses_cr_75.model$mod)$concordance[1]
  output_cindex_SES_75split[2,5]<-rcorrcens(surv2~I(-1*lp.ses_cr_75.model))[1]
  output_cindex_SES_75split[2,6]<-rcorrcens(surv2~I(-1*lp.ses_cr_75.model))[6]

  output_cindex_SES_75split[3,1]<-paste(round(summary(ses_death_75.model$mod)$conf.int[1,1],digits=2),";",round(summary(ses_death_75.model$mod)$conf.int[2,1],digits=2),sep="")
  output_cindex_SES_75split[3,2]<-paste(round(summary(ses_death_75.model$mod)$conf.int[1,3],digits=2),"-",
                                      round(summary(ses_death_75.model$mod)$conf.int[1,4], digits=2),";",
                                      round(summary(ses_death_75.model$mod)$conf.int[2,3],digits=2),"-",
                                      round(summary(ses_death_75.model$mod)$conf.int[2,4], digits=2),sep="")
  output_cindex_SES_75split[3,3]<-paste(round(coef(summary(ses_death_75.model$mod))[1,5],digits=5),";",round(coef(summary(ses_death_75.model$mod))[2,5],digits=5),sep="")
  output_cindex_SES_75split[3,4]<-summary(ses_death_75.model$mod)$concordance[1]
  output_cindex_SES_75split[3,5]<-rcorrcens(surv3~I(-1*lp.ses_death_75.model))[1]
  output_cindex_SES_75split[3,6]<-rcorrcens(surv3~I(-1*lp.ses_death_75.model))[6]

  #write out results of SES model on downsampled data - Table S9
  write.csv(output_cindex_SES_75split, file =paste("MBPS/results/Methylation_multivariate_cox_cindex_",covnames[i],".csv",sep=""))
}
