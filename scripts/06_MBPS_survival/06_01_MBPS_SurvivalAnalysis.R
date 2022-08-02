#Script for performing survival analysis of MBPS data

#Generates plots for Figures 2 & 3 and Figures S10 & S9
#Generates data for Table S6 & S7

#load libraries
library(survival) #survival_3.1-8
library(survminer) #survminer_0.4.6
library(ggplot2) #ggplot2_3.3.5
library(dplyr) #dplyr_1.0.7
library(Hmisc) #Hmisc_4.3-0
library(MXM) #MXM_1.4.5
library(forestplot)
library(risksetROC)

#set user defined working directory and create new folders to output results
BASEDIR<-*user_defined_file_path*
setwd(paste(BASEDIR,"/MBPS/plots/",sep=""))
dir.create("survivalanalysis")
setwd(BASEDIR)


#load data

#read in validation cohort clinical information (Available upon request. Downloaded into metadata folder in '05_02_PrepareMBPSdata.R')
clin<-read.csv("metadata/Validation_Clinical.csv",header=T)

#read in multiplex data, generated from '05_02_PrepareMBPSdata.R'
load("MBPS/results/ProcessedMultiplexMeans.RData")
load("MBPS/results/ProcessedMultiplex.RData")
#contains M, genes, namelist, R_pos, meancovtab, meancovtab1

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
# remove sample in row 37 as no PSA level information so uninformative for this analysis

clin<-clin[-37,]

#1. Univariable clinicopathological analysis

##Kaplan-Meier plots and log-rank p-vals for univariable clinicopathological analysis (Figure 2)

###a. Dichotomised ISUP groups (Gleason 3+4 / ISUP 2 vs Gleason 4+3,8,9 / ISUP 3,4,5)

####BCR
plot_list1=list()

fit1<-survfit(Surv(time_PSAm,PSA_status)~ISUP_groups,data=clin)
plot_list1[[1]]<-ggsurvplot(fit1,data=clin,palette=c("mediumblue","firebrick1"),
                            pval=TRUE,xlab="Time to BCR (months)",title="Biochemical recurrence - Gleason/ISUP grade",
                            break.time.by=50, ggtheme=theme_classic(),risk.table= T,fontsize=3,#tables.theme=clean_theme(),
                            risk.table.y.text.col=T,
                            legend.labs = 
                              c("Group 1", "Group 2"))
####MR
fit2<-survfit(Surv(time_clinicalm,clinical_status)~ISUP_groups,data=clin)
plot_list1[[2]]<-ggsurvplot(fit2,data=clin,palette=c("mediumblue","firebrick1"),
                            pval=TRUE,xlab="Time to clinical relapse (months)",title="Clinical relapse - Gleason/ISUP grade",
                            break.time.by=50, ggtheme=theme_classic(),risk.table= T,fontsize=3,#tables.theme=clean_theme(),
                            risk.table.y.text.col=T,
                            legend.labs = 
                              c("Group 1", "Group 2"))
####PCSM
fit3<-survfit(Surv(time_FUm,Death_PCa_status)~ISUP_groups,data=clin)
plot_list1[[3]]<-ggsurvplot(fit3,data=clin,palette=c("mediumblue","firebrick1"),
                            pval=TRUE,xlab="Survival time (months)",title="Prostate cancer mortality - Gleason/ISUP grade",
                            break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                            risk.table.y.text.col=T,
                            legend.labs = 
                              c("Group1", "Group 2"))
#export plots for Figure 2
pdf("MBPS/plots/survivalanalysis/Gleasongrades_survival_KM.pdf",useDingbats = FALSE)
for (j in 1:3){
  print(plot_list1[[j]])
}
dev.off()



### b. Dichotomised pathological stage (<=pT2 vs >=pT3)

####BCR
plot_list2=list()

fit1<-survfit(Surv(time_PSAm,PSA_status)~Stage_groups,data=clin)
plot_list2[[1]]<-ggsurvplot(fit1,data=clin,palette=c("mediumblue","firebrick1"),
                            pval=TRUE,xlab="Time to BCR (months)",title="Biochemical recurrence - Clinical stage",
                            break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                            risk.table.y.text.col=T,
                            legend.labs = 
                              c("Group 1", "Group 2"))

####MR
fit2<-survfit(Surv(time_clinicalm,clinical_status)~Stage_groups,data=clin)
plot_list2[[2]]<-ggsurvplot(fit2,data=clin,palette=c("mediumblue","firebrick1"),
                            pval=TRUE,xlab="Time to clinical relapse (months)",title="Clinical relapse - Clinical stage",
                            break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                            risk.table.y.text.col=T,
                            legend.labs = 
                              c("Group 1", "Group 2"))

####PCSM
fit3<-survfit(Surv(time_FUm,Death_PCa_status)~Stage_groups,data=clin)
plot_list2[[3]]<-ggsurvplot(fit3,data=clin,palette=c("mediumblue","firebrick1"),
                            pval=TRUE,xlab="Survival time (months)",title="Prostate cancer mortality - Clinical stage",
                            break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                            risk.table.y.text.col=T,
                            legend.labs = 
                              c("Group 1", "Group 2"))
#export plots for Figure 2
pdf("MBPS/plots/survivalanalysis/Pathologicalstage_survival_KM.pdf",useDingbats=FALSE)
for (j in 1:3){
  print(plot_list2[[j]])
}
dev.off()


### c. Dichotamised PSA level (<10 vs >10ng/mL)

clin$PSA_factor<-ifelse(clin$Pro_PreOpPSAvalue<10,0,1)
clin$PSA_factor<-factor(clin$PSA_factor)

####BCR
plot_list3=list()

fit1<-survfit(Surv(time_PSAm,PSA_status)~PSA_factor,data=clin)
plot_list3[[1]]<-ggsurvplot(fit1,data=clin,palette=c("mediumblue","firebrick1"),
                            pval=TRUE,xlab="Time to BCR (months)",title="Biochemical recurrence - PSA levels",
                            break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                            risk.table.y.text.col=T,
                            legend.labs = 
                              c("< 10ng/ml", ">= 10ng/ml"))

####MR
fit2<-survfit(Surv(time_clinicalm,clinical_status)~PSA_factor,data=clin)
plot_list3[[2]]<-ggsurvplot(fit2,data=clin,palette=c("mediumblue","firebrick1"),
                            pval=TRUE,xlab="Time to clinical relapse (months)",title="Clinical relapse - PSA levels",
                            break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                            risk.table.y.text.col=T,
                            legend.labs = 
                              c("< 10ng/ml", ">= 10ng/ml"))

####PCSM
fit3<-survfit(Surv(time_FUm,Death_PCa_status)~PSA_factor,data=clin)
plot_list3[[3]]<-ggsurvplot(fit3,data=clin,palette=c("mediumblue","firebrick1"),
                            pval=TRUE,xlab="Survival time (months)",title="Prostate cancer mortality - PSA levels",
                            break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                            risk.table.y.text.col=T,
                            legend.labs = 
                              c("< 10ng/ml", ">= 10ng/ml"))
#export plots for Figure 2
pdf("MBPS/plots/survivalanalysis/PSAlevels_survival_KM.pdf",useDingbats = FALSE)
for (j in 1:3){
  print(plot_list3[[j]])
}
dev.off()


### d. Dichotamised Margin Status (Negative vs Positive)

clin$margin_group<-ifelse(clin$PPr1.Margin_Status=="Not reported" | clin$PPr1.Margin_Status == "Negative",0,1)

# only 1 with "not reported"

####BCR
plot_list4=list()

fit1<-survfit(Surv(time_PSAm,PSA_status)~margin_group,data=clin)
plot_list4[[1]]<-ggsurvplot(fit1,data=clin,palette=c("mediumblue","firebrick1"),
                            pval=TRUE,xlab="Time to BCR (months)",title="Biochemical recurrence - Margin status",
                            break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                            risk.table.y.text.col=T,
                            legend.labs = 
                              c("Negative", "Positive"))

####MR
fit2<-survfit(Surv(time_clinicalm,clinical_status)~margin_group,data=clin)
plot_list4[[2]]<-ggsurvplot(fit2,data=clin,palette=c("mediumblue","firebrick1"),
                            pval=TRUE,xlab="Time to clinical relapse (months)",title="Clinical relapse - Margin status",
                            break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                            risk.table.y.text.col=T,
                            legend.labs = 
                              c("Negative", "Positive"))

####PCSM
fit3<-survfit(Surv(time_FUm,Death_PCa_status)~margin_group,data=clin)
plot_list4[[3]]<-ggsurvplot(fit3,data=clin,palette=c("mediumblue","firebrick1"),
                            pval=TRUE,xlab="Survival time (months)",title="Prostate cancer mortality - Margin status",
                            break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                            risk.table.y.text.col=T,
                            legend.labs = 
                              c("Negative", "Positive"))

#export plots for Figure 2
pdf("MBPS/plots/survivalanalysis/marginstatus_survival_KM.pdf",useDingbats = FALSE)
for (j in 1:3){
  print(plot_list4[[j]])
}
dev.off()


##Cox analysis and c-index for univariable clinopathological analysis (Table S6)

#create empty matrix for results
output_cindex<-matrix(data=NA,nrow=4,ncol=18)
rownames(output_cindex) <- c("Gleason/ISUP grade groups", "Clinical Stage groups",
                             "PSA levels", "Margin status")
colnames(output_cindex) <- c("Biochemical relapse: HR", "Biochemical relapse: HR 95%CI",
                             "Biochemical relapse: HR pval","Biochemical relapse: coxph C-index",
                             "Biochemical relapse: rcorrcens C-index", "Biochemical relapse: C-index pval",
                             "Clinical relapse: HR", "Clinical relapse: HR 95%CI",
                             "Clinical relapse: HR pval","Clinical relapse: coxph C-index",
                             "Clinical relapse: rcorrcens C-index", "Clinical relapse: C-index pval",
                             "Prostate cancer mortality: HR", "Prostate cancer mortality: HR 95%CI",
                             "Prostate cancer mortality: HR pval","Prostate cancer mortality: coxph C-index",
                             "Prostate cancer mortality: rcorrcens C-index", "Prostate cancer mortality: C-index pval"
)

###a. Dichotomised ISUP groups (Gleason 3+4 / ISUP 2 vs Gleason 4+3,8,9 / ISUP 3,4,5)
output_cindex[1,1]<-round(summary(coxph(Surv(time_PSA,PSA_status)~ISUP_groups,data=clin))$conf.int[,1],digits=2)
output_cindex[1,2]<-paste( round(summary(coxph(Surv(time_PSA,PSA_status)~ISUP_groups,data=clin))$conf.int[,3],digits=2),"-",
                           round(summary(coxph(Surv(time_PSA,PSA_status)~ISUP_groups,data=clin))$conf.int[,4], digits=2),sep="")
output_cindex[1,3]<-coef(summary(coxph(Surv(time_PSA,PSA_status)~ISUP_groups,data=clin)))[,5]
output_cindex[1,4]<-summary(coxph(Surv(time_PSA,PSA_status)~ISUP_groups,data=clin))$concordance[1]
output_cindex[1,5]<-rcorrcens(Surv(time_PSA,PSA_status)~I(-1*as.numeric(ISUP_groups)),data=clin)[1]
output_cindex[1,6]<-rcorrcens(Surv(time_PSA,PSA_status)~I(-1*as.numeric(ISUP_groups)),data=clin)[6]

output_cindex[1,7]<-round(summary(coxph(Surv(time_clinical,clinical_status)~ISUP_groups,data=clin))$conf.int[,1],digits=2)
output_cindex[1,8]<-paste( round(summary(coxph(Surv(time_clinical,clinical_status)~ISUP_groups,data=clin))$conf.int[,3],digits=2),"-",
                           round(summary(coxph(Surv(time_clinical,clinical_status)~ISUP_groups,data=clin))$conf.int[,4], digits=2),sep="")
output_cindex[1,9]<-coef(summary(coxph(Surv(time_clinical,clinical_status)~ISUP_groups,data=clin)))[,5]
output_cindex[1,10]<-summary(coxph(Surv(time_clinical,clinical_status)~ISUP_groups,data=clin))$concordance[1]
output_cindex[1,11]<-rcorrcens(Surv(time_clinical,clinical_status)~I(-1*as.numeric(ISUP_groups)),data=clin)[1]
output_cindex[1,12]<-rcorrcens(Surv(time_clinical,clinical_status)~I(-1*as.numeric(ISUP_groups)),data=clin)[6]

output_cindex[1,13]<-round(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~ISUP_groups,data=clin))$conf.int[,1],digits=2)
output_cindex[1,14]<-paste( round(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~ISUP_groups,data=clin))$conf.int[,3],digits=2),"-",
                            round(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~ISUP_groups,data=clin))$conf.int[,4], digits=2),sep="")
output_cindex[1,15]<-coef(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~ISUP_groups,data=clin)))[,5]
output_cindex[1,16]<-summary(coxph(Surv(Time_last_FU,Death_PCa_status)~ISUP_groups,data=clin))$concordance[1]
output_cindex[1,17]<-rcorrcens(Surv(Time_last_FU,Death_PCa_status)~I(-1*as.numeric(ISUP_groups)),data=clin)[1]
output_cindex[1,18]<-rcorrcens(Surv(Time_last_FU,Death_PCa_status)~I(-1*as.numeric(ISUP_groups)),data=clin)[6]

### b. Dichotomised pathological stage (<=pT2 vs >=pT3)
output_cindex[2,1]<-round(summary(coxph(Surv(time_PSA,PSA_status)~Stage_groups,data=clin))$conf.int[,1],digits=2)
output_cindex[2,2]<-paste( round(summary(coxph(Surv(time_PSA,PSA_status)~Stage_groups,data=clin))$conf.int[,3],digits=2),"-",
                           round(summary(coxph(Surv(time_PSA,PSA_status)~Stage_groups,data=clin))$conf.int[,4], digits=2),sep="")
output_cindex[2,3]<-coef(summary(coxph(Surv(time_PSA,PSA_status)~Stage_groups,data=clin)))[,5]
output_cindex[2,4]<-summary(coxph(Surv(time_PSA,PSA_status)~Stage_groups,data=clin))$concordance[1]
output_cindex[2,5]<-rcorrcens(Surv(time_PSA,PSA_status)~I(-1*as.numeric(Stage_groups)),data=clin)[1]
output_cindex[2,6]<-rcorrcens(Surv(time_PSA,PSA_status)~I(-1*as.numeric(Stage_groups)),data=clin)[6]

output_cindex[2,7]<-round(summary(coxph(Surv(time_clinical,clinical_status)~Stage_groups,data=clin))$conf.int[,1],digits=2)
output_cindex[2,8]<-paste( round(summary(coxph(Surv(time_clinical,clinical_status)~Stage_groups,data=clin))$conf.int[,3],digits=2),"-",
                           round(summary(coxph(Surv(time_clinical,clinical_status)~Stage_groups,data=clin))$conf.int[,4], digits=2),sep="")
output_cindex[2,9]<-coef(summary(coxph(Surv(time_clinical,clinical_status)~Stage_groups,data=clin)))[,5]
output_cindex[2,10]<-summary(coxph(Surv(time_clinical,clinical_status)~Stage_groups,data=clin))$concordance[1]
output_cindex[2,11]<-rcorrcens(Surv(time_clinical,clinical_status)~I(-1*as.numeric(Stage_groups)),data=clin)[1]
output_cindex[2,12]<-rcorrcens(Surv(time_clinical,clinical_status)~I(-1*as.numeric(Stage_groups)),data=clin)[6]

output_cindex[2,13]<-round(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~Stage_groups,data=clin))$conf.int[,1],digits=2)
output_cindex[2,14]<-paste( round(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~Stage_groups,data=clin))$conf.int[,3],digits=2),"-",
                            round(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~Stage_groups,data=clin))$conf.int[,4], digits=2),sep="")
output_cindex[2,15]<-coef(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~Stage_groups,data=clin)))[,5]
output_cindex[2,16]<-summary(coxph(Surv(Time_last_FU,Death_PCa_status)~Stage_groups,data=clin))$concordance[1]
output_cindex[2,17]<-rcorrcens(Surv(Time_last_FU,Death_PCa_status)~I(-1*as.numeric(Stage_groups)),data=clin)[1]
output_cindex[2,18]<-rcorrcens(Surv(Time_last_FU,Death_PCa_status)~I(-1*as.numeric(Stage_groups)),data=clin)[6]

### c. Dichotamised PSA level (<10 vs >10ng/mL)
output_cindex[3,1]<-round(summary(coxph(Surv(time_PSA,PSA_status)~PSA_factor,data=clin))$conf.int[,1],digits=2)
output_cindex[3,2]<-paste( round(summary(coxph(Surv(time_PSA,PSA_status)~PSA_factor,data=clin))$conf.int[,3],digits=2),"-",
                           round(summary(coxph(Surv(time_PSA,PSA_status)~PSA_factor,data=clin))$conf.int[,4], digits=2),sep="")
output_cindex[3,3]<-coef(summary(coxph(Surv(time_PSA,PSA_status)~PSA_factor,data=clin)))[,5]
output_cindex[3,4]<-summary(coxph(Surv(time_PSA,PSA_status)~PSA_factor,data=clin))$concordance[1]
output_cindex[3,5]<-rcorrcens(Surv(time_PSA,PSA_status)~I(-1*as.numeric(PSA_factor)),data=clin)[1]
output_cindex[3,6]<-rcorrcens(Surv(time_PSA,PSA_status)~I(-1*as.numeric(PSA_factor)),data=clin)[6]

output_cindex[3,7]<-round(summary(coxph(Surv(time_clinical,clinical_status)~PSA_factor,data=clin))$conf.int[,1],digits=2)
output_cindex[3,8]<-paste( round(summary(coxph(Surv(time_clinical,clinical_status)~PSA_factor,data=clin))$conf.int[,3],digits=2),"-",
                           round(summary(coxph(Surv(time_clinical,clinical_status)~PSA_factor,data=clin))$conf.int[,4], digits=2),sep="")
output_cindex[3,9]<-coef(summary(coxph(Surv(time_clinical,clinical_status)~PSA_factor,data=clin)))[,5]
output_cindex[3,10]<-summary(coxph(Surv(time_clinical,clinical_status)~PSA_factor,data=clin))$concordance[1]
output_cindex[3,11]<-rcorrcens(Surv(time_clinical,clinical_status)~I(-1*as.numeric(PSA_factor)),data=clin)[1]
output_cindex[3,12]<-rcorrcens(Surv(time_clinical,clinical_status)~I(-1*as.numeric(PSA_factor)),data=clin)[6]

output_cindex[3,13]<-round(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~PSA_factor,data=clin))$conf.int[,1],digits=2)
output_cindex[3,14]<-paste( round(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~PSA_factor,data=clin))$conf.int[,3],digits=2),"-",
                            round(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~PSA_factor,data=clin))$conf.int[,4], digits=2),sep="")
output_cindex[3,15]<-coef(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~PSA_factor,data=clin)))[,5]
output_cindex[3,16]<-summary(coxph(Surv(Time_last_FU,Death_PCa_status)~PSA_factor,data=clin))$concordance[1]
output_cindex[3,17]<-rcorrcens(Surv(Time_last_FU,Death_PCa_status)~I(-1*as.numeric(PSA_factor)),data=clin)[1]
output_cindex[3,18]<-rcorrcens(Surv(Time_last_FU,Death_PCa_status)~I(-1*as.numeric(PSA_factor)),data=clin)[6]

### d. Dichotamised Margin Status (Negative vs Positive)
output_cindex[4,1]<-round(summary(coxph(Surv(time_PSA,PSA_status)~margin_group,data=clin))$conf.int[,1],digits=2)
output_cindex[4,2]<-paste( round(summary(coxph(Surv(time_PSA,PSA_status)~margin_group,data=clin))$conf.int[,3],digits=2),"-",
                           round(summary(coxph(Surv(time_PSA,PSA_status)~margin_group,data=clin))$conf.int[,4], digits=2),sep="")
output_cindex[4,3]<-coef(summary(coxph(Surv(time_PSA,PSA_status)~margin_group,data=clin)))[,5]
output_cindex[4,4]<-summary(coxph(Surv(time_PSA,PSA_status)~margin_group,data=clin))$concordance[1]
output_cindex[4,5]<-rcorrcens(Surv(time_PSA,PSA_status)~I(-1*as.numeric(margin_group)),data=clin)[1]
output_cindex[4,6]<-rcorrcens(Surv(time_PSA,PSA_status)~I(-1*as.numeric(margin_group)),data=clin)[6]

output_cindex[4,7]<-round(summary(coxph(Surv(time_clinical,clinical_status)~margin_group,data=clin))$conf.int[,1],digits=2)
output_cindex[4,8]<-paste( round(summary(coxph(Surv(time_clinical,clinical_status)~margin_group,data=clin))$conf.int[,3],digits=2),"-",
                           round(summary(coxph(Surv(time_clinical,clinical_status)~margin_group,data=clin))$conf.int[,4], digits=2),sep="")
output_cindex[4,9]<-coef(summary(coxph(Surv(time_clinical,clinical_status)~margin_group,data=clin)))[,5]
output_cindex[4,10]<-summary(coxph(Surv(time_clinical,clinical_status)~margin_group,data=clin))$concordance[1]
output_cindex[4,11]<-rcorrcens(Surv(time_clinical,clinical_status)~I(-1*as.numeric(margin_group)),data=clin)[1]
output_cindex[4,12]<-rcorrcens(Surv(time_clinical,clinical_status)~I(-1*as.numeric(margin_group)),data=clin)[6]

output_cindex[4,13]<-round(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~margin_group,data=clin))$conf.int[,1],digits=2)
output_cindex[4,14]<-paste( round(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~margin_group,data=clin))$conf.int[,3],digits=2),"-",
                            round(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~margin_group,data=clin))$conf.int[,4], digits=2),sep="")
output_cindex[4,15]<-coef(summary(coxph(Surv(Time_last_FU,Death_PCa_status)~margin_group,data=clin)))[,5]
output_cindex[4,16]<-summary(coxph(Surv(Time_last_FU,Death_PCa_status)~margin_group,data=clin))$concordance[1]
output_cindex[4,17]<-rcorrcens(Surv(Time_last_FU,Death_PCa_status)~I(-1*as.numeric(margin_group)),data=clin)[1]
output_cindex[4,18]<-rcorrcens(Surv(Time_last_FU,Death_PCa_status)~I(-1*as.numeric(margin_group)),data=clin)[6]

#write out output_cindex for second part of Table S6 and Table 3c_ii
write.csv(output_cindex, file = "MBPS/results/Clinicopathological_univariate_cox_cindex.csv")


#2. Univariable methylation analysis


##Kaplan-Meier plots and log-rank p-vals for univariable methylation analysis (Figure 3b, Figure S10 and Table S6)

# Prepare methylation data to dichotamise each methylation variable using its 75% percentile (third quartile) as cut-off

#remove sample with no PSA data to match clinicopathological analysis
m<-match(clin$Sample_ID,colnames(avM_samp))
avM_samp<-avM_samp[,m]

avM_factor75<-avM_samp
for(i in 1:nrow(avM_samp)) {
  avM_factor75[i,]<-ifelse(avM_samp[i,]<=quantile(avM_samp[i,],0.75,na.rm=T),0,1)
}

tavM_factor75<-t(avM_factor75)
clin_75<-cbind(tavM_factor75,clin)

for(i in 1:18){
  
  ####BCR
  fit1<-survfit(Surv(time_PSAm,PSA_status)~clin_75[,i],data=clin_75)
  
  ####MR
  fit2<-survfit(Surv(time_clinicalm,clinical_status)~clin_75[,i],data=clin_75)
  
  ####PCSM
  fit3<-survfit(Surv(time_FUm,Death_PCa_status)~clin_75[,i],data=clin_75)
  
  plot_list=list()
  ####BCR
  plot_list[[1]]<-ggsurvplot(fit1,data=clin_75,palette=c("mediumblue","firebrick1"),
                             pval=TRUE,xlab="Time to BCR (months)",title=paste("Biochemical recurrence - ",namelist[i],sep=""),
                             break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                             risk.table.y.text.col=T,
                             legend.labs = 
                               c("< 75th percentile", ">= 75th percentile"))

  ####MR
  plot_list[[2]]<-ggsurvplot(fit2,data=clin_75,palette=c("mediumblue","firebrick1"),
                             pval=TRUE,xlab="Time to clinical relapse (months)",title=paste("Clinical relapse - ",namelist[i],sep=""),
                             break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                             risk.table.y.text.col=T,
                             legend.labs = 
                               c("< 75th percentile", ">= 75th percentile"))
  
  ####PCSM
  plot_list[[3]]<-ggsurvplot(fit3,data=clin_75,palette=c("mediumblue","firebrick1"),
                             pval=TRUE,xlab="Survival time (months)",title=paste("Prostate Cancer mortality - ",namelist[i],sep=""),
                             break.time.by=50, ggtheme=theme_classic(),risk.table=T,fontsize=3,#tables.theme=clean_theme(),
                             risk.table.y.text.col=T,
                             legend.labs = 
                               c("< 75th percentile", ">= 75th percentile"))
  
  
  pdf(paste("MBPS/plots/survivalanalysis/",namelist[i],"_survival_KM.pdf",sep=""),useDingbats = FALSE)
  for (j in 1:3){
    print(plot_list[[j]])
  }
  dev.off()
  

 }



##Cox analysis and c-index for univariable methylation analysis (Table S6)
output_cindexM<-matrix(data=NA,nrow=18,ncol=18)
rownames(output_cindexM) <- namelist
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

for(i in 1:length(namelist)){
	temp<-namelist[i]
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

#write out output_cindex for second part of Supp Table 6
write.csv(output_cindexM, file = "MBPS/results/Methylation_univariate_cox_cindex.csv")

##Repeat Cox analysis and c-index for univariable methylation analysis (not included in manuscript)

clin_cont<-cbind(t(avM_samp),clin) 

output_cindexMcont<-matrix(data=NA,nrow=18,ncol=18)
rownames(output_cindexMcont) <- namelist
colnames(output_cindexMcont) <- c("Biochemical relapse: HR", "Biochemical relapse: HR 95%CI", 
                             "Biochemical relapse: HR pval","Biochemical relapse: coxph C-index",
                             "Biochemical relapse: rcorrcens C-index", "Biochemical relapse: C-index pval",
                             "Clinical relapse: HR", "Clinical relapse: HR 95%CI", 
                             "Clinical relapse: HR pval","Clinical relapse: coxph C-index",
                             "Clinical relapse: rcorrcens C-index", "Clinical relapse: C-index pval",
                             "Prostate cancer mortality: HR", "Prostate cancer mortality: HR 95%CI", 
                             "Prostate cancer mortality: HR pval","Prostate cancer mortality: coxph C-index",
                             "Prostate cancer mortality: rcorrcens C-index", "Prostate cancer mortality: C-index pval"
)

for(i in 1:length(namelist)){
	temp<-namelist[i]
	surv<-Surv(clin_75$time_PSA,clin_75$PSA_status)
	output_cindexMcont[i,1]<-round(summary(coxph(surv~clin_cont[,temp]))$conf.int[,1],digits=2)
	output_cindexMcont[i,2]<-paste( round(summary(coxph(surv~clin_cont[,temp]))$conf.int[,3],digits=2),"-",
                           round(summary(coxph(surv~clin_cont[,temp]))$conf.int[,4], digits=2),sep="")
	output_cindexMcont[i,3]<-coef(summary(coxph(surv~clin_cont[,temp])))[,5]
	output_cindexMcont[i,4]<-summary(coxph(surv~clin_cont[,temp]))$concordance[1]
	output_cindexMcont[i,5]<-rcorrcens(surv~I(-1*clin_cont[,temp]))[1]
	output_cindexMcont[i,6]<-rcorrcens(surv~I(-1*clin_cont[,temp]))[6]

	surv2<-Surv(clin_75$time_clinical,clin_75$clinical_status)
	output_cindexMcont[i,7]<-round(summary(coxph(surv2~clin_cont[,temp]))$conf.int[,1],digits=2)
	output_cindexMcont[i,8]<-paste( round(summary(coxph(surv2~clin_cont[,temp]))$conf.int[,3],digits=2),"-",
                           round(summary(coxph(surv2~clin_cont[,temp]))$conf.int[,4], digits=2),sep="")
	output_cindexMcont[i,9]<-coef(summary(coxph(surv2~clin_cont[,temp])))[,5]
	output_cindexMcont[i,10]<-summary(coxph(surv2~clin_cont[,temp]))$concordance[1]
	output_cindexMcont[i,11]<-rcorrcens(surv2~I(-1*clin_cont[,temp]))[1]
	output_cindexMcont[i,12]<-rcorrcens(surv2~I(-1*clin_cont[,temp]))[6]

	surv3<-Surv(clin_75$Time_last_FU,clin_75$Death_PCa_status)
	output_cindexMcont[i,13]<-round(summary(coxph(surv3~clin_cont[,temp]))$conf.int[,1],digits=2)
	output_cindexMcont[i,14]<-paste( round(summary(coxph(surv3~clin_cont[,temp]))$conf.int[,3],digits=2),"-",
                           round(summary(coxph(surv3~clin_cont[,temp]))$conf.int[,4], digits=2),sep="")
	output_cindexMcont[i,15]<-coef(summary(coxph(surv3~clin_cont[,temp])))[,5]
	output_cindexMcont[i,16]<-summary(coxph(surv3~clin_cont[,temp]))$concordance[1]
	output_cindexMcont[i,17]<-rcorrcens(surv3~I(-1*clin_cont[,temp]))[1]
	output_cindexMcont[i,18]<-rcorrcens(surv3~I(-1*clin_cont[,temp]))[6]
}
write.csv(output_cindexMcont, file = "MBPS/results/MethylationContinuous_univariate_cox_cindex.csv")



#Multivariable clinopathological and methylation analysis

##Using SES for model selection and then extracting results for manuscript
 varInt<-c("AC074091.13","BHLHB9","CACNA2D4","CALM1","CCDC67","CD34","CD01","EFCAB4B","EPHB3","LECT1","PAH","PARP6_CELF6","PRDM8","RP11-927P21.4","TBX1","Y_RNA","ZNF655","ISUP_groups","Stage_groups","margin_group","PSA_factor")
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

##Collate results above in table for Table 3
output_cindex_SES_75split<-matrix(data=NA,nrow=3,ncol=6)
rownames(output_cindex_SES_75split) <- c("CACNA2D4 (Biochemical recurrence)+ Grade + PSA",
                                               "Grade + Margin Status (Clinical relapse)",
                                               "CACNA2D4 (Prostate Cancer Mortality) + Grade")
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

#write out results of SES model
write.csv(output_cindex_SES_75split, file = "MBPS/results/SES_cindex_75split.csv")


#Extract clinicopathological variables only for Table 3a_ii

PSA.isup.psa<-coxph(surv~PSA_factor+ISUP_groups,data=clin_75)
lp.PSA.isup.psa<-predict(PSA.isup.psa,type="lp")

#. Gleason grade + PSA levels
output_cindex_multi<-vector(length=6)

output_cindex_multi[1]<-paste(round(summary(PSA.isup.psa)$conf.int[1,1],digits=3),";",round(summary(PSA.isup.psa)$conf.int[2,1],digits=3),sep="")
output_cindex_multi[2]<-paste(round(summary(PSA.isup.psa)$conf.int[1,3],digits=2),"-",
                                round(summary(PSA.isup.psa)$conf.int[1,4], digits=2),";",
                                round(summary(PSA.isup.psa)$conf.int[2,3],digits=2),"-",
                                round(summary(PSA.isup.psa)$conf.int[2,4], digits=2),sep="")
output_cindex_multi[3]<-paste(round(coef(summary(PSA.isup.psa))[1,5],digits=5),";",round(coef(summary(PSA.isup.psa))[2,5],digits=5),sep="")
output_cindex_multi[4]<-summary(PSA.isup.psa)$concordance[1]
output_cindex_multi[5]<-rcorrcens(surv~I(-1*lp.PSA.isup.psa))[1]
output_cindex_multi[6]<-rcorrcens(surv~I(-1*lp.PSA.isup.psa))[6]

names(output_cindex_multi) <- c("PSA relapse: HR", "PSA relapse: HR 95%CI", 
                                   "PSA relapse: HR pval","PSA relapse: coxph C-index","PSA relapse: rcorrcens C-index","PSA relapse: C-index pval")

write.csv(output_cindex_multi,file="MBPS/results/CompareClinPathOnly_BCR.csv")


# Compare models with and without CACNA2D4 with timeROC package - Table S7

#set time points for assessing ROC
times<-c(365,5*365,10*365,15*365)
times2<-c(1,5,10,15)

#changeISUP grade group variable into 0 and 1
clinnewer<-vector()
clinnewer[which(clin_75$ISUP_groups=="Group 1 (3+4)")]<-0
clinnewer[which(clin_75$ISUP_groups=="Group 2 (4+3, 8, 9)")]<-1

#fit multivariable cox models, with and without CACNA2D4, for BCR
fit0 <- coxph( surv~ clin_75$CACNA2D4 + clin_75$PSA_factor + clinnewer, na.action=na.omit )
eta <- fit0$linear.predictor

fit1 <- coxph( surv~ clin_75$PSA_factor + clinnewer, na.action=na.omit )
eta1 <- fit1$linear.predictor

#calculate ROC curves for cox models at each specific time point
ROCwMeth<-ROCnoMeth<-list()
for(i in 1:length(times)){
ROCwMeth[[i]]=risksetROC(Stime=clin_75$time_PSA, status=clin_75$PSA_status,
                    marker=eta, predict.time=times[i], method="Cox",
                    main="ROC Curve", lty=2, col="red") 
                   
ROCnoMeth[[i]]=risksetROC(Stime=clin_75$time_PSA, status=clin_75$PSA_status,
                    marker=eta1, predict.time=times[i], method="Cox",
                    main="ROC Curve", lty=2, col="red")                     
                    
}

#extract AUC
AUC_with<-unlist(lapply(ROCwMeth,function(x){x$AUC}))
AUC_without<-unlist(lapply(ROCnoMeth,function(x){x$AUC}))

#plot ROC curves on same plot
pdf("MBPS/results/CompareCoxBCR_timeROC.pdf",height=9,width=9)
par(mfrow=c(2,2))

for(i in 1:length(times)){

plot(ROCwMeth[[i]]$FP, ROCwMeth[[i]]$TP, type = "l", xlab = "FP", ylab = "TP",main=paste(times2[i],"years"))
abline(c(0, 0), c(1, 1))

lines(ROCnoMeth[[i]]$FP, ROCnoMeth[[i]]$TP, type = "l", xlab = "FP", ylab = "TP",lty=3)

legend("bottomright",lty=c(1,3),legend=c(round(ROCwMeth[[i]]$AUC,2),round(ROCnoMeth[[i]]$AUC,2)))
}
dev.off()


#fit multivariable cox models, with and without CACNA2D4, for PCSM
fit0D <- coxph( surv3~ clin_75$CACNA2D4 + clinnewer, na.action=na.omit )
etaD <- fit0D$linear.predictor

fit1D <- coxph( surv3~ clinnewer, na.action=na.omit )
eta1D <- fit1D$linear.predictor

#calculate ROC curves for cox models at each specific time point
ROCwMethD<-ROCnoMethD<-list()
for(i in 1:length(times)){
ROCwMethD[[i]]=risksetROC(Stime=clin_75$Time_last_FU, status=clin_75$Death_PCa_status,
                    marker=etaD, predict.time=times[i], method="Cox",
                    main="ROC Curve", lty=2, col="red") 
                   
ROCnoMethD[[i]]=risksetROC(Stime=clin_75$Time_last_FU, status=clin_75$Death_PCa_status,
                    marker=eta1D, predict.time=times[i], method="Cox",
                    main="ROC Curve", lty=2, col="red")                     
                    
}

#extract AUC
AUC_withD<-unlist(lapply(ROCwMethD,function(x){x$AUC}))
AUC_withoutD<-unlist(lapply(ROCnoMethD,function(x){x$AUC}))

pdf("MBPS/results/CompareCoxPCSM_timeROC.pdf",height=9,width=9)
par(mfrow=c(2,2))

for(i in 1:length(times)){

plot(ROCwMethD[[i]]$FP, ROCwMethD[[i]]$TP, type = "l", xlab = "False Positive", ylab = "True Positive",main=paste(times2[i],"years"))
abline(c(0, 0), c(1, 1))

lines(ROCnoMethD[[i]]$FP, ROCnoMethD[[i]]$TP, type = "l", xlab = "False Positive", ylab = "True Positive",lty=3)

legend("bottomright",lty=c(1,3),legend=c(round(ROCwMethD[[i]]$AUC,2),round(ROCnoMethD[[i]]$AUC,2)))
}
dev.off()

#fit multivariable cox models for MR
fit0M <- coxph( surv2~ clin_75$CACNA2D4+clinnewer + clin_75$margin_group, na.action=na.omit )
etaM <- fit0M$linear.predictor
ROCnoMethM<-list()
for(i in 1:length(times)){                   
ROCnoMethM[[i]]=risksetROC(Stime=clin_75$time_clinical, status=clin_75$clinical_status,
                    marker=etaM, predict.time=times[i], method="Cox",
                    main="ROC Curve", lty=2, col="red")                     
                    
}

#extract AUC
AUC_withoutM<-unlist(lapply(ROCnoMethM,function(x){x$AUC}))

pdf("MBPS/results/CompareCoxMR_timeROC.pdf",height=9,width=9)
par(mfrow=c(2,2))

for(i in 1:length(times)){

plot(ROCnoMethM[[i]]$FP, ROCnoMethM[[i]]$TP, type = "l", xlab = "False Positive", ylab = "True Positive",main=paste(times2[i],"years"))
abline(c(0, 0), c(1, 1))


legend("bottomright",lty=c(1,3),legend=round(ROCnoMethM[[i]]$AUC,2))
}
dev.off()

alltimes<-rbind(AUC_with,AUC_without,AUC_withoutM,AUC_withD,AUC_withoutD)
colnames(alltimes)<-times

#export as Table S7
write.csv(alltimes,"MBPS/results/CompareCox_timeROC.csv")


# Generate forest plots for univariable Cox analyses -- separate by BCR, MR and PCSM -- Figure 3a and Figure S9

#read in univariable Cox results for multiplex methylation data
meth_res<-read.csv("MBPS/results/Methylation_univariate_cox_cindex.csv",row.names=1)

#read in univariable Cox results for clinicopathological data
clin_res<-read.csv("MBPS/results/Clinicopathological_univariate_cox_cindex.csv",row.names=1)

#Generate full Forest plots for Figure S9

# 1. BCR: make table using HRs and 95% CI

BCR_data <-
  structure(list(
    HR  = meth_res$Biochemical.relapse..HR,
    lower = as.numeric(unlist(lapply(meth_res$Biochemical.relapse..HR.95.CI,function(x){strsplit(as.character(x),"-")[[1]][1]}))),
    upper = as.numeric(unlist(lapply(meth_res$Biochemical.relapse..HR.95.CI,function(x){strsplit(as.character(x),"-")[[1]][2]})))),
    .Names = c("HR", "lower", "upper"),
    row.names = namelist,
    class = "data.frame")

BCR_tabletext <- cbind(namelist,round(as.numeric(meth_res$Biochemical.relapse..HR.pval),3))
colnames(BCR_tabletext)<-NULL

#manually reorder amplicons to order in genome
reord<-c(6,1,10,14,17, 8, 18, 5, 3, 9, 12, 11, 4, 13,15, 16,2)

BCR_data2<-BCR_data[reord,]
BCR_tabletext2<-BCR_tabletext[reord,]

#append clinical data
BCR_dataC <-
  structure(list(
    HR  = clin_res$Biochemical.relapse..HR,
    lower = as.numeric(unlist(lapply(clin_res$Biochemical.relapse..HR.95.CI,function(x){strsplit(as.character(x),"-")[[1]][1]}))),
    upper = as.numeric(unlist(lapply(clin_res$Biochemical.relapse..HR.95.CI,function(x){strsplit(as.character(x),"-")[[1]][2]})))),
    .Names = c("HR", "lower", "upper"),
    row.names = rownames(clin_res),
    class = "data.frame")

BCR_tabletextC <- cbind(rownames(clin_res),round(as.numeric(clin_res$Biochemical.relapse..HR.pval),3))
colnames(BCR_tabletextC)<-NULL

BCR_dataC2<-BCR_dataC[c(1,3,2,4),]
BCR_tabletextC2<-BCR_tabletextC[c(1,3,2,4),]

BCR_dat<-rbind(BCR_data2,BCR_dataC2)
BCR_tt<-rbind(BCR_tabletext2,BCR_tabletextC2)

pdf("MBPS/plots/survivalanalysis/forestplot_BCR_methylation.pdf",width=8,useDingbats = FALSE)
forestplot(BCR_tt,
           BCR_dat,
           graph.pos = 2,
           boxsize = 0.1,
           zero = 1,
           lwd.zero = 2,
           xlog=FALSE,
           xticks=c(0,1,2,3,4),
           col=fpColors(box="blue",line="darkblue",summary="blue"),
           xlab = "Hazard Ratio (95% CI)",
           title = "Biochemical recurrence")
dev.off()

# 2. MR: make table using HRs and 95% CI

MR_data <-
  structure(list(
    HR  = meth_res$Clinical.relapse..HR,
    lower = as.numeric(unlist(lapply(meth_res$Clinical.relapse..HR.95.CI,function(x){strsplit(as.character(x),"-")[[1]][1]}))),
    upper = as.numeric(unlist(lapply(meth_res$Clinical.relapse..HR.95.CI,function(x){strsplit(as.character(x),"-")[[1]][2]})))),
    .Names = c("HR", "lower", "upper"),
    row.names = namelist,
    class = "data.frame")

MR_tabletext <- cbind(namelist,round(as.numeric(meth_res$Clinical.relapse..HR.pval),3))
colnames(MR_tabletext)<-NULL

#manually reorder amplicons to order in genome
MR_data2<-MR_data[reord,]
MR_tabletext2<-MR_tabletext[reord,]

#append clinical data
MR_dataC <-
  structure(list(
    HR  = clin_res$Clinical.relapse..HR,
    lower = as.numeric(unlist(lapply(clin_res$Clinical.relapse..HR.95.CI,function(x){strsplit(as.character(x),"-")[[1]][1]}))),
    upper = as.numeric(unlist(lapply(clin_res$Clinical.relapse..HR.95.CI,function(x){strsplit(as.character(x),"-")[[1]][2]})))),
    .Names = c("HR", "lower", "upper"),
    row.names = rownames(clin_res),
    class = "data.frame")

MR_tabletextC <- cbind(rownames(clin_res),round(as.numeric(clin_res$Clinical.relapse..HR.pval),3))
colnames(MR_tabletextC)<-NULL

MR_dataC2<-MR_dataC[c(1,3,2,4),]
MR_tabletextC2<-MR_tabletextC[c(1,3,2,4),]

MR_dat<-rbind(MR_data2,MR_dataC2)
MR_tt<-rbind(MR_tabletext2,MR_tabletextC2)

pdf("MBPS/plots/survivalanalysis/forestplot_MR_methylation.pdf",width=8,useDingbats = FALSE)
forestplot(MR_tt,
           MR_dat,
           graph.pos = 2,
           boxsize = 0.1,
           zero = 1,
           lwd.zero = 2,
           xlog=FALSE,
           xticks=0:15,
           col=fpColors(box="blue",line="darkblue",summary="blue"),
           xlab = "Hazard Ratio (95% CI)",
           title = "Metastatic relapse")
dev.off()

# 3. PCSM: make table using HRs and 95% CI

PC_data <-
  structure(list(
    HR  = meth_res$Prostate.cancer.mortality..HR,
    lower = as.numeric(unlist(lapply(meth_res$Prostate.cancer.mortality..HR.95.CI,function(x){strsplit(as.character(x),"-")[[1]][1]}))),
    upper = as.numeric(unlist(lapply(meth_res$Prostate.cancer.mortality..HR.95.CI,function(x){strsplit(as.character(x),"-")[[1]][2]})))),
    .Names = c("HR", "lower", "upper"),
    row.names = namelist,
    class = "data.frame")

PC_tabletext <- cbind(namelist,round(as.numeric(meth_res$Prostate.cancer.mortality..HR.pval),3))
colnames(PC_tabletext)<-NULL

#manually reorder amplicons to order in genome
PC_data2<-PC_data[reord,]
PC_tabletext2<-PC_tabletext[reord,]

#append clinical data
PC_dataC <-
  structure(list(
    HR  = clin_res$Prostate.cancer.mortality..HR,
    lower = as.numeric(unlist(lapply(clin_res$Prostate.cancer.mortality..HR.95.CI,function(x){strsplit(as.character(x),"-")[[1]][1]}))),
    upper = as.numeric(unlist(lapply(clin_res$Prostate.cancer.mortality..HR.95.CI,function(x){strsplit(as.character(x),"-")[[1]][2]})))),
    .Names = c("HR", "lower", "upper"),
    row.names = rownames(clin_res),
    class = "data.frame")

PC_tabletextC <- cbind(rownames(clin_res),round(as.numeric(clin_res$Prostate.cancer.mortality..HR.pval),3))
colnames(PC_tabletextC)<-NULL

PC_dataC2<-PC_dataC[c(1,3,2,4),]
PC_tabletextC2<-PC_tabletextC[c(1,3,2,4),]

PC_dat<-rbind(PC_data2,PC_dataC2)
PC_tt<-rbind(PC_tabletext2,PC_tabletextC2)

pdf("MBPS/plots/survivalanalysis/forestplot_PCSM_methylation.pdf",width=8,useDingbats = FALSE)
forestplot(PC_tt,
PC_dat,
           graph.pos = 2,
           boxsize = 0.1,
           zero = 1,
           lwd.zero = 2,
           xlog=FALSE,
           xticks=0:16,
           col=fpColors(box="blue",line="darkblue",summary="blue"),
           xlab = "Hazard Ratio (95% CI)",
           title = "Prostate cancer specific mortality")
dev.off()

#Export versions of forest plots with just significant HR - Figure 3a

# 1. BCR
sigBCR<-which(BCR_tt[,2]<0.05)
BCR_tt_sig<-BCR_tt[sigBCR,]
BCR_dat_sig<-BCR_dat[sigBCR,]
pdf("MBPS/plots/survivalanalysis/forestplot_BCR_methylation_sig.pdf",width=8,useDingbats = FALSE)
forestplot(BCR_tt_sig,
            BCR_dat_sig,
           graph.pos = 2,
           boxsize = 0.1,
           zero = 1,
           lwd.zero = 2,
           xlog=FALSE,
           xticks=0:16,
           col=fpColors(box="blue",line="darkblue",summary="blue"),
           xlab = "Hazard Ratio (95% CI)",
           title = "Biochemical Recurrence")
dev.off()

# 2. MR
sigMR<-which(MR_tt[,2]<0.05)
MR_tt_sig<-MR_tt[sigMR,]
MR_dat_sig<-MR_dat[sigMR,]
pdf("MBPS/plots/survivalanalysis/forestplot_MR_methylation_sig.pdf",width=8,useDingbats = FALSE)
forestplot(MR_tt_sig,
MR_dat_sig,
           graph.pos = 2,
           boxsize = 0.1,
           zero = 1,
           lwd.zero = 2,
           xlog=FALSE,
           xticks=0:16,
           col=fpColors(box="blue",line="darkblue",summary="blue"),
           xlab = "Hazard Ratio (95% CI)",
           title = "Metastatic Relapse")
dev.off()

# 3. PCSM
sigPC<-which(PC_tt[,2]<0.05)
PC_tt_sig<-PC_tt[sigPC,]
PC_dat_sig<-PC_dat[sigPC,]
pdf("MBPS/plots/survivalanalysis/forestplot_PCSM_methylation_sig.pdf",width=8,useDingbats = FALSE)
forestplot(PC_tt_sig,
PC_dat_sig,
           graph.pos = 2,
           boxsize = 0.1,
           zero = 1,
           lwd.zero = 2,
           xlog=FALSE,
           xticks=0:16,
           col=fpColors(box="blue",line="darkblue",summary="blue"),
           xlab = "Hazard Ratio (95% CI)",
           title = "Prostate cancer specific mortality")
dev.off()
