###Script for selecting TCGA PRAD samples with 450K methylation for analysis and finding raw idat barcodes ready for loading in data

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#make functions for extracting sample ID information from clinical and sample data matrices

#Patient identifier
patient_identifier<-function(x){
  paste(strsplit(as.character(x),"-")[[1]][1:3],collapse="-")
  
}

#Sample identifier
sample_identifier<-function(x){
  substr(as.character(x),1,15)
  
}

#Sample identifier including portion
sample_identifier_modif<-function(x){
  substr(as.character(x),1,16)
  
}

#Sample identifier including portion 2
site_identifier<-function(x){
  substr(as.character(x),6,7)
  
}

#read in TCGA sample, annotation and clinical data, originally downloaded from TCGA using the EGC.tools R package
(https://github.com/uscepigenomecenter/EGC.tools)

#read in TCGA PRAD Sample and Data Relationship Form (SDRF), which describes the relationship between samples, arrays, data, and other objects used or produced in the TCGA PRAD study

sdrf<-read.csv("metadata/jhu-usc.edu_PRAD.HumanMethylation450.1.13.0.sdrf.txt",sep="\t",header=TRUE)

#read in biospecimen data for TCGA PRAD samples
biospecimen<-read.csv("metadata/nationwidechildrens.org_biospecimen_sample_prad.txt",sep="\t",header=TRUE)

#remove rows without data
biospecimen2<-biospecimen[2:nrow(biospecimen),]

#read in clinical data for TCGA PRAD samples
clinical<-read.csv("metadata/nationwidechildrens.org_clinical_cqcf_prad.txt",sep="\t",header=TRUE)

#trim off NA rows without data
clinical2<-clinical[3:nrow(clinical),]

#Subset SDRF file for Solid Tissue Normal samples with 450K data only
sdrf_sampleID<-sapply(sdrf$Comment..TCGA.Barcode.,sample_identifier)
biospecimen_sampleID<-sapply(biospecimen2$bcr_sample_barcode,sample_identifier)
biospecimen_patientID<-sapply(biospecimen2$bcr_sample_barcode,patient_identifier)

STN<-which(biospecimen2$sample_type=="Solid Tissue Normal")
biospecimen_sampleID_STN<-biospecimen_sampleID[STN]
biospecimen_patientID_STN<-biospecimen_patientID[STN]

matchup<-match(biospecimen_sampleID_STN,sdrf_sampleID)

sdrf_STN<-sdrf[matchup,]
sdrf_STN2<-sdrf_STN[-which(is.na(sdrf_STN$Protocol.REF)==TRUE),]
#sdrf_STN2: array details for all solid tissue normal with 450K data

#extract patientIDs of those with normal tissue and 450K data
sdrf_STN_patientID<-sapply(sdrf_STN2$Comment..TCGA.Barcode,patient_identifier)

#Subset clinical tables for patients with normal tissue and 450K data
clinical_STN<-clinical[match(sdrf_STN_patientID,clinical$bcr_patient_barcode),]

#Subset SDRF file for Primary Tumour samples with 450K data only
PT<-which(biospecimen2$sample_type=="Primary Tumor")
biospecimen_sampleID_PT<-biospecimen_sampleID[PT]
biospecimen_patientID_PT<-biospecimen_patientID[PT]
matchup2<-match(biospecimen_sampleID_PT,sdrf_sampleID)
sdrf_PT<-sdrf[na.omit(matchup2),]

#extract patient IDs of those with primary tumour tissue
sdrf_PT_patientID<-sapply(sdrf_PT$Comment..TCGA.Barcode.,patient_identifier)

#Filter SDRF and clinical tables to include only samples from caucasian patients with Prostate Adenocarcinoma Acinar Type, as this best reflects the ethnicity and cancer type of the Discovery Cohort, thus reducing chances of confounding of epigenetic effects due to extraneous factors in comparisons between the datasets
clinicalPT<-clinical2[match(sdrf_PT_patientID,clinical2$bcr_patient_barcode),]

keep<-which(clinicalPT$race=="WHITE" & clinicalPT$histologic_diagnosis=="Prostate Adenocarcinoma Acinar Type" )
clinical_PTk<-clinicalPT[keep,]
sdrf_PTk<-sdrf_PT[keep,]

#Filter SDRF and clinical tables to remove duplicated samples (i.e. only include each sample once)
dups<-which(duplicated(sdrf_PTk$Comment..TCGA.Barcode.))
sdrf_PTk2<-sdrf_PTk[-dups,]
clinical_PTk2<-clinical_PTk[-dups,]
identical(as.character(sdrf_PTk2_patientID),as.character(clinical_PTk2$bcr_patient_barcode)) #TRUE

#extract patient IDs of those with primary tumour tissue remaining after filtering
sdrf_PTk2_patientID<-sapply(sdrf_PTk2$Comment..TCGA.Barcode.,patient_identifier)

#create file that contains TCGA PRAD sample name, 450K array barcode and basename, ready for reading in idat files

targets<-data.frame(Sample_Name=sdrf_PTk2$Comment..TCGA.Barcode.,
Array.Data.File=gsub("_Grn.idat","",sdrf_PTk2$Array.Data.File),
Comment..TCGA.Archive.Name.=sdrf_PTk2$Comment..TCGA.Archive.Name.,
Basename=paste("WGBS/data/TCGA_Level1/",gsub("_Grn.idat","",sdrf_PTk2$Array.Data.File),sep=""))

targets2<-data.frame(Sample_Name=sdrf_STN2$Comment..TCGA.Barcode.,
Array.Data.File=gsub("_Grn.idat","",sdrf_STN2$Array.Data.File),
Comment..TCGA.Archive.Name.=sdrf_STN2$Comment..TCGA.Archive.Name.,
Basename=paste("WGBS/data/TCGA_Level1/",gsub("_Grn.idat","",sdrf_STN2$Array.Data.File),sep=""))

targetsboth2<-rbind(targets2[1:50,],targets)

#read additional biospecimen data for TCGA PRAD normal samples using file originally downloaded from TCGA using the EGC.tools R package
(https://github.com/uscepigenomecenter/EGC.tools)
biospecimen_proximity<-read.csv("metadata/biospecimen_normal_control_prad.txt",sep="\t",header=TRUE)

#extract data for proximity of normal sample biopsy relative to tumour for normal samples with 450K data
biospecimen_proximity_STN<-biospecimen_proximity[match(sdrf_STN_patientID,biospecimen_proximity$bcr_patient_barcode),]

#Filter normal samples to only include samples from caucasian patients, as this best reflects the ethnicity and cancer type of the Discovery Cohort. And also filter the tables to only include normal samples that are more than 2cm from the tumour to minimise the chance that they may be showing an epigenetic field effect (and thus be more similar to the tumour than to normal tissue)
removesamp<-as.character(clinical_STN[which(clinical_STN$race!="WHITE"),2])

removesamp2<-as.character(biospecimen_proximity_STN[grep("Adjacent",biospecimen_proximity_STN$normal_tissue_proximity),1])

removesamp
#[1] "TCGA-EJ-7782" "TCGA-EJ-7789" "TCGA-HC-7737"

removesamp2
#[1] "TCGA-EJ-7782" "TCGA-EJ-7789" "TCGA-EJ-7792" "TCGA-EJ-7794"

rs<-unique(c(removesamp,removesamp2))

sample_identifier2<-function(x){
  substr(as.character(x),1,12)
  
}

targetsboth2_SN<-sample_identifier2(targetsboth2$Sample_Name)
mSN<-match(rs,targetsboth2_SN)

targetsboth3<-targetsboth2[-mSN,]

#write out csv file containing final sample IDs and idat barcodes of PRAD TCGA 450K samples used in analysis for this paper
write.csv(targetsboth3, "WGBS/data/targetsboth.csv",header=T)



