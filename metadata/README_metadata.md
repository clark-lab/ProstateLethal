# Metadata 

##### **Input files required to run the scripts in this repository:** 

**AnnotateRegions_updown_update.R** - R code for annotating DMRs, adapted from a function in the aaRon package (https://github.com/astatham/aaRon/)

**GH_Interactions_Double_Elite_hg38.txt** - Genehancer Interactions Double Elite regions (downloaded from UCSC Table Browser)

**LNCaP_15_segments.bed** - bed file of ChromHMM segmentation data for LNCaP cells, generated as part of Giles 2019 Epigenetics and Chromatin (doi.org/10.1186/s13072-019-0258-9)  

**NoisyProbes.csv** - csv file containing details of "noisy probes" that should be excluded from the HM450 analysis (downloaded as Additional File 2 from http://www.biomedcentral.com/1471-2164/15/51)

**PrEC_15_segments.bed** - bed file of ChromHMM segmentation data for PrEC cells, generated as part of Giles 2019 Epigenetics and Chromatin (doi.org/10.1186/s13072-019-0258-9)

**SelectedRegions.bed** - bed file of 18 DMRs selected for validation

**WGBS_samplegroups.csv** - csv file of samples in Discovery Cohort with group information

**hg38ToHg19.over.chain.gz** - chain file for converting genomic coordinates from hg38 to hg19  (downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz)

**jhu-usc.edu_PRAD.HumanMethylation450.1.13.0.sdrf.txt** - TCGA PRAD Sample and Data Relationship Form (SDRF), which describes the relationship between samples, arrays, data, and other objects used or produced in the TCGA PRAD study

**nationwidechildrens.org_biospecimen_normal_control_prad.txt** - TCGA PRAD table giving information about the normal biospecimens assayed

**nationwidechildrens.org_biospecimen_sample_prad.txt** - TCGA PRAD table giving information about the biospecimens assayed

**nationwidechildrens.org_clinical_cqcf_prad.txt** - TCGA PRAD table giving clinical information about the patients corresponding to the samples assayed

All TCGA files can be downloaded using the GitHub ECG.tools R package
(https://github.com/uscepigenomecenter/EGC.tools)
