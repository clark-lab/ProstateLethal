# Scripts to perform analysis described in Pidsley et al [Journal TBC]

**01\_01\_DownloadCpGislandData.R** - used to download genomic coordinates of CpG island data for use in analysis

**01\_02\_GlobalCharacterisation.R** - contains details of WGBS data processing and plots WGBS PCA plot and methylation at CpG islands and repetitive elements (Figures 1a-b & S2)

**02\_01\_DMRcalling_plot.R** - calls DMRs between lethal and non-lethal patient groups from WGBS data and generates associated plots (Figures 1d & S1) and data for tables (Tables 2 & S2)

**02\_02\_DMRmean_heatmap.R** - extracts WGBS data to calculate mean methylation for lethal, non-lethal and normal samples for each of the 1420 DMRs and generates heatmap (Figure 1c) and data for tables (Tables S2 & S4)

**02\_03\_DownloadTranscriptData.R** - used to download hg19 transcript data for use in analysis

**02\_04\_AnnotateDMRs.R** - annotates DMRs with CpG and genic context and PrEC/LNCaP ChromHMM state. Generates data for tables and plots (for Tables 2, S2 & S4 and Figure 1e)

**03\_01\_PlotDMRcontext.R** - summarises overlap of DMRs with CpG islands and genic regions and makes plots (Figure 1e)

**03\_02\_01\_MakeDMRBedfiles_forGAT.R** - makes bed files of DMR genomic coordinates for use in GAT analysis

**03\_02\_02\_MakeWGBSBackgroundBedfile_forGAT.R** - makes bed file of WGBS CpG site genomic coordinates for use as background dataset in GAT analysis

**03\_02\_03\_RunGAT_DMRChromHMM.py** - used to run GAT software to find overlaps between DMRs and ChromHMM 

**03\_02\_04\_PlotGATresults.R** - plots overlap between DMRs and LNCaP and PrEC ChromHMM segmentation data (Figures 1f & S3)

**03\_03\_DMRMeth_LNCaPPrEC.R** - extracts methylation at all DMRs from LNCaP and PrEC WGBS data (Figure 1c and Table S4)

**03\_04\_GeneOntology.R** - performs gene ontology analysis on DMR associated genes (Figure 1g and Table S3)

**04\_01\_01\_CreateBedFilesDMRs.R** - creates bedfiles for WGBS data for lethal, non-lethal and normal groups for IGV tracks (Figure S1)

**04\_01\_02\_PlotWGBSselectedDMRs.R** - makes boxplots of WGBS methylation at 18 DMRs selected for validation (Figure S4a)

**04\_02\_ExtractingBlood_DMRmeth.R** - extracts DNA methylation at selected DMRs in blood using public data (for Table S4)

**04\_03\_01\_SelectingSamplesTCGA_HM450.R** - selects TCGA PRAD samples with 450K methylation data for analysis and extracts raw idat barcodes ready for loading in data

**04\_03\_02\_ProcessingTCGA_HM450.R** - reads in TCGA PRAD HM450 idat files, process data and performs basic QC

**04\_03\_03\_ExtractingTCGA_DMRmeth.R** - removes noisy probes from TCGA PRAD HM450 data and extracts methylation data at DMRs (Figure S4b & S1 and Table S4).  Compares DMR methylation and expression of nearest protein-coding gene (Figures S5 & S11b)

**04\_04\_PlotWGBS_LNCaPPrEC_selectedDMRs.R** - creates boxplot of LNCaP and PrEC WGBS data at 18 selected DMRs (Figure S4c)

**05\_01\_ProcessMBPSdata.sh** - processes MBPS data through MethPanel 

**05\_02\_PrepareMBPSdata.R** - prepares MBPS data from MethPanel ready for visualisation and survival analysis

**05\_03\_ExtractMBPS_QCdata.R** - extracts QC data for MBPS data (Table S5)

**05\_04\_VisualiseMBPS_QCdata.R** - creates plots of MBPS QC data (Figures S6-S8)

**06\_01\_MBPS_SurvivalAnalysis.R** - performs survival analysis of MBPS data and generates plots (Figures 2 & 3 and S9 & S10) and data (Tables S6 & S7)

**06\_02\_01\_MBPS_Downsample.sh** - downsamples raw MBPS data from bam files to  1K, 10K, 50K and 100K 

**06\_02\_02\_MBPSDownsampled_Processing.R** - prepares downsampled CACNA2D4 MBPS data for survival analysis

**06\_02\_03\_MBPSDownsampled_SurvivalAnalysis.R** - performs log rank, univariable Cox and multivariable Cox analysis of CACNA2D4 with downsampled MBPS data (Tables S8 & S9)
