#Script calls DMRs between lethal and non-lethal patient groups from WGBS data
#Generates data for Table 2 and Table S2
#Generates plots for Figure 1d and Figure S1

#note, this script was run on an HPC rather than a personal computer as the BigTable data is large

#load libraries
library(data.table)
library(bsseq)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(DSS)
library(DMRcate)
library(Gviz)

#set user defined working directory and create new folders to output results
BASEDIR<-*user_defined_file_path*
setwd(paste(BASEDIR,"/WGBS",sep=""))
dir.create("results")
setwd(BASEDIR)

#read in WGBS methylated C and coverage data (bigTable)
tab <- fread("WGBS/data/GSE158927_BigTable.tsv.gz")

#convert WGBS data table into a Genomic Ranges object
tab <- tab[-grep("_", tab$chr),]
CpGs <- GRanges(tab$chr, IRanges(tab$position, width=1))
values(CpGs) <- as.data.frame(tab)[,-c(1:3)]
names(values(CpGs)) <- paste("X", names(values(CpGs)), sep='')
rm(tab)

build <- get("BSgenome.Hsapiens.UCSC.hg19")
seqlengths(CpGs) <- seqlengths(build)[seqlevels(CpGs)]

theseCpGs <- CpGs

#read in sample table
samples <- read.csv("metadata/WGBS_samplegroups.csv")

samples <- data.frame(Sample=samples$Sample,
                      C=paste("X", samples$Sample, ".C", sep=""),
                      cov=paste("X", as.character(samples$Sample), ".cov", sep=""),
                      Tissue=samples$Group,
                      row.names=samples$Sample,
                      stringsAsFactors=FALSE)

#define tumour samples
tum<-as.character(samples$Sample[which(samples$Tissue== "Cancer_Dead" | samples$Tissue== "Cancer_Alive")])

#define lethal and non-lethal samples
lethal<-as.character(samples$Sample[which(samples$Tissue== "Cancer_Dead" )])
nonlethal<-as.character(samples$Sample[which(samples$Tissue== "Cancer_Alive")])

#define normal samples
normals<-as.character(samples$Sample[which(samples$Tissue== "Normal")])

vectorstatus<-c(tum,normals)

#copy and rename sample table
samples2<-samples

#add colour column to samples table
samples2$col <- ""
samples2$col[samples2$Tissue %in% "Cancer_Dead"] <- "firebrick"
samples2$col[samples2$Tissue %in% "Cancer_Alive"] <- "forestgreen"
samples2$col[samples2$Tissue %in% "Normal"] <- "blue"
                     
#make separate data frame for coverage and methylated Cs
coverage <- as.data.frame(CpGs[,samples2$cov])
meth <- as.data.frame(CpGs[,samples2$C])

#make separate data frame of coverage and methylated Cs for each sample (required for DSS, make BSseqData)


X1601C <- data.frame(chr=coverage$seqnames, pos=coverage$start,
                    N=coverage$X1601C.cov, X=meth$X1601C.C)
X514C <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X514C.cov, X=meth$X514C.C)
X46C <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X46C.cov, X=meth$X46C.C)
X139C <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X139C.cov, X=meth$X139C.C)
X379C <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X379C.cov, X=meth$X379C.C)
X564C <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X564C.cov, X=meth$X564C.C)
X349C <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X349C.cov, X=meth$X349C.C)


X5237C <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X5237C.cov, X=meth$X5237C.C)
X174C <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X174C.cov, X=meth$X174C.C)
X1579C <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X1579C.cov, X=meth$X1579C.C)
X506C <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X506C.cov, X=meth$X506C.C)
X34C <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X34C.cov, X=meth$X34C.C)
X202C <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X202C.cov, X=meth$X202C.C)
X361C <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X361C.cov, X=meth$X361C.C)
X120C <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X120C.cov, X=meth$X120C.C)
                    
X448N <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X448N.cov, X=meth$X448N.C)
X1601N <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                     N=coverage$X1601N.cov, X=meth$X1601N.C)
X564N <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X564N.cov, X=meth$X564N.C)
X508N <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$X508N.cov, X=meth$X508N.C)
                                        


#make vector of which samples are in which group
lethal_idx<-which(samples$Tissue=="Cancer_Dead")
nonlethal_idx<-which(samples$Tissue=="Cancer_Alive")
normal_idx<-which(samples$Tissue=="Normal")

#samples for DML ('differentially methylated loci') analysis
tumsamps <- paste("X",samples$Sample[c(lethal_idx,nonlethal_idx)],sep="")

#DML analysis using DSS package
data_bsseq_dss <- makeBSseqData(lapply(tumsamps, get), tumsamps)
results_dml <- DMLtest(data_bsseq_dss, group1=tumsamps[lethal_idx], group2=tumsamps[nonlethal_idx], smoothing=F)

#calculate p-values and Benjamini-Hochberg correction
results_dml$pval <- 2*pnorm(-abs(results_dml$stat))
results_dml$fdr <- p.adjust(results_dml$pval, method="BH")

#save out DSS results
save(results_dml,file="WGBS/results/DSSresults.RData")

sig <- paste(results_dml$chr[results_dml$fdr < 0.05], results_dml$pos[results_dml$fdr < 0.05], sep = "_")

#use DML data to identify DMRs using DMRcate pacakge
anno <- cpg.annotate("sequencing", results_dml, analysis.type = "differential")
dmrs <- dmrcate(anno, C=50, min.cpgs = 5)
data(dmrcatedata)
res.ranges <- extractRanges(dmrs, genome="hg19")
res.ranges <- res.ranges[res.ranges$Stouffer < 0.05]

#export table of DMRs - used in Table 2 and Table S2
write.csv(as.data.frame(res.ranges), file="WGBS/results/Analysis_DMRs.csv")

res.ranges <- as.data.frame(res.ranges)

#export bedgraph file of location of DMRs - used as track in IGV in Figure S1
write.table(cbind(as.character(res.ranges$seqnames), res.ranges$start, res.ranges$end, res.ranges$meanbetafc), sep="\t",
            file="WGBS/results/analysisDMRs.bedGraph", quote=F, row.names=F, col.names=F)

#export bedgraph file of location and differential methylation of DMRs - used as track in IGV in Figure S1
#write.table(cbind(as.character(dmrs$input$CHR), dmrs$input$pos, dmrs$input$pos, dmrs$input$betafc), sep="\t",
#            file="WGBS/results/analysisdiffsites.bedGraph", quote=F, row.names=F, col.names=F)

#create DMR plots of top 50 DMRs
res.ranges <- extractRanges(dmrs, genome="hg19")
res.ranges <- res.ranges[res.ranges$Stouffer < 0.05]

#reformat theseCpGs ready for DMRplot
analysisCpGs <- GRanges(theseCpGs, values=values(theseCpGs)[,grep(paste(analysis, collapse="|"), colnames(values(theseCpGs)))])
colnames(values(analysisCpGs)) <- sub("values.", "", colnames(values(analysisCpGs)))

cols <- samples2[unique(sub("X", "", sub("\\..*", "", colnames(values(analysisCpGs))))), "col"]
names(cols) <- sub("Cancer_", "", samples2[unique(sub("X", "", sub("\\..*", "", colnames(values(analysisCpGs))))), "Tissue"])

#save out DMRplots (Figure 1d and Figure S1)
pdf("WGBS/plots/analysis_top50_DMRs.pdf", height=10, width=12)
sapply(1:50, function (x) {DMR.plot(res.ranges, dmr = x, CpGs=analysisCpGs, phen.col = cols, genome="hg19")})
dev.off()

