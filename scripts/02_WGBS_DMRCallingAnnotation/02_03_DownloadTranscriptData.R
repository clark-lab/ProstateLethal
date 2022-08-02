#Script used to download hg19 transcript data for use in analysis

#user to define file path for base directory
BASEDIR<-*user_defined_file_path*

#set working directory
setwd(BASEDIR)

#load libraries
library(aaRon)
library(BSgenome.Hsapiens.UCSC.hg19)

#make tx object from an Ensembl GTF, required for annotating DMRs
tx <- makeTx("Homo_sapiens.GRCh37.75.gtf", Hsapiens, "source")
save(tx, file="WGBS/annotationdata/hg19.tx.Rdata")
