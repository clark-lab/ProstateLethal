#Script used to download genomic coordinates of CpG island data for use in analysis

#user to define file path for base directory
BASEDIR<-*user_defined_file_path*

#set working directory
setwd(BASEDIR)

#make folders
dir.create(file.path("WGBS","annotationdata"), recursive = TRUE)

#load library
library(rtracklayer)

session <- browserSession()
genome(session) <- "hg19"
query <- ucscTableQuery(session, "cpgIslandExt")
CpGislands <- track(query)
genome(CpGislands) <- NA

#CpG island shores & 5kb
CpGshores <- setdiff(resize(CpGislands, width(CpGislands)+4000, fix="center"), CpGislands)
CpG5kb <- resize(CpGislands, 5000)

save(CpGislands, CpGshores, CpG5kb, file="WGBS/annotationdata/CpGislands.Rdata")
export(CpGislands, "WGBS/annotationdata/CpGislands.bed")


