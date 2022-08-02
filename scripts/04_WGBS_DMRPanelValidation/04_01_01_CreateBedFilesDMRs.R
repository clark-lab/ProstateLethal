#Script for creating bedgraph files for WGBS data for lethal, non-lethal and normal groups for IGV tracks
#Generates bedgraph files for use in IGV to create Figure S1

#set user defined working directory
BASEDIR<-*user_defined_file_path*
setwd(BASEDIR)

#load in DMR methylation data generated in '02_02_DMRmean_heatmap.R'
load("WGBS/results/methGR_analysisCpGsInDMRs_10000bp.RData") #an_res_CpGs
load("WGBS/results/methGR_analysisDMRs.RData") #an_DMRs

#read in DMRs that were selected for validation
selecDMR<-read.table("metadata/SelectedRegions.bed")

#identify the 18 DMRs used in validation
matchup<-match(selecDMR$V2,an_DMRs@ranges@start)

#restrict methylation data to just the 18 DMRs used in validation
res_selec<-an_res_CpGs[matchup]

#define order of samples used in analysis
tum <- c("X139C", "X1601C", "X349C", "X379C", "X46C", "X514C", "X564C", "X120C", "X1579C", "X174C", "X202C", "X34C", "X361C", "X506C", "X5237C")

#define lethal and non-lethal samples
lethal<-tum[1:7]
nonlethal<-tum[8:15]

#define normal samples
normals<-c("X448N","X1601N","X564N","X508N")

#separate methylation data by sample group
idx_leth<-vector()
for(i in 1:length(lethal)){
  idx_leth[i]<-grep(lethal[i],colnames(res_selec[[1]]@elementMetadata))
}

idx_nonleth<-vector()
for(i in 1:length(nonlethal)){
  idx_nonleth[i]<-grep(nonlethal[i],colnames(res_selec[[1]]@elementMetadata))
}

idx_norm<-vector()
for(i in 1:length(normals)){
  idx_norm[i]<-grep(normals[i],colnames(res_selec[[1]]@elementMetadata))
}

#create matrix of methylation at CpG in all DMRs, with columns ordered by patient group
DMRlistmat<-list()
for(i in 1:length(res_selec)){
  
  DMRlistmat[[i]]<-as.matrix(res_selec[[i]]@elementMetadata)
}

DMRlistmat_ord<-list()
for(i in 1:length(DMRlistmat)){
  DMRlistmat_ord[[i]]<-DMRlistmat[[i]][,c(idx_leth,idx_nonleth,idx_norm)]
}

allCpGs<-do.call(rbind,DMRlistmat_ord)
dim(allCpGs)

allres<-do.call("c", res_selec)

allres@elementMetadata<-data.frame(allCpGs)

#calculate the average methylation of the lethal group at 18 DMRs
df_lethal2<- data.frame(seqnames=seqnames(allres),
                             starts=start(allres)-1,
                             ends=end(allres)+1,
                            scores=rowMeans(allCpGs[,1:7],na.rm=T),
                             strands=strand(allres))



#calculate the average methylation of the non-lethal group at 18 DMRs
df_nonlethal2<- data.frame(seqnames=seqnames(allres),
                       starts=start(allres)-1,
                       ends=end(allres)+1,
                       scores=rowMeans(allCpGs[,8:15],na.rm=T),
                       strands=strand(allres))

#calculate the average methylation of the normal group at 18 DMRs
df_normal2<- data.frame(seqnames=seqnames(allres),
                          starts=start(allres)-1,
                          ends=end(allres)+1,
                            scores=rowMeans(allCpGs[,16:19],na.rm=T),
                          strands=strand(allres))

# export bedgraph files of WGBS data for lethal, non-lethal and normal groups for IGV tracks to create Figure S1
write.table(df_lethal2,
            file="WGBS/results/LethalAv_NArm.bedGraph",
            quote=F, sep="\t", row.names=F, col.names=F,na="0")

write.table(df_nonlethal2,
            file="WGBS/results/NonlethalAv_NArm.bedGraph",
            quote=F, sep="\t", row.names=F, col.names=F,na="0")

write.table(df_normal2,
            file="WGBS/results/NormalAv_NArm.bedGraph",
            quote=F, sep="\t", row.names=F, col.names=F,na="0")

