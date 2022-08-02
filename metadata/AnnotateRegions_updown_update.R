#' annotateRegions_updown
#'
#' Annotate regions with respect to gene and CpG island annotations, and include whether DMRs are up or downstream of the closest TSS. Other ways to add this step may be to take advantage of the functions precede and follow in GenomicRanges (as well as distance ToNearest)
#'
#' @param reg The \code{GRanges} to be annotated
#' @param tx \code{GRanges} of transcripts, annotated with \code{gene_name}, \code{gene_id}, \code{tx_id} and \code{tx_type}
#' @param CpGislands \code{GRanges} of CpG island positions
#' @return \code{reg} with additional per region metadata:
#'  \item{nGeneTSS}{Number of TSSs overlapping}
#'  \item{nProtGeneTSS}{Number of TSSs overlapping (protein coding only)}
#'  \item{distanceTSS}{Distance to closest TSS}
#'  \item{reg$TSS_updown}{Whether DMR is upstream, downstream or overlapping which the closest TSS}
#'  \item{TSS}{Index of the closest TSS in \code{tx}}
#'  \item{tx_name}{Gencode tx_name}
#'  \item{tx_type}{Gencode tx_type}
#'  \item{gene_id}{Gencode gene_id}
#'  \item{gene_name}{Gencode gene_name}
#'  \item{distanceTSS_prot}{Distance to closest TSS (protein coding only)}
#'  \item{reg$TSS_prot_updown}{Whether DMR is upstream, downstream or overlapping which the closest TSS(protein coding only)}
#'  \item{TSS_prot}{Index of the closest TSS in tx (protein coding only)}
#'  \item{tx_name_prot}{Gencode tx_name (protein coding only)}
#'  \item{tx_type_prot}{Gencode tx_type (protein coding only)}
#'  \item{gene_id_prot}{Gencode gene_id (protein coding only)}
#'  \item{gene_name_prot}{ (Gencode gene_name protein coding only)}
#'  \item{distanceCpGi}{Distance to the closest CpG island}
#'  \item{CpGi}{Index of the closest CpG island in \code{CpGislands}}
#'  \item{promoter}{Proportion overlapping a \code{tx} promoter (+/- 2kb from TSS)}
#'  \item{genebody}{Proportion overlapping a \code{tx} genebody}
#'  \item{intergenic}{Proportion not overlapping \code{promoter} or \code{genebody}}
#'  \item{CpGisland}{Proportion overlapping CpG islands}
#'  \item{CpGshores}{Proportion overlapping CpG shores (+/- 2kb flanking CpG islands)}
#'  \item{nonCpG}{Proportion not overlaping \code{CpGisland} or \code{CpGshores}}
#'
#' @export
#'
#' @importFrom GenomeInfoDb seqlevels seqlengths
#' @importFrom GenomicRanges findOverlaps distanceToNearest resize setdiff flank
#' @importFrom IRanges IRanges
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
annotateRegions_updown <- function(reg, tx, CpGislands) {
    if (!all(c("gene_name", "gene_id", "tx_id", "tx_type") %in% names(values(tx)))) stop("Supplied tx does not contain all required columns, was it created by makeTx?")
    
    if (length(reg)==0) {
        warning("Supplied 'reg' is of length zero!")
        return(reg)
    }
    
    all.levels <- unique(c(seqlevels(reg), seqlevels(tx), seqlevels(CpGislands)))
    seqlevels(reg) <- seqlevels(tx) <- seqlevels(CpGislands) <- all.levels
    
    stopifnot(all(!is.na(seqlengths(tx))))
    # Want to do most analysis for "just" protein coding as well as all genes
    tx2 <- tx[tx$tx_type=="protein_coding"]
    
    # Number of different genes TSSs within the DMR
    reg.TSS <- as.matrix(findOverlaps(reg, resize(tx, 1, fix="start")))
    reg.TSS <- tapply(tx$gene_name[reg.TSS[,2]], reg.TSS[,1], function(x) length(unique(x)))
    reg$nGeneTSS <- 0
    reg$nGeneTSS[as.integer(names(reg.TSS))] <- unname(reg.TSS)
    #just for protein coding genes
    reg.TSS <- as.matrix(findOverlaps(reg, resize(tx2, 1, fix="start")))
    reg.TSS <- tapply(tx2$gene_name[reg.TSS[,2]], reg.TSS[,1], function(x) length(unique(x)))
    reg$nProtGeneTSS <- 0
    reg$nProtGeneTSS[as.integer(names(reg.TSS))] <- unname(reg.TSS)
    
    
    # Distance to closest TSS
    d2Nlist <- (distanceToNearest(reg, resize(tx, 1, fix="start")))
    nearest<-tx[d2Nlist@to]
    OVERLAP<-findOverlaps(reg,resize(nearest, 1, fix="start"))@to
    absdist<-(resize(reg, 1, fix="start")@ranges@start)-(resize(nearest, 1, fix="start")@ranges@start)
    upstream1<-which(absdist<0 & strand(nearest)=="+")
    upstream2<-which(absdist>0 & strand(nearest)=="-")
    upstream<-c(upstream1,upstream2)
    updown<-vector(length=length(absdist))
    updown[upstream]<-"UP"
    updown[-upstream]<-"DOWN"
    updown[OVERLAP]<-"Over"
    
    reg.dist <- as.data.frame(d2Nlist)
    reg$distanceTSS <- reg.dist$distance
    reg$TSS <- reg.dist$subjectHits
    reg$TSSupdown<-updown
    reg$tx_id <- tx$tx_id[reg$TSS]
    reg$tx_type <- tx$tx_type[reg$TSS]
    reg$gene_id <- tx$gene_id[reg$TSS]
    reg$gene_name <- tx$gene_name[reg$TSS]
    
    # Distance to closest protein-coding TSS
    d2Nlist <- (distanceToNearest(reg, resize(tx2, 1, fix="start")))
    nearest<-tx2[d2Nlist@to]
    OVERLAP<-findOverlaps(reg,resize(nearest, 1, fix="start"))@to
    absdist<-(resize(reg, 1, fix="start")@ranges@start)-(resize(nearest, 1, fix="start")@ranges@start)
    upstream1<-which(absdist<0 & strand(nearest)=="+")
    upstream2<-which(absdist>0 & strand(nearest)=="-")
    upstream<-c(upstream1,upstream2)
    updown<-vector(length=length(absdist))
    updown[upstream]<-"UP"
    updown[-upstream]<-"DOWN"
    updown[OVERLAP]<-"Over"
    
    reg.dist <- as.data.frame(d2Nlist)
    reg$distanceTSS_prot <- reg.dist$distance
    reg$TSS_prot <- reg.dist$subjectHits
    reg$TSS_prot_updown<-updown
    reg$tx_id_prot <- tx2$tx_id[reg$TSS_prot]
    reg$gene_id_prot <- tx2$gene_id[reg$TSS_prot]
    reg$gene_name_prot <- tx2$gene_name[reg$TSS_prot]
    
    # Distance to closest CpG island
    reg.dist <- as.data.frame(distanceToNearest(reg, CpGislands))
    reg$distanceCpGi <- reg.dist$distance
    reg$CpGi <- reg.dist$subjectHits
    
    # % promoter/genebody/intergenic (2kb up and down)
    promotersGR <- suppressWarnings(reduce(strip(resize(resize(tx, 1, fix="start"), 4000, fix="center"))))
    genebodyGR <- suppressWarnings(setdiff(reduce(strip(tx)), promotersGR))
    genomeGR <- GRanges(seqlevels(tx), IRanges(1, seqlengths(tx)))
    intergenicGR <- suppressWarnings(setdiff(setdiff(genomeGR, genebodyGR), promotersGR))
    
    reg$promoter <- coverageRatio(reg, promotersGR)
    reg$genebody <- coverageRatio(reg, genebodyGR)
    reg$intergenic <- coverageRatio(reg, intergenicGR)
    
    # % protein coding promoter/genebody/intergenic (2kb up and down)
    promotersGR <- suppressWarnings(reduce(strip(resize(resize(tx2, 1, fix="start"), 4000, fix="center"))))
    genebodyGR <- suppressWarnings(setdiff(reduce(strip(tx2)), promotersGR))
    genomeGR <- suppressWarnings(GRanges(seqlevels(tx2), IRanges(1, seqlengths(tx2))))
    intergenicGR <- suppressWarnings(setdiff(setdiff(genomeGR, genebodyGR), promotersGR))
    
    reg$promoter_prot <- coverageRatio(reg, promotersGR)
    reg$genebody_prot <- coverageRatio(reg, genebodyGR)
    reg$intergenic_prot <- coverageRatio(reg, intergenicGR)
    
    # % CpG island/CpG shore/nonCpG
    reg$CpGisland <- coverageRatio(reg, CpGislands)
    CpGshores <- suppressWarnings(setdiff(reduce(resize(CpGislands, width(CpGislands)+4000, fix="center")), CpGislands))
    reg$CpGshores <- coverageRatio(reg, CpGshores)
    reg$nonCpG <- suppressWarnings(coverageRatio(reg, setdiff(setdiff(genomeGR, CpGislands), CpGshores)))
    
    reg
}