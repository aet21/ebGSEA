#' @name selEIDfromSelCpG
#'
#' @aliases selEIDfromSelCpG
#'
#' @title Select significant genes from a selected set of CpGs
#'
#' @description
#' Select significant genes from a selected set of CpGs with binomial test.
#'
#' @usage selEIDfromSelCpG(selCpG.v, allCpG.v, pvth = 0.3/length(selCpG.v), array = c("450k","850k"))
#'
#' @param selCpG.v
#' A vector of user selected CpGs.
#'
#' @param allCpG.v
#' A vector of all CpGs the user select the CpGs from.
#'
#' @param pvth
#' P-value threshold to infer the number of selected CpGs mapped to a gene is significant or not in a binomial test. (default = 0.3/length(selCpG.v))
#'
#' @param array
#' Array type for the input CpGs. "450k" for Illumina HumanMethylation450 data and "850k" for Illumina MethylationEPIC data.
#'
#' @details
#' This function maps user specified CpGs to genes and return genes whose number of mapped CpGs is significant in binomial test.
#'
#' @return \item{selEID}{A vector of significant genes in Entrez Gene ID.}
#'
#' @return \item{selPV}{P-values of significant genes.}
#'
#' @return \item{allPV}{P-values of all mapped genes.}
#'
#' @references
#' Dong D, Tian Y, Zheng SC, Teschendorff AE.
#' \emph{ebGSEA: an improved Gene Set Enrichment Analysis method for Epigenome-Wide-Association Studies.}
#' BMC Bioinformatics (2019) 35(18):3514-3516.
#' doi:\href{https://doi.org/10.1093/bioinformatics/btz073}{
#' 10.1093/bioinformatics/btz073}.
#'
#' @author Andrew E. Teschendorff, Tianyu Zhu
#'
#' @examples
#' # data("SampleCpG")
#' # sigEID.ls <- selEIDfromSelCpG(selCpG.v = sampleCpG.v, allCpG.v = allCpG.v, array = "450k")
#'
#' @import parallel
#' @importFrom stats pbinom
#'
#' @export
#'
#'
selEIDfromSelCpG <- function(selCpG.v,allCpG.v,pvth=0.3/length(selCpG.v),array=c("450k","850k")){

  if(array=="450k"){
    message(" Mapping 450k probes to genes... ")
    data("dualmap450kEID");
    unqEID.v <- unique(unlist(map450ktoEID.lv[match(selCpG.v,names(map450ktoEID.lv))]));
    tmp.idx <- match(unqEID.v,names(mapEIDto450k.lv));
    tmp.l <- mapEIDto450k.lv[tmp.idx];
    message(" Done ")
  }else {
    message(" Mapping EPIC probes to genes... ")
    data("dualmap850kEID");
    unqEID.v <- unique(unlist(map850ktoEID.lv[match(selCpG.v,names(map850ktoEID.lv))]));
    tmp.idx <- match(unqEID.v,names(mapEIDto850k.lv));
    tmp.l <- mapEIDto850k.lv[tmp.idx];
    message(" Done ")
  }

  ncpg.v <- unlist(lapply(tmp.l,length));
  obsN.v <- unlist(lapply(lapply(tmp.l,intersect,selCpG.v),length));
  prob <- length(selCpG.v)/length(allCpG.v);

  pv.v <- pbinom(obsN.v,size=ncpg.v,prob,lower.tail=FALSE);
  names(pv.v) <- unqEID.v;
  sel.idx <- which(pv.v < pvth);
  selEID.v <- unqEID.v[sel.idx];
  return(list(selEID=selEID.v,selPV=pv.v[sel.idx],allPV=pv.v));
}
