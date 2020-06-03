#' @name doGT
#'
#' @aliases doGT
#'
#' @title Empirical Bayes Global Test
#'
#' @description
#' Function to assess the overall level of differential methylation using all the probes mapping to a gene.
#'
#' @usage doGT(pheno.v, data.m, model = c("linear"), array = c("450k", "850k"), ncores = 4)
#'
#' @param pheno.v
#' A vector of phenotype information, must be matched to columns of the input beta matrix.
#'
#' @param data.m
#' A matrix of beta values with probes by row and samples by column. Missing values shoud be excluded.
#'
#' @param model
#' The regression model for global test. Default is "linear".
#'
#' @param array
#' Array type for the input data. "450k" for Illumina HumanMethylation450 data and "850k" for Illumina MethylationEPIC data.
#'
#' @param ncores
#' Number of cores used for parallel running. (default = 4)
#'
#' @details
#' Global test is the first step of ebGSEA algorithm. ebGSEA ranks genes according to their overall level of differential methylation by adapting the global test from \emph{Geoman et al(2006)}, which can be interpreted as an empirical Bayes generalized regression model. The global test evaluates whether DNA methylation patterns of CpGs mapping to a given gene \emph{g} differ significantly between two phenotypes.
#'
#' @return A matrix with genes in row ranked by statistic from global test.
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
#' # sgt.m <- doGT(pheno.v, data.m, array = c("450k"), ncores = 10)
#' # topGSEA.lm <- doGSEAwt(rankEID.m = sgt.m, ptw.ls = listEZ.lv, ncores = 10, minN = 5, adjPVth = 0.05)
#'
#'
#' @import parallel
#' @import globaltest
#' @importFrom utils data
#'
#' @export
#'
doGT <- function(pheno.v,data.m,model=c("linear"),array=c("450k","850k"),ncores=4){
  if(array=="450k"){
    message(" Mapping 450k probes to genes... ")
    data("dualmap450kEID");
    subsets <- mclapply(mapEIDto450k.lv,intersect,rownames(data.m),mc.cores = ncores);
    message(" Done ")
  }else {
    message(" Mapping EPIC probes to genes... ")
    data("dualmap850kEID");
    subsets <- mclapply(mapEIDto850k.lv,intersect,rownames(data.m),mc.cores = ncores);
    message(" Done ")
  }
  nrep.v <- unlist(lapply(subsets,length));
  selG.idx <- which(nrep.v>0);
  message(" Running Global Test... ")
  gt.o <- gt(response=pheno.v,alternative=t(data.m),model=model,directional = FALSE, standardize = FALSE, permutations = 0, subsets=subsets[selG.idx],trace=F);
  message(" Done ")
  resGT.m <- as.matrix(result(gt.o));
  tmp.s <- sort(resGT.m[,2],decreasing=TRUE,index.return=TRUE);
  sresGT.m <- resGT.m[tmp.s$ix,];
  return(sresGT.m);
}
