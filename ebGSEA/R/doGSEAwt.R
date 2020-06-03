#' @name doGSEAwt
#'
#' @aliases doGSEAwt
#'
#' @title GSEA with Wilcoxon Rank Sum Test and the Known-Population Median Test
#'
#' @description
#' Perform GSEA with wilcoxon rank sum test and known-population test using the ranked gene list from global test.
#'
#' @usage doGSEAwt(rankEID.m, ptw.ls, ncores = 4, minN = 5, adjPVth = 0.05)
#'
#' @param rankEID.m
#' The resulted matrix from doGT function, with genes by row and ranked by statistics from global test. Rownames of the matrix should be gene EntrezID.
#'
#' @param ptw.ls
#' Lists of Gene EntrezID in each pathway of interest. You can get the 8567 biological terms from Molecular Signatures Database by `data("MSigDB-28Feb14-data")`.
#'
#' @param ncores
#' Number of cores used for parallel running. (default = 4)
#'
#' @param minN
#' For each pathway, the minium number of genes(i.e. available in the ranked gene list) to conduct GSEA. If less than this value, the p value of this pathway would be set 1. (default = 5)
#'
#' @param adjPVth
#' Adjusted p value threshold to infer a pathway to be significantly enriched or not. P value was derived from Wilcoxon rank sum test and adjusted with BH method. (default = 0.05)
#'
#' @details GSEA with Wilcoxon Rank Sum Test and the Known-Population Median Test is the second step of ebGSEA algorithm. Once the ranks of genes were derived from \emph{doGT}, enrichment of biological terms can be performed using either a standard one-tailed Wilcoxon rank sum test (WT), or a recently introduced more powerful version called Known-Population Median Test (KPMT). A group of enriched pathways can then be defined over adjusted pvalue from Wilcoxon rank sum test.
#'
#' @return \item{Rank(P)}{A matrix showing enriched pathways ranked by adjusted Wilcox test p values.
#' "nREP" is the number of mapped genes in the pathway,
#' "AUC" is the Area under curve of wilcox test,
#' "P(WT)" is the p-value of wilcox test,
#' "P(KPMT)" is the p-value of the Known-Population Median Test,
#' "adjP" is the adjusted p-value of wilcox test.}
#'
#' @return \item{Rank(AUC)}{A matrix showing enriched pathways ranked by AUC. The columns are defined samely as Rank(P).}
#'
#' @return \item{Genestat}{Lists of gene symbols in each enriched pathway. Each object contains the statistic and p-value from global test of each gene.}
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
#' # data("MSigDB-28Feb14-data")
#' # data("sgtm")
#' # topGSEA.lm <- doGSEAwt(rankEID.m = sgt.m, ptw.ls = listEZ.lv, ncores = 10, minN = 5, adjPVth = 0.05)
#'
#' @import parallel
#'
#' @export
#'
doGSEAwt <- function(rankEID.m,ptw.ls,ncores=4,minN=5,adjPVth=0.05){
  rankEID.v<-rownames(rankEID.m)
  message(" Running Wilcox Test and Known Population Median Test... ")
  gseaWT.m <- matrix(unlist(mclapply(ptw.ls,gseaWTfn,rankEID.v,mc.cores=ncores,minN=minN)),ncol=4,byrow=TRUE)
  message(" Done ")
  colnames(gseaWT.m) <- c("nREP","AUC","P(WT)","P(KPMT)");
  rownames(gseaWT.m) <- names(ptw.ls);

  tmp.s <- sort(gseaWT.m[,3],decreasing=FALSE,index.return=TRUE);
  sgseaWT.m <- gseaWT.m[tmp.s$ix,];
  padj.v <- p.adjust(sgseaWT.m[,3],method="BH");

  sel.idx <- which(padj.v <= adjPVth);
  topGSEAwt.lm <- list();
  sym.v<-convertIDs(rankEID.v, 'ENTREZID', 'SYMBOL', org.Hs.eg.db, ifMultiple="useFirst")
  del.idx<-which(is.na(sym.v))
  sym.v<-sym.v[-del.idx]
  rankEID.m<-rankEID.m[-del.idx,]
  rankEID.v<-rankEID.v[-del.idx]
  if(length(sel.idx)>1){

    topGSEAwt.m <- cbind(sgseaWT.m[sel.idx,],padj.v[sel.idx]);
    colnames(topGSEAwt.m) <- c("nREP","AUC","P(WT)","P(KPMT)","adjP");

    topGSEAwt.lm[[1]] <- topGSEAwt.m;
    tmp.s <- sort(topGSEAwt.m[,2],decreasing=TRUE,index.return=TRUE);
    topGSEAwt.lm[[2]] <- topGSEAwt.m[tmp.s$ix,];
    topGSEAwt.lm[[3]] <- list()
    for (i in 1:nrow(topGSEAwt.m)){
      EID.v<-intersect(ptw.ls[[match(rownames(topGSEAwt.m)[i],names(ptw.ls))]],rankEID.v)
      pathgene.m<-matrix(NA,nrow=length(EID.v),ncol=2)
      rownames(pathgene.m)<-sym.v[match(EID.v,rankEID.v)]
      colnames(pathgene.m)<-c('Pvalue','Statistic')
      pathgene.m[,1:2]<-rankEID.m[match(EID.v,rankEID.v),1:2]
      topGSEAwt.lm[[3]][[i]]<-pathgene.m
    }
    names(topGSEAwt.lm[[3]])<-rownames(topGSEAwt.m)
    names(topGSEAwt.lm) <- c("Rank(P)","Rank(AUC)","Genestat");
  }
  if (length(sel.idx)==1) {
    topGSEAwt.v <- as.vector(c(sgseaWT.m[sel.idx,],padj.v[sel.idx]));
    names(topGSEAwt.v) <- c("nREP","AUC","P(WT)","P(KPMT)","adjP");
    EID.v<-intersect(ptw.ls[[match(rownames(sgseaWT.m)[sel.idx],names(ptw.ls))]],rankEID.v)
    pathgene.m<-matrix(NA,nrow=length(EID.v),ncol=2)
    rownames(pathgene.m)<-sym.v[match(EID.v,rankEID.v)]
    colnames(pathgene.m)<-c('Pvalue','Statistic')
    pathgene.m[,1:2]<-rankEID.m[match(EID.v,rankEID.v),1:2]
    topGSEAwt.lm <- list("Rank(P)"=topGSEAwt.v,"Rank(AUC)"=topGSEAwt.v,"POI"=rownames(sgseaWT.m)[sel.idx],
                         "Genestat"=pathgene.m);
  }

  return(topGSEAwt.lm);
}


### Function to perform GSEA using wilcox-test and the known population median test (threshold independent)
#' @import kpmt
#' @importFrom stats wilcox.test
### gseaWTfn

gseaWTfn <- function(termEID.v,rankEID.v,minN=5){
  commonEID.v <- intersect(termEID.v,rankEID.v);
  nrep <- length(commonEID.v);
  if(length(commonEID.v)>=minN){
    otherEID.v <- setdiff(rankEID.v,termEID.v);
    match(commonEID.v,rankEID.v) -> rank1.idx;
    match(otherEID.v,rankEID.v) -> rank2.idx;
    wt.o <- wilcox.test(rank1.idx,rank2.idx,alt="less");
    pv <- wt.o$p.value;
    n1 <- length(rank1.idx);
    n2 <- length(rank2.idx);
    auc <- 1 - wt.o$stat/(n1*n2);
    ### now do kpmt
    pop.v <- 1:length(rankEID.v);
    names(pop.v) <- rankEID.v;
    obs.v <- commonEID.v;
    pvKPMT <- kpmt(pop=pop.v,obs=obs.v,tail="lower")[[4]];
    out.v <- c(nrep,auc,pv,pvKPMT);
  }else {
    out.v <- c(nrep,0,1,1);
  }
  return(out.v);
}

