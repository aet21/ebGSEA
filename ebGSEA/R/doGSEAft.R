#' @name doGSEAft
#'
#' @aliases doGSEAft
#'
#' @title GSEA with Fisher's Exact Test
#'
#' @description
#' Perform GSEA with Fisher's Exact Test to a group of selected genes
#'
#' @usage doGSEAft(selEID.v, ptw.ls, allEID.v, ncores = 4, minN = 5, adjPVth = 0.05)
#'
#' @param selEID.v
#' A vector of selected Entrez Gene ID.
#'
#' @param ptw.ls
#' Lists of Gene EntrezID in each pathway of interest. You can get the 8567 biological terms from Molecular Signatures Database by `data("MSigDB-28Feb14-data")`.
#'
#' @param allEID.v
#' A vector of the universal set of Entrez Gene ID which you select genes from.
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
#' @details
#' GSEA with Fisher's Exact Test is an extended enrichement test for user specified genes that doesn't require for ranks.
#'
#' @return \item{Rank(P)}{A matrix showing enriched pathways ranked by adjusted Fisher's Exact Test p values.
#' "nREP" is the number of genes in the pathway,
#' "nOVL" is the number of selected genes in the pathway,
#' "OR" is the odds ratio of Fisher's Exact Test,
#' "P" is the p value of Fisher's Exact Test,
#' "adjP" is the adjusted p value of Fisher's Exact Test (method='BH'),
#' "Genes" is all the selected genes in the pathway.}
#'
#' @return \item{Rank(OR)}{A matrix showing enriched pathways ranked by odds ratio. The columns are samely defined as in Rank(P).}
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
#' # topGSEAft.lm <- doGSEAft(selEID.v = sigEID.ls$selEID, ptw.ls = listEZ.lv, allEID.v = names(mapEIDto450k.lv), ncores = 1, adjPVth = 0.05)
#'
#' # Details can be found in tutorial
#'
#' @import parallel
#'
#' @export
#'
doGSEAft <- function(selEID.v,ptw.ls,allEID.v,ncores=4,minN=5,adjPVth=0.05){

  message(" Running Fisher's Exact Test... ")
  gseaFT.m <- matrix(unlist(mclapply(ptw.ls,gseaFTfn,selEID.v,allEID.v,minN,mc.cores=ncores)),ncol=5,byrow=TRUE)
  message(" Done ")
  colnames(gseaFT.m) <- c("nREP","nOVL","OR","P","Genes");
  rownames(gseaFT.m) <- names(ptw.ls);

  tmp.s <- sort(as.numeric(gseaFT.m[,4]),decreasing=FALSE,index.return=TRUE);
  sgseaFT.m <- gseaFT.m[tmp.s$ix,];
  padj.v <- p.adjust(as.numeric(sgseaFT.m[,4]),method="BH");
  sel.idx <- which(padj.v <= adjPVth);

  topGSEAft.m <- cbind(sgseaFT.m[sel.idx,1:4],padj.v[sel.idx],sgseaFT.m[sel.idx,5]);
  try(colnames(topGSEAft.m) <- c("nREP","nOVL","OR","P","adjP","Genes"))

  topGSEAft.lm <- list();
  topGSEAft.lm[[1]] <- topGSEAft.m;
  tmp.s <- sort(as.numeric(topGSEAft.m[,3]),decreasing=TRUE,index.return=TRUE);
  topGSEAft.lm[[2]] <- topGSEAft.m[tmp.s$ix,];
  names(topGSEAft.lm) <- c("Rank(P)","Rank(OR)");

  return(topGSEAft.lm);
}


### Function to perform GSEA using Fisher-test
#' @importFrom stats fisher.test
gseaFTfn <- function(termEID.v,selEID.v,allEID.v,minN=5){

  enrEID.v <- as.character(intersect(termEID.v,selEID.v));
  presEID.v <- intersect(termEID.v,allEID.v);
  nrep <- length(enrEID.v);
  if(nrep >= minN){
    #print(enrEID.v);
    tmp.m <- matrix(nrow=2,ncol=2);
    tmp.m[1,1] <- nrep;
    tmp.m[1,2] <- length(selEID.v)-nrep;
    tmp.m[2,1] <- length(intersect(allEID.v,termEID.v))-nrep;
    tmp.m[2,2] <- length(allEID.v) - tmp.m[1,1] - tmp.m[1,2] - tmp.m[2,1];
    ft.o <- fisher.test(tmp.m,alt="gr");

    mapSYM.v <- convertIDs(enrEID.v, "ENTREZID", "SYMBOL", org.Hs.eg.db,ifMultiple="useFirst");
    #      print(mapSYM.v);
    out.v <- c(sum(tmp.m[,1]),nrep,ft.o$est,ft.o$p.value,PasteVector(mapSYM.v));
  }else {
    out.v <- c(length(presEID.v),nrep,1,1,"NoGenes");
  }
  return(out.v);
}


### Auxiliary functions
PasteVector <- function(v){

  vt <- v[1];
  if(length(v) > 1){
    for(g in 2:length(v)){
      vt <- paste(vt,v[g],sep=" ")

    }
  }
  vt <- paste(vt," EnD",sep="");
  out.v <- sub(" EnD","",vt);
  out.v <- sub("NA , ","",out.v);
  out.v <- sub(" , NA","",out.v);
  out.v <- sub(" , NA , "," , ",out.v);
  return(out.v);
}

