#' @name convertIDs
#'
#' @aliases convertIDs
#'
#' @title Convert between different gene IDs
#'
#' @description
#' Convert gene IDs based on `select` function in AnnotationDbi package
#'
#' @usage convertIDs(ids, from, to, db, ifMultiple=c("putNA", "useFirst"))
#'
#' @param ids
#' A vector of the gene IDs to convert from.
#'
#' @param from
#' The name of the gene identifier to convert from. All possible keys are returned by using the `keys` method.
#'
#' @param to
#' The name of the gene identifier to convert to. All possible keys are returned by using the `keys` method.
#'
#' @param db
#' The AnnotationDb object with the annoation of gene IDs to convert from and to.
#'
#' @param ifMultiple
#' If there are multiple hits for an input gene ID, whether to return 'NA' (ifMultiple = "putNA") or the first hit (ifMultiple = "useFirst").
#'
#' @return A vector of converted gene IDs.
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
#' #data("sgtm")
#' #rankEID.v <- rownames(sgt.m)
#' #sym.v <- convertIDs(rankEID.v, 'ENTREZID', 'SYMBOL', org.Hs.eg.db, ifMultiple="useFirst")
#'
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi select
#'
#' @export

convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )

  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }

  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}
