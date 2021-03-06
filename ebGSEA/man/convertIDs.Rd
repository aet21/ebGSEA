% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convertIDs.R
\name{convertIDs}
\alias{convertIDs}
\title{Convert between different gene IDs}
\usage{
convertIDs(ids, from, to, db, ifMultiple=c("putNA", "useFirst"))
}
\arguments{
\item{ids}{A vector of the gene IDs to convert from.}

\item{from}{The name of the gene identifier to convert from. All possible keys are returned by using the `keys` method.}

\item{to}{The name of the gene identifier to convert to. All possible keys are returned by using the `keys` method.}

\item{db}{The AnnotationDb object with the annoation of gene IDs to convert from and to.}

\item{ifMultiple}{If there are multiple hits for an input gene ID, whether to return 'NA' (ifMultiple = "putNA") or the first hit (ifMultiple = "useFirst").}
}
\value{
A vector of converted gene IDs.
}
\description{
Convert gene IDs based on `select` function in AnnotationDbi package
}
\examples{
#data("sgtm")
#rankEID.v <- rownames(sgt.m)
#sym.v <- convertIDs(rankEID.v, 'ENTREZID', 'SYMBOL', org.Hs.eg.db, ifMultiple="useFirst")

}
\references{
Dong D, Tian Y, Zheng SC, Teschendorff AE.
\emph{ebGSEA: an improved Gene Set Enrichment Analysis method for Epigenome-Wide-Association Studies.}
BMC Bioinformatics (2019) 35(18):3514-3516.
doi:\href{https://doi.org/10.1093/bioinformatics/btz073}{
10.1093/bioinformatics/btz073}.
}
\author{
Andrew E. Teschendorff, Tianyu Zhu
}
