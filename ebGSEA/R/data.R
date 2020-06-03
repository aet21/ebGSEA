#' Biological terms from Molecular Signatures Database 28Feb14
#'
#' The 8567 biological terms from \href{https://www.gsea-msigdb.org/gsea/msigdb/}{Molecular Signatures Database} 28Feb14
#'
#' \itemize{
#'   \item listEZ.lv : Gene sets in NCBI (Entrez) Gene IDs for each biological term
#'   \item listG.lv : Gene stes in gene symbols for each biological term
#'   \item listclassALL.lv : The type of biological term defined by MSiDB
#' }
#'
#' @docType data
#' @keywords biological terms
#' @name MSigDB-28Feb14-data
#' @usage data("MSigDB-28Feb14-data")
#' @format Two lists with 8567 items and a vector with 8567 items
#' @references
#' Dong D, Tian Y, Zheng SC, Teschendorff AE.
#' \emph{ebGSEA: an improved Gene Set Enrichment Analysis method for Epigenome-Wide-Association Studies.}
#' BMC Bioinformatics (2019) 35(18):3514-3516.
#' doi:\href{https://doi.org/10.1093/bioinformatics/btz073}{
#' 10.1093/bioinformatics/btz073}.
#'
#' Subramanian A, Tamayo P, Mootha VK, et al.
#' \emph{Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles.}
#' Proc Natl Acad Sci U S A (2005) 102(43):15545‐15550.
#' doi:\href{https://doi.org/10.1073/pnas.0506580102}{
#' 10.1073/pnas.0506580102}.
#'
NULL

#' Illumina HM450k probes annotation to Entrez Gene ID
#'
#' This annotation file is derived from Illumina HM450k annotation
#'
#' \itemize{
#'   \item map450ktoEID.lv : A list mapping 450k probes to Entrez Gene ID
#'   \item mapEIDto450k.lv : A list mapping Entrez Gene ID to 450k probes
#' }
#'
#' @docType data
#' @keywords annotation
#' @name dualmap450kEID
#' @usage data("dualmap450kEID")
#' @format A list with 485577 items and a list with 18649 items
#' @references
#' Dong D, Tian Y, Zheng SC, Teschendorff AE.
#' \emph{ebGSEA: an improved Gene Set Enrichment Analysis method for Epigenome-Wide-Association Studies.}
#' BMC Bioinformatics (2019) 35(18):3514-3516.
#' doi:\href{https://doi.org/10.1093/bioinformatics/btz073}{
#' 10.1093/bioinformatics/btz073}.
#'
NULL

#' Illumina EPIC probes annotation to Entrez Gene ID
#'
#' This annotation file is derived from Illumina annotation file
#'
#' \itemize{
#'   \item map850ktoEID.lv : A list mapping 850k probes to Entrez Gene ID
#'   \item mapEIDto850k.lv : A list mapping Entrez Gene ID to 850k probes
#' }
#'
#' @docType data
#' @keywords annotation
#' @name dualmap850kEID
#' @usage data("dualmap850kEID")
#' @format A list with 867531 items and a list with 24357 items
#' @references
#' Dong D, Tian Y, Zheng SC, Teschendorff AE.
#' \emph{ebGSEA: an improved Gene Set Enrichment Analysis method for Epigenome-Wide-Association Studies.}
#' BMC Bioinformatics (2019) 35(18):3514-3516.
#' doi:\href{https://doi.org/10.1093/bioinformatics/btz073}{
#' 10.1093/bioinformatics/btz073}.
#'
NULL

#' Intermediate result of global test in Tutorial
#'
#' The result of doing global test to buccal swab samples
#'
#' \itemize{
#'   \item sgt.m : The resulting matrix from `doGT` function with the tutorial sample
#' }
#'
#' @docType data
#' @keywords tutorial
#' @name sgtm
#' @usage data("sgtm")
#' @format A matrix with 18618 rows and 5 columns
#' @references
#' Dong D, Tian Y, Zheng SC, Teschendorff AE.
#' \emph{ebGSEA: an improved Gene Set Enrichment Analysis method for Epigenome-Wide-Association Studies.}
#' BMC Bioinformatics (2019) 35(18):3514-3516.
#' doi:\href{https://doi.org/10.1093/bioinformatics/btz073}{
#' 10.1093/bioinformatics/btz073}.
#'
NULL

#' Sample CpGs for fisher's exact test in Tutorial
#'
#' The differentialy methylated cytosines associated with smoking pack-years identified from buccal swab dataset
#'
#' \itemize{
#'   \item sampleCpG.v : The pre-selected differentially methylated cytosines
#'   \item allCpG.v : All the CpGs from buccal swab dataset
#' }
#'
#' @docType data
#' @keywords tutorial
#' @name SampleCpG
#' @usage data("SampleCpG")
#' @format A vector with length 40626 and a vector with length 484272
#' @references
#' Dong D, Tian Y, Zheng SC, Teschendorff AE.
#' \emph{ebGSEA: an improved Gene Set Enrichment Analysis method for Epigenome-Wide-Association Studies.}
#' BMC Bioinformatics (2019) 35(18):3514-3516.
#' doi:\href{https://doi.org/10.1093/bioinformatics/btz073}{
#' 10.1093/bioinformatics/btz073}.
#'
NULL
