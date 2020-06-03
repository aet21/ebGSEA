#' @title Empirical Bayes Gene Set Enrichment Analysis
#'
#' @description
#' \pkg{ebGSEA} (Empirical Bayes Gene Set Enrichment Analysis) implements a GSEA designed for
#' Illumina Infinium Methylation beadchips, based on an empirical Bayes method to rank genes
#' based on their level of differential methylation, subsequently assessing enrichment of biological
#' terms using this ranked list.
#'
#' @details
#' \pkg{ebGSEA} leverages the evidence of differential methylation from all CpGs/probes
#' mapping to a given gene, to rank genes according to their overall level of differential
#' methylation. A key property of ebGSEA is that, like it does not favour genes
#' with high or low CpG/probe representation, thus avoiding the bias, whilst also
#' rendering the method sensitive enough to detect true biological enrichment. With genes
#' ranked by this empirical Bayes regression model, GSEA can subsequently be performed
#' using a non parametric Wilcoxon rank sum test or the known population median test
#' (KPMT), thus allowing GSEA to be performed in a threshold independent manner.
#'
#' @name ebGSEA-package
#'
#' @docType package
#'
#' @author Andrew E Teschendorff, Tianyu Zhu
#'
#' @references
#' Dong D, Tian Y, Zheng SC, Teschendorff AE.
#' \emph{ebGSEA: an improved Gene Set Enrichment Analysis method for Epigenome-Wide-Association Studies.}
#' BMC Bioinformatics (2019) 35(18):3514-3516.
#' doi:\href{https://doi.org/10.1093/bioinformatics/btz073}{
#' 10.1093/bioinformatics/btz073}.
#'
#' @examples
#' ### see example in tutorial
#'
NULL
