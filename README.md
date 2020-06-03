---
title: "Introduction to ebGSEA"
author:
- name: "Andrew E. Teschendorff"
  affiliation: 
  - CAS Key Lab of Computational Biology, PICB, SINH
- name: "Tianyu Zhu"
  affiliation:
  - CAS Key Lab of Computational Biology, PICB, SINH
date: "2020-06-03"
package: ebGSEA
output:
  BiocStyle::html_document:
    theme: readable
bibliography: ebGSEA.bib
vignette: >
  %\VignetteIndexEntry{Empirical Bayes Gene Set Enrichment Analysis}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
  

```
## Loading required package: BiocStyle
```

```
## 
## Attaching package: 'BiocStyle'
```

```
## The following objects are masked from 'package:rmarkdown':
## 
##     html_document, md_document, pdf_document
```

# Introduction

Gene Set Enrichment Analysis is one of the most common tasks in the analysis of omic data, and is critical for biological interpretation. In the context of Epigenome Wide Association Studies (EWAS), which typically rank individual cytosines according to the level of differential methylation, enrichment analysis of biological pathways is challenging due to differences in CpG/probe density between genes. Here we use an empirical Bayes Gene Set Enrichment Analysis (ebGSEA) algorithm, which does not rank CpGs but genes according to the overall level of differential methylation of its CpGs/probes, allowing unbiased and sensitive detection of enriched pathways. ebGSEA can be applied to any EWAS that uses Illumina HM450k and EPIC beadarrays. For details please refer to our publication listed at the end of this tutorial[@dong2019ebgsea].

# Tutorial Example

## Differentially methylated genes in a buccal swab EWAS of smoking

To illustrate the functions, we use a subset of a HM450k buccal swab dataset[@teschendorff2015correlation], containing 325 buccal swab samples with smoking-pack-years as phenotype, and 7933 CpG probes. The beta value matrix has been processed, with samples by column and probes by row.


```r
load("ebgseaDATA.Rd");
ls();
```

```
## [1] "dataSMK.m"  "phenoSMK.v"
```

```r
dim(dataSMK.m)
```

```
## [1] 7933  325
```

To rank genes according to the levels of differential methylation of probes mapping to each gene, we apply the global test[@goeman2004a] function. Here are the inputs for the `doGT` function:

* pheno.v: A vector containing the phenotype information, matched to the columns of the input DNAm data matrix. In our case this contains the smoking pack years information.
* data.m: The matrix of DNAm beta values with probes along rows and samples along columns. Missing values are not allowed and should be imputed beforehand.
* model: The regression model for the global test. Default is “linear”.
* array: Array type for the input data. “450k” for Illumina HumanMethylation450 data and “850k” for Illumina MethylationEPIC data.
* ncores: Number of cores used for parallel running. (default = 4)

To run `doGT` :

```r
library(ebGSEA)
```

```
## 
```

```r
sgt.m <-doGT(phenoSMK.v,dataSMK.m,array="450k",ncores=4)
```

```
##  Mapping 450k probes to genes...
```

```
##  Done
```

```
##  Running Global Test...
```

```
##  Done
```

```r
head(sgt.m)
```

```
##            p-value Statistic Expected   Std.dev #Cov
## 8140  2.897869e-14 15.856766 0.308642 0.4210446    3
## 2877  3.086939e-12 13.997658 0.308642 0.4358116    1
## 4171  3.168428e-10 10.804754 0.308642 0.4080425    2
## 9639  8.863525e-10 10.160649 0.308642 0.4025246    4
## 84623 5.816328e-09  9.976524 0.308642 0.4358116    1
## 26166 1.346266e-08  9.519885 0.308642 0.4358116    1
```

The output `sgt.m` is a matrix with rows labeling genes, ordered according to their overall level of differential methylation. The last column gives us the number of CpGs/probes mapping to the given gene.

## Pathway enrichment analysis with Wilcox Rank Sum Test and Known-Population Meian Test

The next step is to apply the function `doGSEAwt` to do a pathway enrichment analysis in a threshold independent manner, using the Wilcox test (WT) and the Known-Population Median Test (KPMT). Here are the input parameters:

* rankEID.m: The matrix output object from `doGT` function, with rows labeling genes. Rownames of the matrix should be Gene Entrez IDs.
* ptw.ls: List of vectors consisting of Gene EntrezIDs of genes pathways of interest. For a pathway-database we use 8567 biological terms from the Molecular Signatures Database[@subramanian2005gene], by invoking `data("MSigDB-28Feb14-data")`.
* ncores: Number of cores used for parallel running. (default = 4)
* minN: For each pathway, the minium number of genes(i.e. available in the ranked gene list) to conduct GSEA. If less than this value, the P-value of this pathway would be set 1. (default = 5)
* adjPVth: Adjusted P-value threshold to declare a pathway to be significantly enriched. P-value was derived from Wilcoxon rank sum test and adjusted with BH method. (default = 0.05)


```r
data("MSigDB-28Feb14-data")
topGSEA.lm <- doGSEAwt(rankEID.m = sgt.m, ptw.ls = listEZ.lv, ncores = 4, minN = 5,adjPVth = 0.05)
```

```
##  Running Wilcox Test and Known Population Median Test...
```

```
##  Done
```

```
## 'select()' returned 1:1 mapping between keys and columns
```
The output is a list of three objects: 

* Rank(P): A matrix showing enriched pathways ranked by adjusted Wilcox test P-values.
* Rank(AUC): A matrix showing enriched pathways ranked by AUC.
* Genestat: Lists of gene symbols in each enriched pathway. Each object contains the statistic and P-value from global test of each gene.


```r
topGSEA.lm[[1]]
```

```
##                                        nREP       AUC        P(WT)      P(KPMT)
## DODD_NASOPHARYNGEAL_CARCINOMA_UP       1465 0.5565279 9.056876e-11 6.538912e-10
## MARTENS_TRETINOIN_RESPONSE_UP           213 0.5922184 2.470191e-06 3.254273e-04
## BENPORATH_ES_WITH_H3K27ME3              360 0.5719432 2.498177e-06 2.218392e-06
## REGULATION_OF_CELL_MIGRATION              9 0.9127553 9.123491e-06 6.725557e-05
## SCHAEFFER_PROSTATE_DEVELOPMENT_48HR_UP  165 0.5967493 1.131655e-05 5.984702e-04
##                                                adjP
## DODD_NASOPHARYNGEAL_CARCINOMA_UP       7.759026e-07
## MARTENS_TRETINOIN_RESPONSE_UP          7.133962e-03
## BENPORATH_ES_WITH_H3K27ME3             7.133962e-03
## REGULATION_OF_CELL_MIGRATION           1.938977e-02
## SCHAEFFER_PROSTATE_DEVELOPMENT_48HR_UP 1.938977e-02
```

The columns of this matrix indicate the following:

* nREP: Number of genes mapped in this pathway and present on array
* AUC: Area under curve from Wilcox test
* P(WT): P-value from Wilcox Test
* P(KPMT): P-value from Known Population Median Test
* adjP: Adjusted P-value for each pathway, using BH method

We can see that the top-ranked biological term is related to nasopharyngeal carcinoma, which is a sensible result because the nasopharynx is exposed to smoke carcinogens, and smoking is indeed a major risk factor for this type of cancer.

The third list Genestat consists of the statistic and P-value from global test, of each gene in a specific enriched pathway. To see the genes in the first pathway:

```r
head(topGSEA.lm$Genestat[[1]])
```

```
##             Pvalue Statistic
## RIBC1    0.3178248 0.3089405
## KIAA1614 0.7876981 0.0224852
## GPRC5C   0.3604448 0.3082687
## PRKAA2   0.2204181 0.4645438
## CLIC6    0.4085130 0.2115970
## TPRG1L   0.2138787 0.4707674
```

We can plot the AUC and adjP for each enriched pathway:

```r
plot(x = topGSEA.lm[[1]][,2], y = -log10(topGSEA.lm[[1]][,5]), xlab ="AUC", ylab = "-log10(adjP)", main = "AUC and adjP for each enriched pathway", pch = 21, bg = "red");
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

## Pathway enrichment analysis with Fisher's Exact Test

Additionally, ebGSEA allows users to do GSEA in the more traditional way on a group of specified CpGs/Genes with Fisher's Exact Test, without the need of ranking genes.

If your input is a set of CpGs (say a list of top-ranked DMCs), you may use `selEIDfromSelCpG` function to first derive a corresponding list of genes. This function implements a binomial test to determine the statistical significance of a gene being selected on account of the number of CpGs mapping to the gene. Here are the input parameters:

* selCpG.v: A vector of user selected CpGs.
* allCpG.v: A vector representing all CpGs in the input DNAm data matrix.
* pvth: P-value threshold to infer the number of selected CpGs mapped to a gene is significant or not in a binomial test. (default = 0.3/length(selCpG.v))
* array: Array type for the input CpGs. "450k" for Illumina HumanMethylation450 data and "850k" for Illumina MethylationEPIC data.

To illustrate this function, we import a set of selected CpGs derived from a more complete buccal swab dataset (i.e. over 400k CpGs) by invoking `data("sampleCpG")`. This set of 40626 CpGs exhibit differential methylation associated with smoking pack-years.


```r
data("SampleCpG")
sigEID.ls <- selEIDfromSelCpG(selCpG.v = sampleCpG.v, allCpG.v = allCpG.v, array = "450k")
```

```
##  Mapping 450k probes to genes...
```

```
##  Done
```

The group of CpGs are significantly mapped to 255 genes. 

```r
summary(sigEID.ls)
```

```
##        Length Class  Mode     
## selEID  255   -none- character
## selPV   255   -none- numeric  
## allPV  9153   -none- numeric
```

Now we can apply the `doGSEAft` function to do a pathway enrichment analysis with Fisher's Exact Test. Here are the input parameters:

* selEID.v: A vector of selected Entrez Gene ID.
* ptw.ls: Lists of Gene EntrezID in each pathway of interest. You can get the 8567 biological terms from Molecular Signatures Database by `data("MSigDB-28Feb14-data")`.
* allEID.v: A vector of the universal set of Entrez Gene ID which you select genes from.
* ncores: Number of cores used for parallel running. (default = 4)
* minN: For each pathway, the minium number of genes(i.e. available in the ranked gene list) to conduct GSEA. If less than this value, the p value of this pathway would be set 1. (default = 5)
* adjPVth: Adjusted p value threshold to infer a pathway to be significantly enriched or not. P value was derived from Wilcoxon rank sum test and adjusted with BH method. (default = 0.05)


```r
topGSEAft.lm <- doGSEAft(selEID.v = sigEID.ls$selEID, ptw.ls = listEZ.lv, allEID.v = names(mapEIDto450k.lv), ncores = 1, adjPVth = 0.05)
```
The output of `topGSEAft.lm` consists of the following items:

* Rank(P): A matrix showing enriched pathways ranked by adjusted Fisher's Exact Test P-values. "nREP" is the number of genes in the pathway, "nOVL" is the number of selected genes in the pathway, "OR" is the odds ratio of Fisher's Exact Test, "P" is the P-value of Fisher's Exact Test, "adjP" is the adjusted P-value of Fisher's Exact Test (method='BH'), "Genes" is all the selected genes in the pathway.
* Rank(OR): A matrix showing enriched pathways ranked by odds ratio. The columns are samely defined as in Rank(P).

There are 79 enriched pathways identified by Fisher's exact test, each with 6 features:

```r
summary(topGSEAft.lm)
```

```
##          Length Class  Mode     
## Rank(P)  474    -none- character
## Rank(OR) 474    -none- character
```


```r
head(topGSEAft.lm$`Rank(P)`)
```

```
##                                               nREP   nOVL OR                
## JAEGER_METASTASIS_DN                          "250"  "25" "8.774214905728"  
## CHARAFE_BREAST_CANCER_BASAL_VS_MESENCHYMAL_UP "113"  "16" "12.6237084427608"
## ONDER_CDH1_TARGETS_2_DN                       "442"  "26" "4.90504233115453"
## HATADA_METHYLATED_IN_LUNG_CANCER_UP           "365"  "22" "4.96807269041407"
## HOX_GENES                                     "64"   "10" "13.8598927020685"
## CAGGTG_V$E12_Q6                               "2314" "64" "2.40406322547187"
##                                               P                     
## JAEGER_METASTASIS_DN                          "8.37864316495949e-15"
## CHARAFE_BREAST_CANCER_BASAL_VS_MESENCHYMAL_UP "3.17153270879009e-12"
## ONDER_CDH1_TARGETS_2_DN                       "4.47046000236504e-10"
## HATADA_METHYLATED_IN_LUNG_CANCER_UP           "6.7330366346862e-09" 
## HOX_GENES                                     "1.5174248583798e-08" 
## CAGGTG_V$E12_Q6                               "1.90038100328723e-08"
##                                               adjP                  
## JAEGER_METASTASIS_DN                          "7.1779835994208e-11" 
## CHARAFE_BREAST_CANCER_BASAL_VS_MESENCHYMAL_UP "1.35852603581023e-08"
## ONDER_CDH1_TARGETS_2_DN                       "1.27661436134204e-06"
## HATADA_METHYLATED_IN_LUNG_CANCER_UP           "1.44204812123392e-05"
## HOX_GENES                                     "2.59995575234795e-05"
## CAGGTG_V$E12_Q6                               "2.71342734252695e-05"
##                                               Genes                                                                                                                                                                                                                                                                                                                                                                                        
## JAEGER_METASTASIS_DN                          "EHF KLK11 LAMB3 PKP3 CALML3 ALDH3B2 EXPH5 SFN S100A14 IRX4 ELMO3 FXYD3 GPX2 TRIM29 PDZK1IP1 RAB25 S100A2 TFAP2B AOPEP PRSS8 KRT5 SOX15 KRT15 LY6D C1orf116"                                                                                                                                                                                                                                 
## CHARAFE_BREAST_CANCER_BASAL_VS_MESENCHYMAL_UP "EHF PDZK1IP1 PRSS8 KRT5 C6orf132 C1orf116 PKP3 RAB25 TRIM29 FXYD3 KRT15 PTPN6 S100A14 PROM2 ELMO3 BLNK"                                                                                                                                                                                                                                                                                     
## ONDER_CDH1_TARGETS_2_DN                       "ELMO3 DDR1 CD9 IL1RN GJB5 S100A14 EXPH5 TRIM29 C1orf116 SLC7A5 KLK11 EHF SFN RPS6KA1 S100A2 ARHGAP25 KRT5 SOX15 FXYD3 PKP3 KRT8 KRT15 PRSS8 RAB25 LAMB3 IRX4"                                                                                                                                                                                                                               
## HATADA_METHYLATED_IN_LUNG_CANCER_UP           "HOXA9 NFAM1 REC8 FAM83F HLX HOXD12 PHKG1 PLCB2 OSM CMTM2 LSP1 CD93 RSPH6A NFIX TAGLN CYP2W1 HOXB4 HOXA7 PNMA1 TMEM204 CCDC140 SECTM1"                                                                                                                                                                                                                                                       
## HOX_GENES                                     "HLX IRX4 HOXA2 HOXC10 HOXA11 HOXD12 HOXB3 HOXA10 HOXA7 HOXA9"                                                                                                                                                                                                                                                                                                                               
## CAGGTG_V$E12_Q6                               "OSR1 PKP3 SPI1 ELMO3 CCDC140 PAX3 PAX7 NRG2 CD9 KRT15 FBXO2 CTBP1 SEMA4B ALX3 DDR1 TRPV3 FMNL1 SIM1 LSP1 MFAP4 OTX2 PRDM16 CDH7 HEYL C1orf210 FAM53B HOXA11 HOXA10 SMAD3 SFN ITPK1 HOXA7 HOXB3 SOX15 KHDRBS2 NFIX KRT8 GPT FGF1 LMO2 NR4A2 RAB25 BLNK EXPH5 GPX2 PLA2G3 LRRN4CL AGAP2 LY6G6C TNNI2 FAM83F TNS4 IRX4 SYT8 GCNT3 S100A16 S100A14 ACOT11 EHF CLDN15 WNT3A PPP1R16B SGIP1 PLCB2"
```
We can see that a lung cancer biological term appears at the top, which again is meaningful since lung cancer is strongly related to smoking.

Finally, we can plot the odds ratio and adjPval for each enriched pathway.

```r
plot(x = log2(as.numeric(topGSEAft.lm[[1]][,3])), y = -log10(as.numeric(topGSEAft.lm[[2]][,5])), xlab ="log2(OR)", ylab = "-log10(adjP)", main = "OR and adjP for each enriched pathway", pch = 21, bg = "red")
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)


# Session information


```r
sessionInfo()
```

```
## R version 3.6.3 (2020-02-29)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 18.04.4 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ebGSEA_0.1.0     BiocStyle_2.12.0 rmarkdown_2.2    roxygen2_7.1.0  
## [5] knitr_1.28      
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.4.6         highr_0.8            compiler_3.6.3      
##  [4] BiocManager_1.30.10  bitops_1.0-6         tools_3.6.3         
##  [7] digest_0.6.25        bit_1.1-15.2         lattice_0.20-41     
## [10] memoise_1.1.0        annotate_1.62.0      evaluate_0.14       
## [13] RSQLite_2.1.2        pkgconfig_2.0.3      rlang_0.4.6         
## [16] Matrix_1.2-18        DBI_1.1.0            yaml_2.2.1          
## [19] parallel_3.6.3       xfun_0.14            stringr_1.4.0       
## [22] xml2_1.3.2           IRanges_2.18.3       S4Vectors_0.22.1    
## [25] vctrs_0.3.0          grid_3.6.3           stats4_3.6.3        
## [28] bit64_0.9-7          globaltest_5.38.0    Biobase_2.44.0      
## [31] R6_2.4.1             AnnotationDbi_1.46.1 survival_3.1-12     
## [34] XML_3.99-0.3         org.Hs.eg.db_3.8.2   purrr_0.3.4         
## [37] blob_1.2.1           magrittr_1.5         splines_3.6.3       
## [40] htmltools_0.4.0      BiocGenerics_0.30.0  kpmt_0.1.0          
## [43] xtable_1.8-4         stringi_1.4.3        RCurl_1.98-1.2
```

# References

