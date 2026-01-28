
* **Shaoxun Yuan**

* **Affiliation**: School of Artificial Intelligence and Information Technology, Nanjing University of Chinese Medicine, China

* **Email**: [yuanshaoxun@njucm.edu.cn](mailto:yuanshaoxun@njucm.edu.cn)

# Table of Contents  

1. [Introduction](#1-introduction)  
2. [Installation](#2-installation)  
3. [Load packages required for this tutorial](#3-load-packages-required-for-this-tutorial)  
4. [Input data](#4-input-data)  
5. [Functions](#5-functions)  
   1. [SE_combine](#51-se_combine)  
   2. [SE_combat](#52-se_combat)  
   3. [SE_impute](#53-se_impute)  
   4. [SE_detectratio](#54-se_detectratio)  
   5. [SE_DEseq2](#55-se_deseq2)  
   6. [SE_distribution](#56-se_distribution)  
   7. [SE_boxplot](#57-se_boxplot)  
   8. [SE_PCAplot](#58-se_pcaplot)  
   9. [SE_filter](#59-se_filter)  
   10. [SE_subset](#510-se_subset)  
   11. [SE_normalize](#511-se_normalize)  
   12. [SE_limma](#512-se_limma)  
   13. [SE_edgeR](#513-se-edger)  
   14. [SE_HVG](#514-se-hvg)  
   15. [SE_GSEA](#515-se-gsea)  
   16. [SE_GO](#516-se-go)  
   17. [SE_KEGG](#517-se-kegg)  
   18. [SE_enrichplot](#518-se-enrichplot)  
   19. [SE_hierarchical](#519-se-hierarchical)  
   20. [SE_kmeans](#520-se-kmeans)  
   21. [SE_clusterplot](#521-se-clusterplot)  
   22. [SE_tSNE](#522-se-tsne)  
   23. [SE_UMAP](#523-se-umap)  
   24. [SE_MDS](#524-se-mds)  
   25. [SE_heatmap](#525-se-heatmap)  
   26. [SE_correlation](#526-se-correlation)  
   27. [SE_WGCNA](#527-se-wgcna)  
   28. [SE_networkplot](#528-se-networkplot)  
   29. [SE_randomforest](#529-se-randomforest)  
   30. [SE_SVM](#530-se-svm)  
   31. [SE_crossvalidation](#531-se-crossvalidation)  
   32. [SE_AUCcalc](#532-se-auccalc)  
   33. [SE_sampleQC](#533-se-sampleqc)  
   34. [SE_featureQC](#534-se-featureqc)  
   35. [SE_QCreport](#535-se-qcreport)  
   36. [SE_export](#536-se-export)  
   37. [SE_metadata](#537-se-metadata)  
   38. [SE_rename](#538-se-rename)  
   39. [SE_convert](#539-se-convert)  
   40. [SE_merge](#540-se-merge)  
   41. [SE_summary](#541-se-summary)  
   42. [SE_volcano](#542-se-volcano)  
   43. [SE_MAplot](#543-se-maplot)  
   44. [SE_density](#544-se-density)  
   45. [SE_trend](#545-se-trend)  
   46. [SE_timeseries](#546-se-timeseries)  
   47. [SE_transform](#547-se-transform)  
   48. [SE_batchdetect](#548-se-batchdetect)  
   49. [SE_silhouette](#549-se-silhouette)  
   50. [SE_featureselection](#550-se-featureselection)  
   51. [SE_circTest](#551-se-circtest)  
   52. [SE_COCONUT](#552-se-coconut)  
   53. [loadSE](#553-loadse)  
   54. [loadSElist](#554-loadselist)  
6. [Acknowledgements](#6-acknowledgements)  
7. [Session Info](#7-session-info)  
8. [References](#8-references)

## 1. Introduction

SEtoolbox is an R package that operates, analyzes and visualizes SummarizedExperiment objects.

## 2. Installation

To install the SEtoolbox package, you first need to install the `devtools` package. 

Run the following command in your R console: 
```r
install.packages("devtools")

devtools::install_github("shaoxunyuan/SEtoolbox")
```

## 3. Load packages required for this tutorial

During this tutorial, we might need to use a few additional packages.

Since we specified dependencies = TRUE when installing SEtoolbox package, these additional packages have already been installed.

We can load them directly. 

```r
library(SummarizedExperiment)

library(tidyverse)

library(plyr)

library(dplyr)

library(reshape2)

library(DESeq2)
```

## 4. Input data

For this tutorial, SEtoolbox will be working with a [`SummarizedExperiment`](https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html) object.

## 5. Functions

Functions in SEtoolbox can be obtain using 

```r
help(package="SEtoolbox")
```

### 5.1. SE_combine

Combine multiple `SummarizedExperiment` objects.

This function merges multiple `SummarizedExperiment` objects based on the specified merge type (intersection or union) for all assays present in the input list.

#### Usage

    SE_combine(se_list, merge_type = "intersection") 

#### Arguments

`se_list`: A list of SummarizedExperiment objects to be combined (se_list = list(SE1,SE2,SE3))

`merge_type`: A character string specifying the type of merge to perform. Options are:

    1.`intersection` (default): Keep only common features across all objects.

    2.`union`: Keep all features.

### 5.2. SE_combat

Batch effect correct using Combat in sva package for `SummarizedExperiment` object.

#### Usage

    SE_combat(SE, col_for_combat, col_for_compare) 

#### Arguments

	`SE`: A `SummarizedExperiment` object containing the data to be imputed.

### 5.3 SE_impute

Fill missing values in a `SummarizedExperiment` object.

This function imputes missing values (NA) in the given `SummarizedExperiment` object using specified methods.

Multiple imputation techniques can be utilized to handle missing values, ensuring the robustness of subsequent analyses.

#### Usage

    SE_impute(SE, assayname = "TPM", group = "group", ZerosAsNA = FALSE, RemoveNA = TRUE,  
              cutoff = 20, method = c("none", "LOD", "half_min", "median", "mean", "min", "knn", "rf", "global_mean", "svd", "QRILC"),  
              LOD = NULL, knum = 10)  

#### Arguments

`SE`: A `SummarizedExperiment` object containing the data to be imputed.

`assayname`: The name of the `SummarizedExperiment` assay, specifying the type of data to be imputed.

`group`: A character string specifying the grouping variable in the sample data.

`ZerosAsN`: A logical value indicating whether to treat zeros as NA. Default is `FALSE`.

`RemoveNA`: A logical value indicating whether to remove samples with a high percentage of NA values based on the cutoff. Default is `TRUE`.

`cutoff`: A numerical value representing the percentage cutoff for NA samples. Default is `20`.

`method`: A character string specifying the imputation method to use. Options include:

    1.`none` (default): No imputation, replace NA with zero.

    2.`LOD`: Replace NA with the limit of detection (LOD).

    3.`half_min`: Replace NA with half of the minimum value.

    4.`median`: Replace NA with the median value.

    5.`mean`: Replace NA with the mean value.

    6.`min`: Replace NA with the minimum value.

    7.`knn`: K-nearest neighbors imputation.

    8.`rf`: Random forest imputation.

    9.`global_mean`: Global mean imputation.

    10.`svd`: Singular value decomposition imputation.

    11.`QRILC`: Quantile regression imputation.

`LOD`: A numerical value representing the limit of detection (used for the LOD imputation method). Default is `NULL`.

`knum`: An integer value representing the number of neighbors in the KNN imputation method. Default is `10`.

### 5.4 SE_detectraio

Calculate detection ratio and update `SummarizedExperiment` object's `rowData`.

This function computes the detection ratio of expression data and updates the `rowData` of the provided `SummarizedExperiment` object with detection sample counts and ratios. It also generates a histogram of detection ratios.

#### Usage

    SE_detectratio(SE, assayname = "TPM")  

#### Arguments

`SE`: A `SummarizedExperiment` object containing expression data.

`assayname`: The name of the assay to be used for calculations. Default is `"TPM"`.

### 5.5 SE_DEseq2

Perform differential expression analysis using `DESeq2`.

This function performs differential expression analysis on count data contained in a `SummarizedExperiment` object using the `DESeq2` package.

#### Usage

    SE_DEseq2(SE, assayname = "Count", groupname = "group")  

#### Arguments

`SE`: An `SummarizedExperiment` object containing count data.

`assayname`: The name of the assay to use for the analysis. Default is `"Count"`.

`groupname`: The name of the column in `colData(SE)` that contains the factor for grouping samples. Default is `"group"`.

### 5.6 SE_distribution  

Generates two plots illustrating the distribution of non-zero entries in a `SummarizedExperiment` object.  

This function creates a bar plot to display the count of non-zero entries for each feature (gene) in the specified expression matrix and a histogram showing the distribution of the fraction of non-zero entries across samples.  

It allows for the option to treat zeros as `NA`, thereby excluding them from the count.  

#### Usage  

SE_distribution(SE, assayname = "TPM", ZeroasNA = TRUE)

#### Arguments

`SE`: An `SummarizedExperiment` object containing the expression data.

`assayname`: A string indicating the assay name to use from the SummarizedExperiment object. The default is "TPM".

`ZeroasNA`: A logical value indicating whether zeros should be treated as NA. The default is TRUE.

### 5.7 SE_boxplot

Generates a boxplot for specified features in a `SummarizedExperiment` object.

This function creates a boxplot for the specified features (genes) within a given `SummarizedExperiment` object.

It supports normalization and grouping of the data.

#### Usage

    SE_boxplot(SE, feature_of_interest, assayname = "TPM", group_col = NA, normalization = "none")  

#### Arguments

`SE`: A SummarizedExperiment object.

`feature_of_interest`: Character vector of gene identifiers.

`assayname`: The assay name in the SummarizedExperiment object.

`group_col`: Column name in colData(SE) for grouping.

`normalization`: Normalization method ("scale", "log", or "none").

### 5.8 SE_PCAplot

Generate PCA plots.

This function takes a `SummarizedExperiment` object, computes PCA, and visualizes the results.

#### Usage

    SE_PCAplot(SE, assayname = "TPM", groupname = "group", outlier_threshold = 2, scale = TRUE, feature_of_interesting = NULL)  

#### Arguments

`SE`: SummarizedExperiment object containing gene expression data.

`assayname`: Name of the expression data, default is "TPM".

`groupname`: Name of the grouping column, default is "group".

`outlier_threshold`: Outlier filtering threshold, default is 2.

`scale`: Whether to standardize the data, default is TRUE.

`feature_of_interesting`: Vector of specific feature names; if NULL, all features are used, default is NULL.

### 5.9 SE_filter

Filter features in SummarizedExperiment object based on various criteria.

This function filters features (genes/proteins) in a SummarizedExperiment object based on various criteria including expression level, detection ratio, variance, and coefficient of variation. It provides flexible filtering options to remove low-quality features before downstream analysis.

#### Usage

    SE_filter(SE, assayname = "TPM", min_expr = 0, min_detectratio = 0, min_variance = 0, min_cv = 0, group_colname = NULL)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use for filtering. The default value is "TPM".

`min_expr`: Numeric value. Minimum expression threshold. Features with mean expression below this value will be removed. Default is 0.

`min_detectratio`: Numeric value between 0 and 1. Minimum detection ratio threshold. Features with detection ratio below this value will be removed. Default is 0.

`min_variance`: Numeric value. Minimum variance threshold. Features with variance below this value will be removed. Default is 0.

`min_cv`: Numeric value. Minimum coefficient of variation (CV) threshold. Features with CV below this value will be removed. Default is 0.

`group_colname`: A string representing the column name in `colData` that contains group information. This is optional; if provided, filtering will be performed within each group separately.

### 5.10 SE_subset

Subset SummarizedExperiment object by features or samples.

This function subsets a SummarizedExperiment object based on feature names, sample names, or conditions from colData. It provides flexible subsetting options to extract specific subsets of data for downstream analysis.

#### Usage

    SE_subset(SE, features = NULL, samples = NULL, condition = NULL, exclude_features = NULL, exclude_samples = NULL)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`features`: A character vector of feature names to keep. If NULL, all features are kept. Default is NULL.

`samples`: A character vector of sample names to keep. If NULL, all samples are kept. Default is NULL.

`condition`: A named list or vector for conditional subsetting. Names should be column names in colData, and values should be the values to keep. Default is NULL.

`exclude_features`: A character vector of feature names to exclude. Default is NULL.

`exclude_samples`: A character vector of sample names to exclude. Default is NULL.

### 5.11 SE_normalize

Normalize expression data in SummarizedExperiment object.

This function normalizes expression data in a SummarizedExperiment object using various normalization methods including TPM, FPKM, RPKM, log2 transformation, quantile normalization, and library size normalization.

#### Usage

    SE_normalize(SE, assayname = "Counts", method = "log2", pseudocount = 1, gene_length = NULL)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to normalize. The default value is "Counts".

`method`: A character string specifying the normalization method. Options include "TPM", "FPKM", "RPKM", "log2", "quantile", "library_size", "median", "upper_quartile". Default is "log2".

`pseudocount`: Numeric value to add before log transformation to avoid log(0). Default is 1.

`gene_length`: A numeric vector of gene lengths in base pairs. Required for TPM/FPKM/RPKM methods.

### 5.12 SE_limma

Perform differential expression analysis using limma.

This function performs differential expression analysis using the limma package, which is suitable for both microarray and RNA-seq data. It uses linear models and empirical Bayes methods to identify differentially expressed features.

#### Usage

    SE_limma(SE, assayname = "log2", group_colname = "group", contrast = NULL, adjust_method = "BH", pvalue_threshold = 0.05, logFC_threshold = 1)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use for analysis. The default value is "log2".

`group_colname`: A string representing the column name in `colData` that contains group information. Default is "group".

`contrast`: A character string specifying the contrast for differential analysis (e.g., "Treatment-Control"). Default is NULL, which will compare the first two groups.

`adjust_method`: A character string specifying the p-value adjustment method. Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default is "BH".

`pvalue_threshold`: Numeric value for p-value threshold. Default is 0.05.

`logFC_threshold`: Numeric value for log2 fold change threshold. Default is 1.

### 5.13 SE_edgeR

Perform differential expression analysis using edgeR.

This function performs differential expression analysis on count data using the edgeR package, which is specifically designed for RNA-seq data analysis.

#### Usage

    SE_edgeR(SE, assayname = "Counts", group_colname = "group", contrast = NULL, adjust_method = "BH", pvalue_threshold = 0.05, logFC_threshold = 1)

#### Arguments

`SE`: A `SummarizedExperiment` object containing count data.

`assayname`: A string indicating which assay to use for analysis. The default value is "Counts".

`group_colname`: A string representing the column name in `colData` that contains group information. Default is "group".

`contrast`: A character string specifying the contrast for differential analysis (e.g., "Treatment-Control"). Default is NULL, which will compare the first two groups.

`adjust_method`: A character string specifying the p-value adjustment method. Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default is "BH".

`pvalue_threshold`: Numeric value for p-value threshold. Default is 0.05.

`logFC_threshold`: Numeric value for log2 fold change threshold. Default is 1.

### 5.14 SE_HVG

Identify highly variable genes in SummarizedExperiment object.

This function identifies highly variable genes (HVGs) in a SummarizedExperiment object using various methods including variance, coefficient of variation (CV), and the method from the scran package.

#### Usage

    SE_HVG(SE, assayname = "TPM", method = "variance", n_top_genes = 1000, min_mean = 0.1, max_mean = 1000, min_dispersion = 0.1)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`method`: A character string specifying the method to use for identifying HVGs. Options include "variance", "cv", "scran". Default is "variance".

`n_top_genes`: Numeric value. Number of top highly variable genes to return. Default is 1000.

`min_mean`: Numeric value. Minimum mean expression threshold for scran method. Default is 0.1.

`max_mean`: Numeric value. Maximum mean expression threshold for scran method. Default is 1000.

`min_dispersion`: Numeric value. Minimum dispersion threshold for scran method. Default is 0.1.

### 5.15 SE_GSEA

Perform Gene Set Enrichment Analysis (GSEA).

This function performs Gene Set Enrichment Analysis (GSEA) on a SummarizedExperiment object using the fgsea package. It identifies gene sets that are significantly enriched in a ranked list of genes.

#### Usage

    SE_GSEA(SE, assayname = "log2", group_colname = "group", gene_sets, min_size = 10, max_size = 500, n_perm = 1000, adjust_method = "BH", pvalue_threshold = 0.05)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "log2".

`group_colname`: A string representing the column name in `colData` that contains group information. Default is "group".

`gene_sets`: A list of gene sets where each element is a character vector of gene symbols.

`min_size`: Numeric value. Minimum size of gene sets to consider. Default is 10.

`max_size`: Numeric value. Maximum size of gene sets to consider. Default is 500.

`n_perm`: Numeric value. Number of permutations to perform. Default is 1000.

`adjust_method`: A character string specifying the p-value adjustment method. Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default is "BH".

`pvalue_threshold`: Numeric value for p-value threshold. Default is 0.05.

### 5.16 SE_GO

Perform Gene Ontology (GO) enrichment analysis.

This function performs Gene Ontology (GO) enrichment analysis on a list of genes using the clusterProfiler package. It identifies GO terms that are significantly enriched in the input gene list.

#### Usage

    SE_GO(SE, genes, ont = "BP", organism = "human", pvalue_cutoff = 0.05, qvalue_cutoff = 0.2, min_gene = 10, max_gene = 500)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`genes`: A character vector of gene symbols.

`ont`: A character string specifying the GO ontology to use. Options include "BP" (Biological Process), "CC" (Cellular Component), "MF" (Molecular Function). Default is "BP".

`organism`: A character string specifying the organism. Default is "human".

`pvalue_cutoff`: Numeric value for p-value cutoff. Default is 0.05.

`qvalue_cutoff`: Numeric value for q-value cutoff. Default is 0.2.

`min_gene`: Numeric value. Minimum number of genes in a GO term. Default is 10.

`max_gene`: Numeric value. Maximum number of genes in a GO term. Default is 500.

### 5.17 SE_KEGG

Perform KEGG pathway enrichment analysis.

This function performs KEGG pathway enrichment analysis on a list of genes using the clusterProfiler package. It identifies KEGG pathways that are significantly enriched in the input gene list.

#### Usage

    SE_KEGG(SE, genes, organism = "hsa", pvalue_cutoff = 0.05, qvalue_cutoff = 0.2, min_gene = 10, max_gene = 500)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`genes`: A character vector of gene symbols.

`organism`: A character string specifying the KEGG organism code. Default is "hsa" (human).

`pvalue_cutoff`: Numeric value for p-value cutoff. Default is 0.05.

`qvalue_cutoff`: Numeric value for q-value cutoff. Default is 0.2.

`min_gene`: Numeric value. Minimum number of genes in a KEGG pathway. Default is 10.

`max_gene`: Numeric value. Maximum number of genes in a KEGG pathway. Default is 500.

### 5.18 SE_enrichplot

Visualize enrichment analysis results.

This function visualizes the results of enrichment analysis (GO/KEGG/GSEA) using various plot types including bar plots, dot plots, and network plots.

#### Usage

    SE_enrichplot(enrich_result, plot_type = "bar", n_terms = 10, pvalue_cutoff = 0.05)

#### Arguments

`enrich_result`: An enrichment analysis result object from clusterProfiler or fgsea.

`plot_type`: A character string specifying the plot type. Options include "bar", "dot", "network", "cnet". Default is "bar".

`n_terms`: Numeric value. Number of top terms to display. Default is 10.

`pvalue_cutoff`: Numeric value for p-value cutoff. Default is 0.05.

### 5.19 SE_hierarchical

Perform hierarchical clustering analysis.

This function performs hierarchical clustering on a SummarizedExperiment object and visualizes the results as a dendrogram. It can cluster either samples or features.

#### Usage

    SE_hierarchical(SE, assayname = "TPM", cluster_by = "samples", method = "complete", distance = "euclidean")

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`cluster_by`: A character string specifying whether to cluster samples or features. Options include "samples", "features". Default is "samples".

`method`: A character string specifying the clustering method. Options include "complete", "single", "average", "ward.D", "ward.D2". Default is "complete".

`distance`: A character string specifying the distance metric. Options include "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski". Default is "euclidean".

### 5.20 SE_kmeans

Perform k-means clustering analysis.

This function performs k-means clustering on a SummarizedExperiment object and adds cluster assignments to the rowData or colData.

#### Usage

    SE_kmeans(SE, assayname = "TPM", n_clusters = 3, cluster_by = "samples", nstart = 10)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`n_clusters`: Numeric value. Number of clusters to form. Default is 3.

`cluster_by`: A character string specifying whether to cluster samples or features. Options include "samples", "features". Default is "samples".

`nstart`: Numeric value. Number of random initializations. Default is 10.

### 5.21 SE_clusterplot

Visualize clustering results.

This function visualizes clustering results using dimensionality reduction techniques such as PCA, t-SNE, or UMAP, with points colored by cluster assignment.

#### Usage

    SE_clusterplot(SE, assayname = "TPM", cluster_colname, dim_reduction = "PCA", group_colname = NULL)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`cluster_colname`: A string representing the column name in rowData or colData that contains cluster assignments.

`dim_reduction`: A character string specifying the dimensionality reduction method. Options include "PCA", "tSNE", "UMAP". Default is "PCA".

`group_colname`: A string representing the column name in colData that contains group information. Default is NULL.

### 5.22 SE_tSNE

Perform t-distributed Stochastic Neighbor Embedding (t-SNE) analysis.

This function performs t-SNE dimensionality reduction on a SummarizedExperiment object and visualizes the results.

#### Usage

    SE_tSNE(SE, assayname = "TPM", group_colname = "group", perplexity = 30, theta = 0.5, max_iter = 1000)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`group_colname`: A string representing the column name in `colData` that contains group information. Default is "group".

`perplexity`: Numeric value. Perplexity parameter for t-SNE. Default is 30.

`theta`: Numeric value. Speed/accuracy trade-off parameter for t-SNE. Default is 0.5.

`max_iter`: Numeric value. Maximum number of iterations for t-SNE. Default is 1000.

### 5.23 SE_UMAP

Perform Uniform Manifold Approximation and Projection (UMAP) analysis.

This function performs UMAP dimensionality reduction on a SummarizedExperiment object and visualizes the results.

#### Usage

    SE_UMAP(SE, assayname = "TPM", group_colname = "group", n_neighbors = 15, min_dist = 0.1, n_components = 2)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`group_colname`: A string representing the column name in `colData` that contains group information. Default is "group".

`n_neighbors`: Numeric value. Number of neighbors for UMAP. Default is 15.

`min_dist`: Numeric value. Minimum distance for UMAP. Default is 0.1.

`n_components`: Numeric value. Number of components for UMAP. Default is 2.

### 5.24 SE_MDS

Perform Multidimensional Scaling (MDS) analysis.

This function performs MDS dimensionality reduction on a SummarizedExperiment object and visualizes the results.

#### Usage

    SE_MDS(SE, assayname = "TPM", group_colname = "group", distance = "euclidean")

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`group_colname`: A string representing the column name in `colData` that contains group information. Default is "group".

`distance`: A character string specifying the distance metric. Options include "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski". Default is "euclidean".

### 5.25 SE_heatmap

Generate heatmap for gene expression data.

This function generates a heatmap for gene expression data in a SummarizedExperiment object, with optional clustering and annotation.

#### Usage

    SE_heatmap(SE, assayname = "TPM", features = NULL, group_colname = "group", scale = "row", clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete")

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`features`: A character vector of feature names to include in the heatmap. If NULL, all features are included. Default is NULL.

`group_colname`: A string representing the column name in `colData` that contains group information. Default is "group".

`scale`: A character string specifying the scaling method. Options include "row", "column", "none". Default is "row".

`clustering_distance_rows`: A character string specifying the distance metric for row clustering. Default is "euclidean".

`clustering_distance_cols`: A character string specifying the distance metric for column clustering. Default is "euclidean".

`clustering_method`: A character string specifying the clustering method. Default is "complete".

### 5.26 SE_correlation

Calculate correlation between features or samples.

This function calculates correlation coefficients between features or samples in a SummarizedExperiment object and visualizes the results as a heatmap.

#### Usage

    SE_correlation(SE, assayname = "TPM", correlation_type = "pearson", correlate_by = "features", features = NULL)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`correlation_type`: A character string specifying the correlation method. Options include "pearson", "spearman", "kendall". Default is "pearson".

`correlate_by`: A character string specifying whether to correlate features or samples. Options include "features", "samples". Default is "features".

`features`: A character vector of feature names to include. If NULL, all features are included. Default is NULL.

### 5.27 SE_WGCNA

Perform Weighted Gene Co-expression Network Analysis (WGCNA).

This function performs WGCNA on a SummarizedExperiment object to identify co-expression modules and their relationships with sample traits.

#### Usage

    SE_WGCNA(SE, assayname = "TPM", power = 6, minModuleSize = 30, mergeCutHeight = 0.25, traitData = NULL)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`power`: Numeric value. Soft-thresholding power for network construction. Default is 6.

`minModuleSize`: Numeric value. Minimum module size. Default is 30.

`mergeCutHeight`: Numeric value. Cut height for merging modules. Default is 0.25.

`traitData`: A data frame containing sample traits. If NULL, colData(SE) is used. Default is NULL.

### 5.28 SE_networkplot

Visualize gene co-expression networks.

This function visualizes gene co-expression networks or protein-protein interaction networks using igraph and ggraph.

#### Usage

    SE_networkplot(network, layout = "fr", n_nodes = 100, n_edges = 500, node_size = 5, node_color = "blue", edge_color = "gray")

#### Arguments

`network`: A network object (igraph or adjacency matrix).

`layout`: A character string specifying the layout algorithm. Options include "fr", "kk", "circle", "grid". Default is "fr".

`n_nodes`: Numeric value. Number of top nodes to include. Default is 100.

`n_edges`: Numeric value. Number of top edges to include. Default is 500.

`node_size`: Numeric value. Node size. Default is 5.

`node_color`: Character string. Node color. Default is "blue".

`edge_color`: Character string. Edge color. Default is "gray".

### 5.29 SE_randomforest

Perform random forest classification or regression.

This function performs random forest analysis on a SummarizedExperiment object for classification or regression tasks.

#### Usage

    SE_randomforest(SE, assayname = "TPM", group_colname = "group", n_trees = 500, mtry = NULL, importance = TRUE, ntree = 500)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`group_colname`: A string representing the column name in `colData` that contains the response variable. Default is "group".

`n_trees`: Numeric value. Number of trees in the forest. Default is 500.

`mtry`: Numeric value. Number of variables randomly sampled as candidates at each split. If NULL, default value is used. Default is NULL.

`importance`: Logical value. Whether to calculate variable importance. Default is TRUE.

`ntree`: Numeric value. Number of trees in the forest (alternative parameter name). Default is 500.

### 5.30 SE_SVM

Perform Support Vector Machine (SVM) classification or regression.

This function performs SVM analysis on a SummarizedExperiment object for classification or regression tasks.

#### Usage

    SE_SVM(SE, assayname = "TPM", group_colname = "group", kernel = "radial", cost = 1, gamma = "auto")

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`group_colname`: A string representing the column name in `colData` that contains the response variable. Default is "group".

`kernel`: A character string specifying the kernel type. Options include "linear", "polynomial", "radial", "sigmoid". Default is "radial".

`cost`: Numeric value. Cost parameter. Default is 1.

`gamma`: Numeric value or "auto". Gamma parameter. Default is "auto".

### 5.31 SE_crossvalidation

Perform cross-validation for machine learning models.

This function performs cross-validation for machine learning models (random forest, SVM) to evaluate model performance.

#### Usage

    SE_crossvalidation(SE, assayname = "TPM", group_colname = "group", model_type = "randomforest", k_folds = 5, n_repeats = 1)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`group_colname`: A string representing the column name in `colData` that contains the response variable. Default is "group".

`model_type`: A character string specifying the model type. Options include "randomforest", "SVM". Default is "randomforest".

`k_folds`: Numeric value. Number of folds for cross-validation. Default is 5.

`n_repeats`: Numeric value. Number of repeats for cross-validation. Default is 1.

### 5.32 SE_AUCcalc

Calculate Area Under the Curve (AUC) for classification models.

This function calculates the AUC for classification models using the pROC package.

#### Usage

    SE_AUCcalc(predictions, labels)

#### Arguments

`predictions`: A numeric vector of predicted values or probabilities.

`labels`: A factor vector of true class labels.

### 5.33 SE_sampleQC

Perform sample quality control.

This function performs quality control on samples in a SummarizedExperiment object, including assessment of library size, number of detected features, and sample correlation.

#### Usage

    SE_sampleQC(SE, assayname = "TPM", group_colname = "group")

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`group_colname`: A string representing the column name in `colData` that contains group information. Default is "group".

### 5.34 SE_featureQC

Perform feature quality control.

This function performs quality control on features in a SummarizedExperiment object, including assessment of expression levels, detection rates, and variability.

#### Usage

    SE_featureQC(SE, assayname = "TPM", group_colname = "group")

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`group_colname`: A string representing the column name in `colData` that contains group information. Default is "group".

### 5.35 SE_QCreport

Generate comprehensive quality control report.

This function generates a comprehensive quality control report for a SummarizedExperiment object, including sample and feature QC metrics, and data visualization.

#### Usage

    SE_QCreport(SE, assayname = "TPM", group_colname = "group", output_dir = ".")

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`group_colname`: A string representing the column name in `colData` that contains group information. Default is "group".

`output_dir`: A string specifying the output directory for the report. Default is ".".

### 5.36 SE_export

Export SummarizedExperiment object to various formats.

This function exports a SummarizedExperiment object to various formats including CSV, TSV, Excel, and RData.

#### Usage

    SE_export(SE, assayname = "TPM", output_format = "csv", output_file = "SE_export")

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to export. The default value is "TPM".

`output_format`: A character string specifying the output format. Options include "csv", "tsv", "xlsx", "rdata". Default is "csv".

`output_file`: A string specifying the output file name (without extension). Default is "SE_export".

### 5.37 SE_metadata

Manage metadata for SummarizedExperiment object.

This function manages metadata for a SummarizedExperiment object, including adding, removing, and updating metadata entries.

#### Usage

    SE_metadata(SE, metadata = NULL, action = "add")

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`metadata`: A named list of metadata entries. Default is NULL.

`action`: A character string specifying the action to perform. Options include "add", "remove", "update", "get". Default is "add".

### 5.38 SE_rename

Rename features or samples in SummarizedExperiment object.

This function renames features or samples in a SummarizedExperiment object based on a mapping provided by the user.

#### Usage

    SE_rename(SE, mapping, rename_by = "features")

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`mapping`: A named vector or data frame where names are old identifiers and values are new identifiers.

`rename_by`: A character string specifying whether to rename features or samples. Options include "features", "samples". Default is "features".

### 5.39 SE_convert

Convert SummarizedExperiment object to other formats.

This function converts a SummarizedExperiment object to other formats including ExpressionSet, DGEList, and data frames.

#### Usage

    SE_convert(SE, assayname = "TPM", to_format = "ExpressionSet")

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use for conversion. The default value is "TPM".

`to_format`: A character string specifying the format to convert to. Options include "ExpressionSet", "DGEList", "data.frame", "matrix". Default is "ExpressionSet".

### 5.40 SE_merge

Merge multiple SummarizedExperiment objects.

This function merges multiple SummarizedExperiment objects based on common features or samples.

#### Usage

    SE_merge(se_list, merge_by = "features", merge_type = "inner")

#### Arguments

`se_list`: A list of SummarizedExperiment objects to merge.

`merge_by`: A character string specifying whether to merge by features or samples. Options include "features", "samples". Default is "features".

`merge_type`: A character string specifying the merge type. Options include "inner", "outer", "left", "right". Default is "inner".

### 5.41 SE_summary

Generate summary statistics for SummarizedExperiment object.

This function generates summary statistics for a SummarizedExperiment object, including number of features and samples, expression statistics, and metadata information.

#### Usage

    SE_summary(SE, assayname = "TPM")

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use for summary statistics. The default value is "TPM".

### 5.42 SE_volcano

Generate volcano plot for differential expression analysis results.

This function generates a volcano plot for differential expression analysis results, showing the relationship between log2 fold change and adjusted p-value.

#### Usage

    SE_volcano(SE, assayname = "log2", logFC_threshold = 1, pvalue_threshold = 0.05)

#### Arguments

`SE`: A `SummarizedExperiment` object containing differential expression results in rowData.

`assayname`: A string indicating which assay to use. The default value is "log2".

`logFC_threshold`: Numeric value for log2 fold change threshold. Default is 1.

`pvalue_threshold`: Numeric value for adjusted p-value threshold. Default is 0.05.

### 5.43 SE_MAplot

Generate MA plot for differential expression analysis results.

This function generates an MA plot for differential expression analysis results, showing the relationship between mean expression and log2 fold change.

#### Usage

    SE_MAplot(SE, assayname = "log2", logFC_threshold = 1, pvalue_threshold = 0.05)

#### Arguments

`SE`: A `SummarizedExperiment` object containing differential expression results in rowData.

`assayname`: A string indicating which assay to use. The default value is "log2".

`logFC_threshold`: Numeric value for log2 fold change threshold. Default is 1.

`pvalue_threshold`: Numeric value for adjusted p-value threshold. Default is 0.05.

### 5.44 SE_density

Generate density plot for gene expression data.

This function generates a density plot for gene expression data in a SummarizedExperiment object, showing the distribution of expression values.

#### Usage

    SE_density(SE, assayname = "TPM", group_colname = "group")

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`group_colname`: A string representing the column name in `colData` that contains group information. Default is "group".

### 5.45 SE_trend

Analyze expression trends across conditions or time points.

This function analyzes expression trends across conditions or time points in a SummarizedExperiment object.

#### Usage

    SE_trend(SE, assayname = "TPM", group_colname = "group", feature_of_interest)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`group_colname`: A string representing the column name in `colData` that contains group information. Default is "group".

`feature_of_interest`: A character vector of feature names to analyze.

### 5.46 SE_timeseries

Analyze time series gene expression data.

This function analyzes time series gene expression data in a SummarizedExperiment object, including identifying temporal patterns and differentially expressed genes across time points.

#### Usage

    SE_timeseries(SE, assayname = "TPM", time_colname = "time", group_colname = "group")

#### Arguments

`SE`: A `SummarizedExperiment` object containing time series gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`time_colname`: A string representing the column name in `colData` that contains time point information. Default is "time".

`group_colname`: A string representing the column name in `colData` that contains group information. Default is "group".

### 5.47 SE_transform

Transform expression data in SummarizedExperiment object.

This function transforms expression data in a SummarizedExperiment object using various transformation methods including log transformation, z-score normalization, and quantile transformation.

#### Usage

    SE_transform(SE, assayname = "TPM", method = "log2", pseudocount = 1)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to transform. The default value is "TPM".

`method`: A character string specifying the transformation method. Options include "log2", "log10", "zscore", "quantile". Default is "log2".

`pseudocount`: Numeric value to add before log transformation to avoid log(0). Default is 1.

### 5.48 SE_batchdetect

Detect batch effects in SummarizedExperiment object.

This function detects batch effects in a SummarizedExperiment object using principal component analysis (PCA) and correlation analysis.

#### Usage

    SE_batchdetect(SE, assayname = "TPM", batch_colname = "batch")

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`batch_colname`: A string representing the column name in `colData` that contains batch information. Default is "batch".

### 5.49 SE_silhouette

Calculate silhouette width for clustering evaluation.

This function calculates silhouette width for clustering results to evaluate cluster quality.

#### Usage

    SE_silhouette(SE, assayname = "TPM", cluster_colname, distance = "euclidean")

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`cluster_colname`: A string representing the column name in rowData or colData that contains cluster assignments.

`distance`: A character string specifying the distance metric. Options include "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski". Default is "euclidean".

### 5.50 SE_featureselection

Perform feature selection for machine learning.

This function performs feature selection for machine learning tasks using various methods including variance threshold, recursive feature elimination, and LASSO.

#### Usage

    SE_featureselection(SE, assayname = "TPM", group_colname = "group", method = "variance", n_features = 100)

#### Arguments

`SE`: A `SummarizedExperiment` object containing gene expression data.

`assayname`: A string indicating which assay to use. The default value is "TPM".

`group_colname`: A string representing the column name in `colData` that contains the response variable. Default is "group".

`method`: A character string specifying the feature selection method. Options include "variance", "rfe", "lasso". Default is "variance".

`n_features`: Numeric value. Number of features to select. Default is 100.

### 5.51 SE_circTest

Perform circular permutation test for enrichment analysis.

This function performs circular permutation test for enrichment analysis, which is useful for testing enrichment of gene sets in ranked lists.

#### Usage

    SE_circTest(ranked_list, gene_set, n_perm = 1000)

#### Arguments

`ranked_list`: A named numeric vector of ranked values.

`gene_set`: A character vector of gene symbols.

`n_perm`: Numeric value. Number of permutations to perform. Default is 1000.

### 5.52 SE_COCONUT

Perform COCONUT analysis for functional annotation.

This function performs COCONUT (COmprehensive COnsortium Network Utility Tool) analysis for functional annotation of genes.

#### Usage

    SE_COCONUT(genes, organism = "human", database = "GO")

#### Arguments

`genes`: A character vector of gene symbols.

`organism`: A character string specifying the organism. Default is "human".

`database`: A character string specifying the database to use. Options include "GO", "KEGG", "Reactome". Default is "GO".

### 5.53 loadSE

Load example SummarizedExperiment object.

This function loads an example SummarizedExperiment object for testing and demonstration purposes.

#### Usage

    loadSE()

#### Arguments

None.

### 5.54 loadSElist

Load list of example SummarizedExperiment objects.

This function loads a list of example SummarizedExperiment objects for testing and demonstration purposes.

#### Usage

    loadSElist()

#### Arguments

None.

## 6. Acknowledgements

The author would like to thank Deepseek.

## 7. Session Info

    > sessionInfo()
    R version 4.2.0 (2022-04-22)
    Platform: x86_64-conda-linux-gnu (64-bit)
    Running under: CentOS Linux 7 (Core)
    
    Matrix products: default
    BLAS/LAPACK: /home/shaoxun/anaconda3/envs/yuanshaoxun/lib/libopenblasp-r0.3.21.so
    
    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
     [9] LC_ADDRESS=C               LC_TELEPHONE=C
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
    
    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base
    
    other attached packages:
    [1] SEtoolbox_0.1.0
    
    loaded via a namespace (and not attached):
     [1] Rcpp_1.0.14                 pillar_1.10.1
     [3] bslib_0.9.0                 compiler_4.2.0
     [5] later_1.4.1                 jquerylib_0.1.4
     [7] GenomeInfoDb_1.34.9         XVector_0.38.0
     [9] MatrixGenerics_1.10.0       bitops_1.0-9
    [11] tools_4.2.0                 zlibbioc_1.44.0
    [13] digest_0.6.37               tibble_3.2.1
    [15] gtable_0.3.6                lattice_0.22-6
    [17] jsonlite_1.9.1              memoise_2.0.1
    [19] lifecycle_1.0.4             pkgconfig_2.0.3
    [21] rlang_1.1.5                 Matrix_1.6-5
    [23] DelayedArray_0.24.0         shiny_1.10.0
    [25] cli_3.6.4.9000              fastmap_1.2.0
    [27] GenomeInfoDbData_1.2.9      dplyr_1.1.4
    [29] generics_0.1.3              vctrs_0.6.5
    [31] S4Vectors_0.36.2            sass_0.4.9
    [33] IRanges_2.32.0              tidyselect_1.2.1
    [35] grid_4.2.0                  stats4_4.2.0
    [37] glue_1.8.0                  Biobase_2.58.0
    [39] R6_2.6.1                    ggplot2_3.5.1
    [41] magrittr_2.0.3              scales_1.3.0
    [43] promises_1.3.2              matrixStats_1.5.0
    [45] htmltools_0.5.8.1           BiocGenerics_0.44.0
    [47] GenomicRanges_1.50.2        SummarizedExperiment_1.28.0
    [49] colorspace_2.1-1            mime_0.12
    [51] xtable_1.8-4                httpuv_1.6.15
    [53] munsell_0.5.1               RCurl_1.98-1.16
    [55] cachem_1.1.0

## 8. References

Will update soon

Last updated: 2026-01-28
