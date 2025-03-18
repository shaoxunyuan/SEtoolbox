
* **Shaoxun Yuan**

* **Affiliation**: School of Artificial Intelligence and Information Technology, Nanjing University of Chinese Medicine, China

* **Email**: [yuanshaoxun@njucm.edu.cn](mailto:yuanshaoxun@njucm.edu.cn)

{% include toc.html %}

# Table of Contents  

1. [Introduction](#1-introduction)  
2. [Installation](#2-installation)  
3. [Load Packages Required for This Tutorial](#3-load-packages-required-for-this-tutorial)  
4. [Input Data](#4-input-data)  
5. [Functions](#5-functions)  
   - [5.1 SE_combine](#51-se_combine)  
   - [5.2 SE_impute](#52-se_impute)  
   - [5.3 SE_detectratio](#53-se_detectraio)  
   - [5.4 SE_DEseq2](#54-se_deseq2)  
   - [5.5 SE_boxplot](#55-se_boxplot)  
   - [5.6 SE_PCAplot](#56-se_pcaplot)  
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

### 5.2 SE_impute

Fill missing values in a `SummarizedExperiment` object.

This function imputes missing values (NA) in the given `SummarizedExperiment` object using specified methods.

Multiple imputation techniques can be utilized to handle missing values, ensuring the robustness of subsequent analyses.

#### Usage

    SE_impute(object, assayname = "TPM", group = "group", ZerosAsNA = FALSE, RemoveNA = TRUE,  
              cutoff = 20, method = c("none", "LOD", "half_min", "median", "mean", "min", "knn", "rf", "global_mean", "svd", "QRILC"),  
              LOD = NULL, knum = 10)  

#### Arguments

`object`: A `SummarizedExperiment` object containing the data to be imputed.

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

### 5.3 SE_detectraio

Calculate detection ratio and update `SummarizedExperiment` object's `rowData`.

This function computes the detection ratio of expression data and updates the `rowData` of the provided `SummarizedExperiment` object with detection sample counts and ratios. It also generates a histogram of detection ratios.

#### Usage

    SE_detectratio(SE, assayname = "TPM")  

#### Arguments

`SE`: A `SummarizedExperiment` object containing expression data.

`assayname`: The name of the assay to be used for calculations. Default is `"TPM"`.

### 5.4 SE_DEseq2

Perform differential expression analysis using `DESeq2`.

This function performs differential expression analysis on count data contained in a `SummarizedExperiment` object using the `DESeq2` package.

#### Usage

    SE_DEseq2(SE, assayname = "Count", groupname = "group")  

#### Arguments

`SE`: A `SummarizedExperiment` object containing count data.

`assayname`: The name of the assay to use for the analysis. Default is `"Count"`.

`groupname`: The name of the column in `colData(SE)` that contains the factor for grouping samples. Default is `"group"`.

### 5.5 SE_boxplot

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

### 5.6 SE_PCAplot

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
