# Introduction to SEtoolbox

* **Shaoxun Yuan**
  
* **Affiliation**: School of Artificial Intelligence and Information Technology, Nanjing University of Chinese Medicine, China
  
* **Email**: [yuanshaoxun@njucm.edu.cn](mailto:yuanshaoxun@njucm.edu.cn)
  

# Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Load packages required for this tutorial](#load-packages-required-for-this-tutorial)
4. [Input data](#input-data)
5. [Main functions](#main-functions)
  1. [functions1](#functions1)
  2. [functions2](#functions2)
  3. [functions3](#functions3)
6. [A coherent example](#a-coherent-example)
7. [Acknowledgements](#acknowledgements)
8. [Session Info](#session-info)
9. [References](#references)

# Introduction

SEtoolbox is an R package that operates, analyzes and visualizes SummarizedExperiment objects.

# Installation

To install the SEtoolbox package, you first need to install the `devtools` package, which provides functions to facilitate package installation from various sources, including GitHub. Run the following command in your R console:

    install.packages("devtools")  
    
    devtools::install_github("shaoxunyuan/SEtoolbox")

# Load packages required for this tutorial

During this tutorial, we might need to use a few additional packages.

Since we specified dependencies = TRUE when installing SEtoolbox package, these additional packages have already been installed.

We can load them directly.

    library(BSgenome.Hsapiens.UCSC.hg19)
    
    library(GenomicRanges)
    
    library(DT)
    
    library(rtracklayer) 

# Input data

For this tutorial, SEtoolbox will be working with a [`SummarizedExperiment`](https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html) object.

# Main functions

Functions in SEtoolbox can be obtain using

    help(package="SEtoolbox")

## SE_combine

Combine multiple `SummarizedExperiment` objects.

This function merges multiple `SummarizedExperiment` objects based on the specified merge type (intersection or union) for all assays present in the input list.

#### #### Usage

    SE_combine(se_list, merge_type = "intersection") 

#### #### Arguments

`se_list`: A list of SummarizedExperiment objects to be combined.

`merge_type`: A character string specifying the type of merge to perform. Options are:

    1.`intersection` (default): Keep only common features across all objects.

    2.`union`: Keep all features.

## SE_impute

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

## functions3

这里是函数3的内容。

# A coherent example

这里是完整示例的内容。

# Acknowledgements

这里是致谢的内容。

# Session Info

这里是会话信息的内容。

# References

这里是参考文献的内容。