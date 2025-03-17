# Introduction to SEtoolbox  

- **Shaoxun Yuan**  

- **Affiliation**: School of Artificial Intelligence and Information Technology, Nanjing University of Chinese Medicine, China  

- **Email**: [yuanshaoxun@njucm.edu.cn](mailto:yuanshaoxun@njucm.edu.cn)

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

 

# Installation  

To install the SEtoolbox package, you first need to install the `devtools` package, which provides functions to facilitate package installation from various sources, including GitHub. Run the following command in your R console:  

```r  
install.packages("devtools")  

devtools::install_github("shaoxunyuan/SEtoolbox")
```

# Load packages required for this tutorial  
During this tutorial, we might need to use a few additional packages. 

Since we specified dependencies = TRUE when installing SEtoolbox package, these additional packages have already been installed.

We can load them directly.
```r  
library(BSgenome.Hsapiens.UCSC.hg19)

library(GenomicRanges)

library(DT)

library(rtracklayer) 
```

# Input data  

For this tutorial, we will be working with a [`SummarizedExperiment`](https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html) object. 

This object is a central data structure in Bioconductor that provides a way to store and work with high-dimensional genomic data, such as RNA-seq counts or other assays.

# Main functions  
Functions in SEtoolbox can be obtain using
```r
help(package="SEtoobox")
```

## SE_combine  

Combine multiple `SummarizedExperiment` objects.  

This function merges multiple `SummarizedExperiment` objects based on the specified merge type (intersection or union) for all assays present in the input list.  

## Usage  

```r  
SE_combine(se_list, merge_type = "intersection") 
```
Arguments
se_list: A list of SummarizedExperiment objects to be combined.

merge_type: A character string specifying the type of merge to perform. Options are:

"intersection" (default): Keep only common features across all objects.
"union": Keep all features.


## functions2  
这里是函数2的内容。  

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