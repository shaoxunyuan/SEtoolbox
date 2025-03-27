# SEtoolbox

SEtoolbox is a comprehensive toolkit for SummarizedExperiment analysis

# Installation

You can install the package directly from GitHub,

```r
# install.packages("devtools")
devtools::install_github("shaoxunyuan/SEtoolbox")
```

To run the sample code in our [vignette](
https://shaoxunyuan.github.io/SEtoolbox/
), set the `dependencies` parameter to `TRUE`,

```r
# install.packages("devtools")

devtools::install_github("shaoxunyuan/SEtoolbox", dependencies = TRUE)
```

# Quick Start

## Load example SE object

To load the example `SummarizedExperiment` object, use the following command:  

```r
library(SEtoolbox)  
library(SummarizedExperiment)

# Single SE object
SE  = loadSE()

# List of SE object. A list contain 3 SE object.
SElist  = loadSElist()
```

## Run functions using example SE object

Now you can use SE with all functions 

### SE_combine

Combine multiple SummarizedExperiment objects into one SummarizedExperiment object

```r
#A list contains multiple SE objects. All SE objects should have same group column name (para: group_col) and health control label (para: label_healthy).

SElist = loadSElist()

SE_combine(SElist,merge_type = "intersection")  # Keep intersect features across SE objects

SE_combine(SElist,merge_type = "union")  # Keep union features across SE objects
```

### SE_detectratio

Calculate the number of non-zero samples for each feature and update the results in rowData  

```r
library(SEtoolbox)

SE = loadSE()

# Output a list contain a SE object and detectio ratio density plot 

SE_detectratio(SE, assayname = "TPM", group_col = "group")  
```

### SE_impute

Imputes missing values in the given SummarizedExperiment object using specified methods.

```r
library(SEtoolbox)

SE = loadSE()

SE_impute(SE, group_colname = "group", method = "knn")  
```

### SE_COCONUT

Combat using COCONUT method.

```r
library(SEtoolbox)

#A list contains multiple SE objects. All SE objects should have same group column name (para: group_col) and health control label (para: label_healthy).

SElist = loadSElist()

# Return a list contains multiple SE objects after combat. Retain the original row and column information.

SElistbatch = SE_COCONUT(SElist, assayname = "TPM", group_col = "group", label_healthy = "HC")
```

### SE_DEseq2

Use the count matrix to calculate differential results, with results updated in rowData  

```r
SE_DEseq2(SE)  
```

### SE_distribution

Distribution plot of missing value

```r
SE_distribution(SE)  
```

### SE_PCAplot

Select specific features for PCA plot  

```r
SE_PCAplot(SE)  
```

### SE_heatmap

Select specific features for heatmap 

```r
SE_heatmap(SE)  
```

# Vignette

For full use of our package, please refer to our [documents](
https://shaoxunyuan.github.io/SEtoolbox/)
(highly recommended). 

# Bug Reports

If you encounter any issues, have questions, or would like to make suggestions, 
feel free to report them at our 

[bug tracker](https://github.com/shaoxunyuan/SEtoolbox/issues).

# Contact

For additional inquiries, contact us at: 
Email: shaoxunyuan@njucm.edu.cn