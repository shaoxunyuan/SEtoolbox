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

To create the example `SummarizedExperiment` object, use the following command:  

```r
library(SEtoolbox)  

# Single SE object
data_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(data_matrix) <- paste0("Gene", 1:100)
colnames(data_matrix) <- paste0("Sample", 1:10)
sample_info <- DataFrame(group = rep(c("A", "B"), each = 5))
SE <- SummarizedExperiment(assays = list(TPM = data_matrix), colData = sample_info) 

# List of SE object
SElist <- vector("list", 3)  
for (i in 1:3) {  
  data_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)  
  rownames(data_matrix) <- paste0("Gene", 1:100)  
  colnames(data_matrix) <- paste0("Sample", 1:10) 
  sample_info <- DataFrame(group = rep(c("A", "B"), each = 5))  
  SElist[[i]] <- SummarizedExperiment(assays = list(TPM = data_matrix), colData = sample_info)   
  names(SElist)[i] <- paste0("SE", i)  
}  
```



Now you can use SE with all functions 

Combine multiple SummarizedExperiment objects into one SummarizedExperiment object

```r
SE_Combine(list(SE1,SE2,SE3))  
```

Calculate the number of non-zero samples for each feature and update the results in rowData  

```r
SE_detectratio(SE)  
```

Imputes missing values (NA) in the given SummarizedExperiment object using specified methods.

```r
SE_impute(SE)  
```

Use the count matrix to calculate differential results, with results updated in rowData  

```r
SE_DEseq2(SE)  
```

Distribution plot of missing value

```r
SE_distribution(SE)  
```

Select specific features for PCA plot  

```r
SE_PCAplot(SE)  
```

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
