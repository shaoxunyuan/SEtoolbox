# SEtoolbox <img src="./vignettes/images/logo.png" align = "right" width = "150" />

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

# Quick use

To load the example SummarizedExperiment object, use the following command:  

```r  
library(SEtoolbox)  

# Load the example SummarizedExperiment object  
SE <- loadSE()  

# Now you can use SE with all functions 

## Calculate the number of non-zero samples for each feature and update the results in rowData  
SE_detectratio(SE)  

## Select specific features to create a boxplot, supports grouping or non-grouping  
SE_boxplot(SE)  

## Use the count matrix to calculate differential results, with results updated in rowData  
SE_DEseq2(SE)  

## Select specific features for PCA plot  
SE_PCAplot(SE)  

## Select specific features for heatmap  
SE_heatmap(SE)  

```