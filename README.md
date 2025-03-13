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

You can check all available functions,
```r
library(SEtoolbox)  
ls("package:SEtoolbox") 
getNamespaceExports("SEtoolbox")
```

# Quick use

To load the example SummarizedExperiment object, use the following command:  

```r  
library(SEtoolbox)  

# Load the example SummarizedExperiment object  
SE <- loadSE()  

# Now you can use SE with all functions 

SE_detectratio(SE)
SE_boxplot(SE)
SE_DEseq2(SE)
SE_PCAplot(SE)
SE_heatmap(SE)