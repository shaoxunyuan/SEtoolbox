# SEtoolbox <img src="./vignettes/images/logo.png" align = "right" width = "150" />

SEtoolbox is an R package for SummarizedExperiment analysis and visualization.

# Installation
### Option 1: Install from GitHub

You can install the package directly from GitHub,
```r
# install.packages("devtools")
devtools::install_github("shaoxunyuan/SEtoolbox")
```

To run the sample code in our [vignette](
https://bioconductor.org/packages/devel/bioc/vignettes/G4SNVHunter/inst/doc/G4SNVHunter.html
), set the `dependencies` parameter to `TRUE`,
```r
# install.packages("devtools")
devtools::install_github("shaoxunyuan/SEtoolbox", dependencies = TRUE)
```

You can check all available functions
```r
ls("package:SEtoolbox") 
getNamespaceExports("SEtoolbox")
```