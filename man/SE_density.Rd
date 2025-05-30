% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SE_density.R
\name{SE_density}
\alias{SE_density}
\title{Plot Density and Boxplots of Gene Expression}
\usage{
SE_density(SE, assayname = "TPM", group_colname = NULL)
}
\arguments{
\item{SE}{A SummarizedExperiment object that contains assay data.}

\item{assayname}{A character string specifying the assay to be used.
Default is "TPM".}

\item{group_colname}{A character string specifying the column name in
the sample metadata that will be used for grouping.   
If NULL (default), a combined density plot for all samples is created.}
}
\value{
A list containing two ggplot2 objects:   
\item{plot_density}{a density plot of gene expression.}  
\item{plot_boxplot}{a boxplot of log2-transformed expression values for each sample.}
}
\description{
This function generates density plots and boxplots of gene expression values   
from a SummarizedExperiment object, allowing for optional grouping of samples.
}
\details{
The expression values are transformed to a log2 scale before   
plotting. Density plots can be generated for individual groups if a   
grouping column is provided; otherwise, a global density plot is shown.   
Boxplots display the distribution of log2-transformed expression values for each sample,   
highlighting any outliers.
}
\examples{
 
# Load necessary libraries  
library(SummarizedExperiment)  
 
# Create a sample SummarizedExperiment object  
# (Sample code for creation)  
# SE <- SummarizedExperiment(assays = list(TPM = matrix(rnorm(100), nrow=10)),   
#                             colData = data.frame(BioSample = 1:10, Group = c("A", "B")))   
 
# Generate density and boxplot  
# plots <- SE_density(SE, assayname = "TPM", group_colname = "Group")  
}
