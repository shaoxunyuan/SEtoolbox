#' @title Impute Missing Values in SummarizedExperiment  
#'  
#' @description  
#' This function imputes missing values in a SummarizedExperiment object  
#' using the k-nearest neighbor method. The function replaces zeros in the   
#' expression matrix with NA to prepare for imputation.  
#'  
#' @param SE A SummarizedExperiment object containing expression data.  
#' @param assayname A string specifying the name of the assay to use (default is "TPM").   
#'                  This should match one of the assay names in the SE object.  
#' @param group A string specifying the grouping variable name for the samples (default is "group").  
#' @param method A string specifying the imputation method to use. The default method is "knn".  
#'  
#' @return A SummarizedExperiment object with imputed expression data.  
#'  
#' @import SummarizedExperiment  
#' @import impute  
#'  
#' @examples  
#' # Load necessary libraries  
#' library(SummarizedExperiment)  
#' library(impute)  
#'  
#' # Create a simulated SummarizedExperiment object  
#' set.seed(123)  # For reproducibility  
#' example_data <- matrix(rnorm(100), nrow = 10)  
#' example_data[sample(1:100, 20)] <- 0  # Introduce some zeros  
#' example_se <- SummarizedExperiment(assays = list(TPM = example_data))  
#'  
#' # Impute missing values  
#' imputed_se <- SE_impute(example_se, assayname = "TPM", group = "group", method = "knn")  
#'  
#' @export  
SE_impute <- function(SE, assayname = "TPM", group = "group", method = "knn") {  
    options(warn = -1)
	library(SummarizedExperiment)  
    library(impute)  
    
    # Retrieve feature and sample information  
    feature_info <- rowData(SE)  
    sample_info <- colData(SE)  
    mdata <- if (length(SE@metadata) == 0) { NULL } else { SE@metadata }   
    
    # Prepare expression data  
    expdata <- assay(SE, assayname)  
    expdata[expdata == 0] <- NA  # Replace zeros with NA for imputation  
    
    # Impute missing values using kNN method  
    expdata_impute <- impute.knn(expdata, k = 10, rowmax = 0.5, colmax = 1, maxp = 1500, rng.seed = 123)  
    
    # Create a new SummarizedExperiment object with imputed data  
    SEimpute <- SummarizedExperiment(assays = expdata_impute[["data"]], rowData = feature_info, colData = sample_info, metadata = mdata)  
    
    return(SEimpute)  
}