#' Impute Missing Values in SummarizedExperiment  
#'  
#' This function imputes missing values in a SummarizedExperiment object using the k-nearest neighbor method.  
#'  
#' @param SE A SummarizedExperiment object with expression data.  
#' @param assayname Name of the assay to use (default: "TPM").  
#' @param group Grouping variable name (default: "group").  
#' @param method Imputation method (default: "knn").  
#'  
#' @return A SummarizedExperiment object with imputed data.  
#'  
#' @import SummarizedExperiment  
#' @import impute  
#'  
SE_impute <- function(SE, assayname = "TPM", group = "group", method = "knn"){  
    library(SummarizedExperiment)  
    library(impute)  
    feature_info <- rowData(SE)  
    sample_info <- colData(SE)  
    mdata <- if (length(SE@metadata) == 0) { NULL } else { SE@metadata }   

    # knn  
    expdata = assay(SE,assayname)  
    expdata[expdata==0] = NA  
    expdata_impute = impute.knn(expdata ,k = 10, rowmax = 0.5, colmax = 1, maxp = 1500, rng.seed=123)  
    
    SEimpute <- SummarizedExperiment(assays = expdata_impute[["data"]], rowData = feature_info, colData = sample_info, metadata = mdata)  
    return(SEimpute)  
}