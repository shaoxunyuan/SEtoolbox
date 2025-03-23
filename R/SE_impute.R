#' Impute missing values in a SummarizedExperiment object
#'
#' This function is used to impute missing values in the expression data of a SummarizedExperiment object.
#' It replaces zero values in the expression data with NA and then performs imputation using the kNN method.
#'
#' @param SE A SummarizedExperiment object containing expression data.
#' @param assayname The name of the expression data matrix to be imputed. Default is "TPM".
#' @param group The column name of sample grouping information. Currently not used in this function. Default is "group".
#' @param method The imputation method. Currently only "knn" is supported. Default is "knn".
#'
#' @return A new SummarizedExperiment object with imputed expression data.
#'
#' @import SummarizedExperiment
#' @import impute
#'
#' @examples
#' # Assume se is a SummarizedExperiment object
#' # SEimpute <- SE_impute(se)
#'
#' @export
SE_impute <- function(SE, assayname = "TPM", group = "group", method = "knn") {  
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