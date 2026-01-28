#' @title SE_transform: Transform Expression Data in SummarizedExperiment Object
#' @description This function applies various transformations to expression data in a SummarizedExperiment object, including log2, log10, natural log, square root, and inverse transformations.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to transform. The default value is \code{"TPM"}.
#' @param method A character string specifying the transformation method. Options include "log2", "log10", "ln", "sqrt", "inverse", "arcsinh". Default is "log2".
#' @param pseudocount Numeric value to add before transformation to avoid log(0) or division by zero. Default is 1.
#' @param new_assayname A string indicating the name for the new transformed assay. If NULL, it will be named as "{assayname}_{method}". Default is NULL.
#' @return A transformed \code{SummarizedExperiment} object with the transformed data stored in a new assay.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Log2 transformation
#' SE_log2 <- SE_transform(SE, assayname = "TPM", method = "log2")
#' 
#' # Log10 transformation
#' SE_log10 <- SE_transform(SE, assayname = "TPM", method = "log10")
#' 
#' # Square root transformation
#' SE_sqrt <- SE_transform(SE, assayname = "TPM", method = "sqrt")
#' 
#' # Arcsinh transformation (useful for data with zeros and negative values)
#' SE_arcsinh <- SE_transform(SE, assayname = "TPM", method = "arcsinh")
#' @export
SE_transform <- function(SE, assayname = "TPM", method = "log2", pseudocount = 1, new_assayname = NULL) {
    
    exp_data <- assay(SE, assayname)
    
    if (method == "log2") {
        transformed_data <- log2(exp_data + pseudocount)
        cat("Log2 transformation applied with pseudocount =", pseudocount, "\n")
    } else if (method == "log10") {
        transformed_data <- log10(exp_data + pseudocount)
        cat("Log10 transformation applied with pseudocount =", pseudocount, "\n")
    } else if (method == "ln") {
        transformed_data <- log(exp_data + pseudocount)
        cat("Natural log transformation applied with pseudocount =", pseudocount, "\n")
    } else if (method == "sqrt") {
        transformed_data <- sqrt(pmax(exp_data, 0))
        cat("Square root transformation applied\n")
    } else if (method == "inverse") {
        transformed_data <- 1 / (exp_data + pseudocount)
        cat("Inverse transformation applied with pseudocount =", pseudocount, "\n")
    } else if (method == "arcsinh") {
        transformed_data <- asinh(exp_data)
        cat("Arcsinh transformation applied\n")
    } else {
        stop("Unknown transformation method. Available methods: log2, log10, ln, sqrt, inverse, arcsinh")
    }
    
    if (is.null(new_assayname)) {
        new_assayname <- paste0(assayname, "_", method)
    }
    
    assay(SE, new_assayname) <- transformed_data
    cat("Transformed data stored in assay:", new_assayname, "\n")
    
    return(SE)
}
