#' @title SE_convert: Convert SummarizedExperiment to Other Data Structures
#' @description This function converts a SummarizedExperiment object to other data structures such as ExpressionSet, data frame, or matrix.
#' @param SE A \code{SummarizedExperiment} object to convert.
#' @param assayname A string indicating which assay to convert. The default value is \code{"TPM"}.
#' @param target_format A character string specifying the target format. Options include "ExpressionSet", "data.frame", "matrix", "list". Default is "data.frame".
#' @param include_metadata Logical value indicating whether to include metadata in the output. Default is TRUE.
#' @return An object in the specified target format.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Convert to data frame
#' df <- SE_convert(SE, assayname = "TPM", target_format = "data.frame")
#' 
#' # Convert to matrix
#' mat <- SE_convert(SE, assayname = "TPM", target_format = "matrix")
#' 
#' # Convert to list
#' lst <- SE_convert(SE, assayname = "TPM", target_format = "list")
#' @export
SE_convert <- function(SE, assayname = "TPM", target_format = "data.frame", 
                      include_metadata = TRUE) {
    
    exp_data <- assay(SE, assayname)
    
    if (target_format == "ExpressionSet") {
        if (!requireNamespace("Biobase", quietly = TRUE)) {
            stop("Biobase package is required for ExpressionSet conversion")
        }
        
        pheno_data <- AnnotatedDataFrame(colData(SE))
        feature_data <- AnnotatedDataFrame(rowData(SE))
        
        eset <- Biobase::ExpressionSet(assayData = exp_data,
                                       phenoData = pheno_data,
                                       featureData = feature_data)
        
        cat("Converted to ExpressionSet\n")
        cat("Features:", nrow(eset), "\n")
        cat("Samples:", ncol(eset), "\n")
        
        return(eset)
        
    } else if (target_format == "data.frame") {
        if (include_metadata) {
            df <- data.frame(exp_data, row.names = rownames(exp_data))
            df <- cbind(df, rowData(SE))
        } else {
            df <- data.frame(exp_data, row.names = rownames(exp_data))
        }
        
        cat("Converted to data.frame\n")
        cat("Dimensions:", dim(df), "\n")
        
        return(df)
        
    } else if (target_format == "matrix") {
        cat("Converted to matrix\n")
        cat("Dimensions:", dim(exp_data), "\n")
        
        return(as.matrix(exp_data))
        
    } else if (target_format == "list") {
        result <- list(
            expression = as.matrix(exp_data),
            rowdata = as.data.frame(rowData(SE)),
            coldata = as.data.frame(colData(SE))
        )
        
        cat("Converted to list\n")
        cat("Expression dimensions:", dim(result$expression), "\n")
        
        return(result)
        
    } else {
        stop("Unknown target_format. Available options: ExpressionSet, data.frame, matrix, list")
    }
}
