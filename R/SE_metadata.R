#' @title SE_metadata: Extract and Manipulate Metadata from SummarizedExperiment Object
#' @description This function extracts, updates, or manipulates metadata (colData and rowData) from a SummarizedExperiment object.
#' @param SE A \code{SummarizedExperiment} object.
#' @param metadata_type A character string specifying which metadata to extract. Options include "coldata", "rowdata", "both". Default is "both".
#' @param add_coldata A named list or data frame to add to colData. Default is NULL.
#' @param add_rowdata A named list or data frame to add to rowData. Default is NULL.
#' @param remove_coldata A character vector of column names to remove from colData. Default is NULL.
#' @param remove_rowdata A character vector of column names to remove from rowData. Default is NULL.
#' @return A list containing the requested metadata or an updated \code{SummarizedExperiment} object if modifications were made.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Extract all metadata
#' metadata <- SE_metadata(SE, metadata_type = "both")
#' 
#' # View colData
#' print(metadata$coldata)
#' 
#' # View rowData
#' print(metadata$rowdata)
#' 
#' # Add new column to colData
#' SE_updated <- SE_metadata(SE, add_coldata = list(new_col = c(1, 2, 3)))
#' 
#' # Remove column from rowData
#' SE_updated <- SE_metadata(SE, remove_rowdata = c("old_col"))
#' @export
SE_metadata <- function(SE, metadata_type = "both", add_coldata = NULL, 
                       add_rowdata = NULL, remove_coldata = NULL, 
                       remove_rowdata = NULL) {
    
    if (!is.null(add_coldata)) {
        if (is.list(add_coldata)) {
            for (col_name in names(add_coldata)) {
                colData(SE)[[col_name]] <- add_coldata[[col_name]]
            }
            cat("Added", length(add_coldata), "columns to colData\n")
        } else {
            colData(SE) <- cbind(colData(SE), add_coldata)
            cat("Added data frame to colData\n")
        }
    }
    
    if (!is.null(add_rowdata)) {
        if (is.list(add_rowdata)) {
            for (col_name in names(add_rowdata)) {
                rowData(SE)[[col_name]] <- add_rowdata[[col_name]]
            }
            cat("Added", length(add_rowdata), "columns to rowData\n")
        } else {
            rowData(SE) <- cbind(rowData(SE), add_rowdata)
            cat("Added data frame to rowData\n")
        }
    }
    
    if (!is.null(remove_coldata)) {
        for (col_name in remove_coldata) {
            if (col_name %in% colnames(colData(SE))) {
                colData(SE)[[col_name]] <- NULL
            }
        }
        cat("Removed", length(remove_coldata), "columns from colData\n")
    }
    
    if (!is.null(remove_rowdata)) {
        for (col_name in remove_rowdata) {
            if (col_name %in% colnames(rowData(SE))) {
                rowData(SE)[[col_name]] <- NULL
            }
        }
        cat("Removed", length(remove_rowdata), "columns from rowData\n")
    }
    
    result <- list()
    
    if (metadata_type == "coldata" || metadata_type == "both") {
        result$coldata <- as.data.frame(colData(SE))
        cat("Extracted colData:", nrow(result$coldata), "rows,", 
            ncol(result$coldata), "columns\n")
    }
    
    if (metadata_type == "rowdata" || metadata_type == "both") {
        result$rowdata <- as.data.frame(rowData(SE))
        cat("Extracted rowData:", nrow(result$rowdata), "rows,", 
            ncol(result$rowdata), "columns\n")
    }
    
    if (length(result) > 0) {
        return(result)
    } else {
        return(SE)
    }
}
