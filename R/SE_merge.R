#' @title SE_merge: Merge rowData or colData from Multiple SummarizedExperiment Objects
#' @description This function merges rowData or colData from multiple SummarizedExperiment objects into a single object. It is useful for combining metadata from different sources.
#' @param SE_list A list of \code{SummarizedExperiment} objects to merge.
#' @param merge_type A character string specifying what to merge. Options include "rowdata", "coldata", "both". Default is "both".
#' @param merge_method A character string specifying the merge method. Options include "inner", "outer", "left", "right". Default is "outer".
#' @param by A character vector specifying column names to merge by. Default is NULL (merge by row/column names).
#' @return A \code{SummarizedExperiment} object with merged metadata.
#' @examples 
#' # Create example SummarizedExperiment objects
#' SE1 <- loadSE()
#' SE2 <- loadSE()
#' 
#' # Merge colData from multiple SE objects
#' SE_merged <- SE_merge(list(SE1, SE2), merge_type = "coldata", 
#'                         merge_method = "outer")
#' 
#' # Merge both rowData and colData
#' SE_merged <- SE_merge(list(SE1, SE2), merge_type = "both", 
#'                         merge_method = "outer")
#' @export
SE_merge <- function(SE_list, merge_type = "both", merge_method = "outer", by = NULL) {
    
    if (length(SE_list) < 2) {
        stop("At least two SummarizedExperiment objects are required")
    }
    
    SE_merged <- SE_list[[1]]
    
    if (merge_type == "rowdata" || merge_type == "both") {
        merged_rowdata <- rowData(SE_merged)
        
        for (i in 2:length(SE_list)) {
            current_rowdata <- rowData(SE_list[[i]])
            
            if (is.null(by)) {
                common_rows <- intersect(rownames(merged_rowdata), rownames(current_rowdata))
                
                if (merge_method == "inner") {
                    merged_rowdata <- merge(merged_rowdata[common_rows, ], 
                                          current_rowdata[common_rows, ], 
                                          by = "row.names", all = FALSE)
                } else if (merge_method == "outer") {
                    merged_rowdata <- merge(merged_rowdata, current_rowdata, 
                                          by = "row.names", all = TRUE)
                } else if (merge_method == "left") {
                    merged_rowdata <- merge(merged_rowdata, current_rowdata, 
                                          by = "row.names", all.x = TRUE, all.y = FALSE)
                } else if (merge_method == "right") {
                    merged_rowdata <- merge(merged_rowdata, current_rowdata, 
                                          by = "row.names", all.x = FALSE, all.y = TRUE)
                } else {
                    stop("Unknown merge_method. Available options: inner, outer, left, right")
                }
                
                rownames(merged_rowdata) <- merged_rowdata$Row.names
                merged_rowdata$Row.names <- NULL
            } else {
                merged_rowdata <- merge(merged_rowdata, current_rowdata, 
                                      by = by, all = merge_method == "outer")
            }
        }
        
        rowData(SE_merged) <- merged_rowdata
        cat("Merged rowData:", nrow(merged_rowdata), "rows,", 
            ncol(merged_rowdata), "columns\n")
    }
    
    if (merge_type == "coldata" || merge_type == "both") {
        merged_coldata <- colData(SE_merged)
        
        for (i in 2:length(SE_list)) {
            current_coldata <- colData(SE_list[[i]])
            
            if (is.null(by)) {
                common_cols <- intersect(rownames(merged_coldata), rownames(current_coldata))
                
                if (merge_method == "inner") {
                    merged_coldata <- merge(merged_coldata[common_cols, ], 
                                          current_coldata[common_cols, ], 
                                          by = "row.names", all = FALSE)
                } else if (merge_method == "outer") {
                    merged_coldata <- merge(merged_coldata, current_coldata, 
                                          by = "row.names", all = TRUE)
                } else if (merge_method == "left") {
                    merged_coldata <- merge(merged_coldata, current_coldata, 
                                          by = "row.names", all.x = TRUE, all.y = FALSE)
                } else if (merge_method == "right") {
                    merged_coldata <- merge(merged_coldata, current_coldata, 
                                          by = "row.names", all.x = FALSE, all.y = TRUE)
                } else {
                    stop("Unknown merge_method. Available options: inner, outer, left, right")
                }
                
                rownames(merged_coldata) <- merged_coldata$Row.names
                merged_coldata$Row.names <- NULL
            } else {
                merged_coldata <- merge(merged_coldata, current_coldata, 
                                      by = by, all = merge_method == "outer")
            }
        }
        
        colData(SE_merged) <- merged_coldata
        cat("Merged colData:", nrow(merged_coldata), "rows,", 
            ncol(merged_coldata), "columns\n")
    }
    
    cat("Metadata merge completed\n")
    cat("Number of SE objects merged:", length(SE_list), "\n")
    
    return(SE_merged)
}
