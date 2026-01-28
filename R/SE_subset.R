#' @title SE_subset: Subset SummarizedExperiment Object by Features or Samples
#' @description This function subsets a SummarizedExperiment object based on feature names, sample names, or conditions from colData. It provides flexible subsetting options to extract specific subsets of data for downstream analysis.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param features A character vector of feature names to keep. If NULL, all features are kept. Default is NULL.
#' @param samples A character vector of sample names to keep. If NULL, all samples are kept. Default is NULL.
#' @param condition A named list or vector for conditional subsetting. Names should be column names in colData, and values should be the values to keep. Default is NULL.
#' @param exclude_features A character vector of feature names to exclude. Default is NULL.
#' @param exclude_samples A character vector of sample names to exclude. Default is NULL.
#' @return A subsetted \code{SummarizedExperiment} object containing only the specified features and samples.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Subset by specific features
#' SE_subset <- SE_subset(SE, features = c("Gene1", "Gene2", "Gene3"))
#' 
#' # Subset by specific samples
#' SE_subset <- SE_subset(SE, samples = c("Sample1", "Sample2"))
#' 
#' # Subset by condition from colData
#' SE_subset <- SE_subset(SE, condition = list(group = "Treatment"))
#' 
#' # Subset by multiple conditions
#' SE_subset <- SE_subset(SE, condition = list(group = c("Treatment1", "Treatment2"), batch = 1))
#' 
#' # Exclude specific features
#' SE_subset <- SE_subset(SE, exclude_features = c("GeneX", "GeneY"))
#' 
#' # Combine features and condition subsetting
#' SE_subset <- SE_subset(SE, features = c("Gene1", "Gene2"), condition = list(group = "Treatment"))
#' @export
SE_subset <- function(SE, features = NULL, samples = NULL, condition = NULL, 
                      exclude_features = NULL, exclude_samples = NULL) {
    
    keep_features <- rep(TRUE, nrow(SE))
    keep_samples <- rep(TRUE, ncol(SE))
    
    if (!is.null(features)) {
        keep_features <- keep_features & (rownames(SE) %in% features)
        cat(paste0("Keeping ", sum(keep_features), " features specified\n"))
    }
    
    if (!is.null(samples)) {
        keep_samples <- keep_samples & (colnames(SE) %in% samples)
        cat(paste0("Keeping ", sum(keep_samples), " samples specified\n"))
    }
    
    if (!is.null(exclude_features)) {
        keep_features <- keep_features & !(rownames(SE) %in% exclude_features)
        cat(paste0("Excluding ", length(exclude_features), " features\n"))
    }
    
    if (!is.null(exclude_samples)) {
        keep_samples <- keep_samples & !(colnames(SE) %in% exclude_samples)
        cat(paste0("Excluding ", length(exclude_samples), " samples\n"))
    }
    
    if (!is.null(condition)) {
        for (col_name in names(condition)) {
            if (col_name %in% colnames(colData(SE))) {
                col_values <- condition[[col_name]]
                keep_condition <- colData(SE)[[col_name]] %in% col_values
                keep_samples <- keep_samples & keep_condition
                cat(paste0("Keeping ", sum(keep_condition), " samples where ", col_name, " is ", 
                           paste(col_values, collapse = ", "), "\n"))
            } else {
                warning(paste0("Column '", col_name, "' not found in colData. Skipping this condition."))
            }
        }
    }
    
    SE_subset <- SE[keep_features, keep_samples]
    cat(paste0("Final subsetted SE: ", nrow(SE_subset), " features, ", ncol(SE_subset), " samples\n"))
    
    return(SE_subset)
}
