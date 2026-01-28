#' @title SE_rename: Rename Features or Samples in SummarizedExperiment Object
#' @description This function renames features (rows) or samples (columns) in a SummarizedExperiment object. It can rename based on a mapping or apply a prefix/suffix.
#' @param SE A \code{SummarizedExperiment} object.
#' @param rename_type A character string specifying what to rename. Options include "features", "samples", "both". Default is "features".
#' @param mapping A named character vector or list for renaming. Names are old names, values are new names. Default is NULL.
#' @param prefix A character string to add as prefix. Default is NULL.
#' @param suffix A character string to add as suffix. Default is NULL.
#' @param pattern A character string pattern to replace. Default is NULL.
#' @param replacement A character string to replace pattern with. Default is NULL.
#' @return A renamed \code{SummarizedExperiment} object.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Rename features using mapping
#' mapping <- c("Gene1" = "NewGene1", "Gene2" = "NewGene2")
#' SE_renamed <- SE_rename(SE, rename_type = "features", mapping = mapping)
#' 
#' # Add prefix to sample names
#' SE_renamed <- SE_rename(SE, rename_type = "samples", prefix = "Sample_")
#' 
#' # Replace pattern in feature names
#' SE_renamed <- SE_rename(SE, rename_type = "features", 
#'                         pattern = "Gene", replacement = "Feature")
#' @export
SE_rename <- function(SE, rename_type = "features", mapping = NULL, 
                    prefix = NULL, suffix = NULL, pattern = NULL, replacement = NULL) {
    
    if (rename_type == "features" || rename_type == "both") {
        old_names <- rownames(SE)
        new_names <- old_names
        
        if (!is.null(mapping)) {
            for (old_name in names(mapping)) {
                if (old_name %in% new_names) {
                    new_names[new_names == old_name] <- mapping[[old_name]]
                }
            }
            cat("Renamed", length(mapping), "features using mapping\n")
        }
        
        if (!is.null(prefix)) {
            new_names <- paste0(prefix, new_names)
            cat("Added prefix '", prefix, "' to features\n")
        }
        
        if (!is.null(suffix)) {
            new_names <- paste0(new_names, suffix)
            cat("Added suffix '", suffix, "' to features\n")
        }
        
        if (!is.null(pattern) && !is.null(replacement)) {
            new_names <- gsub(pattern, replacement, new_names)
            cat("Replaced pattern '", pattern, "' with '", replacement, "'\n")
        }
        
        rownames(SE) <- new_names
        cat("Renamed", length(old_names), "features\n")
    }
    
    if (rename_type == "samples" || rename_type == "both") {
        old_names <- colnames(SE)
        new_names <- old_names
        
        if (!is.null(mapping)) {
            for (old_name in names(mapping)) {
                if (old_name %in% new_names) {
                    new_names[new_names == old_name] <- mapping[[old_name]]
                }
            }
            cat("Renamed", length(mapping), "samples using mapping\n")
        }
        
        if (!is.null(prefix)) {
            new_names <- paste0(prefix, new_names)
            cat("Added prefix '", prefix, "' to samples\n")
        }
        
        if (!is.null(suffix)) {
            new_names <- paste0(new_names, suffix)
            cat("Added suffix '", suffix, "' to samples\n")
        }
        
        if (!is.null(pattern) && !is.null(replacement)) {
            new_names <- gsub(pattern, replacement, new_names)
            cat("Replaced pattern '", pattern, "' with '", replacement, "'\n")
        }
        
        colnames(SE) <- new_names
        cat("Renamed", length(old_names), "samples\n")
    }
    
    return(SE)
}
