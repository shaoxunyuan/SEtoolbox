#' @title SE_correlation: Correlation Analysis
#' @description This function performs correlation analysis on a SummarizedExperiment object. It can compute correlations between features, between samples, or between features and sample traits.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for correlation analysis. The default value is \code{"log2"}.
#' @param method A character string specifying what to correlate. Options include "features", "samples", "feature_trait". Default is "features".
#' @param correlation_method A character string specifying the correlation method. Options include "pearson", "spearman", "kendall". Default is "pearson".
#' @param trait_col A string indicating the column name in colData containing trait information (only used if method = "feature_trait"). Default is NULL.
#' @return A correlation matrix or a list containing correlation results.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Compute feature-feature correlations
#' cor_result <- SE_correlation(SE, assayname = "log2", method = "features")
#' 
#' # Compute sample-sample correlations
#' cor_result <- SE_correlation(SE, assayname = "log2", method = "samples")
#' 
#' # Compute feature-trait correlations
#' cor_result <- SE_correlation(SE, assayname = "log2", method = "feature_trait", 
#'                              trait_col = "group")
#' @export
SE_correlation <- function(SE, assayname = "log2", method = "features", 
                          correlation_method = "pearson", trait_col = NULL) {
    
    exp_data <- assay(SE, assayname)
    
    if (method == "features") {
        cor_matrix <- cor(exp_data, method = correlation_method, use = "pairwise.complete.obs")
        cat("Feature-feature correlation analysis completed\n")
        cat("Correlation method:", correlation_method, "\n")
        cat("Correlation matrix dimensions:", dim(cor_matrix), "\n")
        return(cor_matrix)
        
    } else if (method == "samples") {
        cor_matrix <- cor(t(exp_data), method = correlation_method, use = "pairwise.complete.obs")
        cat("Sample-sample correlation analysis completed\n")
        cat("Correlation method:", correlation_method, "\n")
        cat("Correlation matrix dimensions:", dim(cor_matrix), "\n")
        return(cor_matrix)
        
    } else if (method == "feature_trait") {
        if (is.null(trait_col)) {
            stop("trait_col must be specified when method = 'feature_trait'")
        }
        
        if (!trait_col %in% colnames(colData(SE))) {
            stop(paste0("Column '", trait_col, "' not found in colData"))
        }
        
        trait <- colData(SE)[[trait_col]]
        
        if (is.numeric(trait)) {
            cor_values <- apply(exp_data, 1, function(x) cor(x, trait, method = correlation_method))
        } else {
            trait_factor <- as.factor(trait)
            cor_values <- apply(exp_data, 1, function(x) {
                group_means <- tapply(x, trait_factor, mean, na.rm = TRUE)
                cor(group_means, as.numeric(names(group_means)), method = correlation_method)
            })
        }
        
        cor_df <- data.frame(
            feature = names(cor_values),
            correlation = cor_values,
            stringsAsFactors = FALSE
        )
        
        cor_df <- cor_df[order(abs(cor_df$correlation), decreasing = TRUE), ]
        rownames(cor_df) <- NULL
        
        cat("Feature-trait correlation analysis completed\n")
        cat("Correlation method:", correlation_method, "\n")
        cat("Trait column:", trait_col, "\n")
        
        return(cor_df)
        
    } else {
        stop("Unknown method. Available options: features, samples, feature_trait")
    }
}
