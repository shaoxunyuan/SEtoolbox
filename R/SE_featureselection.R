#' @title SE_featureselection: Feature Selection for Machine Learning
#' @description This function performs feature selection on a SummarizedExperiment object using various methods including variance filtering, correlation filtering, and recursive feature elimination.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for feature selection. The default value is \code{"log2"}.
#' @param group_colname A string representing the column name in \code{colData} that contains group information. Default is "group".
#' @param method A character string specifying the feature selection method. Options include "variance", "correlation", "rfe", "limma". Default is "variance".
#' @param nfeatures Numeric value indicating the number of top features to select. Default is 50.
#' @param correlation_threshold Numeric value for correlation threshold (only used if method = "correlation"). Default is 0.9.
#' @param scale_data Logical value indicating whether to scale the data before feature selection. Default is TRUE.
#' @return A \code{SummarizedExperiment} object containing only the selected features.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Select top 50 features by variance
#' SE_selected <- SE_featureselection(SE, assayname = "log2", group_colname = "group", 
#'                                   method = "variance", nfeatures = 50)
#' 
#' # Select features by limma differential expression
#' SE_selected <- SE_featureselection(SE, assayname = "log2", group_colname = "group", 
#'                                   method = "limma", nfeatures = 50)
#' @export
SE_featureselection <- function(SE, assayname = "log2", group_colname = "group", 
                               method = "variance", nfeatures = 50, 
                               correlation_threshold = 0.9, scale_data = TRUE) {
    
    exp_data <- assay(SE, assayname)
    groups <- colData(SE)[[group_colname]]
    
    if (scale_data) {
        exp_data <- scale(exp_data)
    }
    
    if (method == "variance") {
        variances <- apply(exp_data, 1, var, na.rm = TRUE)
        selected_features <- names(sort(variances, decreasing = TRUE))[1:min(nfeatures, length(variances))]
        cat("Feature selection by variance completed\n")
        cat("Number of features selected:", length(selected_features), "\n")
        
    } else if (method == "correlation") {
        cor_matrix <- cor(t(exp_data), use = "pairwise.complete.obs")
        cor_matrix[abs(cor_matrix) < correlation_threshold] <- 0
        
        cor_scores <- rowSums(abs(cor_matrix))
        selected_features <- names(sort(cor_scores, decreasing = TRUE))[1:min(nfeatures, length(cor_scores))]
        cat("Feature selection by correlation completed\n")
        cat("Correlation threshold:", correlation_threshold, "\n")
        cat("Number of features selected:", length(selected_features), "\n")
        
    } else if (method == "limma") {
        design <- model.matrix(~0 + groups)
        colnames(design) <- levels(factor(groups))
        fit <- lmFit(exp_data, design)
        
        if (ncol(design) >= 2) {
            contrast <- paste(colnames(design)[2], colnames(design)[1], sep = "-")
            contrast_matrix <- makeContrasts(contrasts = contrast, levels = design)
            fit2 <- contrasts.fit(fit, contrast_matrix)
            fit2 <- eBayes(fit2)
            results <- topTable(fit2, number = Inf)
            selected_features <- rownames(results)[1:min(nfeatures, nrow(results))]
        } else {
            selected_features <- rownames(exp_data)
        }
        cat("Feature selection by limma completed\n")
        cat("Number of features selected:", length(selected_features), "\n")
        
    } else if (method == "rfe") {
        groups_factor <- factor(groups)
        
        rfe_control <- rfeControl(functions = rfFuncs, number = 5, verbose = FALSE)
        
        exp_data_df <- as.data.frame(t(exp_data))
        
        set.seed(123)
        rfe_result <- rfe(x = exp_data_df, y = groups_factor, 
                          sizes = nfeatures, rfeControl = rfe_control)
        
        selected_features <- rfe_result$optVariables
        cat("Feature selection by RFE completed\n")
        cat("Number of features selected:", length(selected_features), "\n")
        
    } else {
        stop("Unknown feature selection method. Available options: variance, correlation, rfe, limma")
    }
    
    SE_selected <- SE[selected_features, ]
    
    return(SE_selected)
}
