#' @title SE_batchdetect: Detect Batch Effects in SummarizedExperiment Object
#' @description This function detects batch effects in a SummarizedExperiment object using PCA visualization and statistical tests. It helps identify whether batch effects exist in the data before applying batch correction methods.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for batch effect detection. The default value is \code{"TPM"}.
#' @param batch_colname A string representing the column name in \code{colData} that contains batch information. Default is "batch".
#' @param group_colname A string representing the column name in \code{colData} that contains group information (e.g., treatment vs control). Default is "group".
#' @param nfeatures Numeric value indicating the number of most variable features to use for PCA. Default is 1000.
#' @param scale Logical value indicating whether to scale the data before PCA. Default is TRUE.
#' @return A list containing:  
#' \item{plot_pca_batch}{A \code{ggplot} object showing PCA colored by batch.}  
#' \item{plot_pca_group}{A \code{ggplot} object showing PCA colored by group.}  
#' \item{plot_pca_both}{A \code{ggplot} object showing PCA with both batch and group information.}  
#' \item{pca_result}{A list containing PCA results including scores, loadings, and variance explained.}  
#' \item{batch_effect_detected}{Logical value indicating whether batch effects are detected.}
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Detect batch effects
#' result <- SE_batchdetect(SE, assayname = "TPM", batch_colname = "batch", group_colname = "group")
#' 
#' # View PCA plots
#' print(result$plot_pca_batch)
#' print(result$plot_pca_group)
#' print(result$plot_pca_both)
#' 
#' # Check if batch effects were detected
#' if (result$batch_effect_detected) {
#'     cat("Batch effects detected. Consider applying batch correction methods.\n")
#' }
#' @export
SE_batchdetect <- function(SE, assayname = "TPM", batch_colname = "batch", group_colname = "group", nfeatures = 1000, scale = TRUE) {
    
    exp_data <- assay(SE, assayname)
    
    variances <- apply(exp_data, 1, var, na.rm = TRUE)
    top_features <- names(sort(variances, decreasing = TRUE))[1:min(nfeatures, length(variances))]
    exp_data_subset <- exp_data[top_features, ]
    
    pca_result <- prcomp(t(exp_data_subset), scale. = scale, center = TRUE)
    
    pca_df <- data.frame(
        PC1 = pca_result$x[, 1],
        PC2 = pca_result$x[, 2],
        PC3 = pca_result$x[, 3],
        Sample = rownames(pca_result$x)
    )
    
    if (batch_colname %in% colnames(colData(SE))) {
        pca_df$Batch <- colData(SE)[[batch_colname]]
    }
    
    if (group_colname %in% colnames(colData(SE))) {
        pca_df$Group <- colData(SE)[[group_colname]]
    }
    
    plot_pca_batch <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
        geom_point(size = 3, alpha = 0.7) +
        labs(title = "PCA Colored by Batch",
             x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
             y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)")) +
        theme_minimal() +
        theme(legend.position = "right")
    
    plot_pca_group <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
        geom_point(size = 3, alpha = 0.7) +
        labs(title = "PCA Colored by Group",
             x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
             y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)")) +
        theme_minimal() +
        theme(legend.position = "right")
    
    if ("Batch" %in% colnames(pca_df) && "Group" %in% colnames(pca_df)) {
        plot_pca_both <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Batch)) +
            geom_point(size = 3, alpha = 0.7) +
            labs(title = "PCA Colored by Group and Shaped by Batch",
                 x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
                 y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)")) +
            theme_minimal() +
            theme(legend.position = "right")
    } else {
        plot_pca_both <- NULL
    }
    
    batch_effect_detected <- FALSE
    if ("Batch" %in% colnames(pca_df) && length(unique(pca_df$Batch)) > 1) {
        batch_var <- by(pca_df$PC1, pca_df$Batch, var)
        batch_effect_detected <- any(batch_var > mean(batch_var) * 2)
    }
    
    cat("Batch effect detection completed\n")
    cat("Variance explained by PC1:", round(summary(pca_result)$importance[2, 1] * 100, 1), "%\n")
    cat("Variance explained by PC2:", round(summary(pca_result)$importance[2, 2] * 100, 1), "%\n")
    cat("Batch effect detected:", batch_effect_detected, "\n")
    
    return(list(
        plot_pca_batch = plot_pca_batch,
        plot_pca_group = plot_pca_group,
        plot_pca_both = plot_pca_both,
        pca_result = pca_result,
        batch_effect_detected = batch_effect_detected
    ))
}
