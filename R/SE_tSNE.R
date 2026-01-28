#' @title SE_tSNE: t-SNE Dimensionality Reduction and Visualization
#' @description This function performs t-Distributed Stochastic Neighbor Embedding (t-SNE) dimensionality reduction on a SummarizedExperiment object and creates visualizations.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for t-SNE. The default value is \code{"log2"}.
#' @param nfeatures Numeric value indicating the number of top variable features to use. Default is 1000.
#' @param perplexity Numeric value for t-SNE perplexity parameter. Typical values are between 5 and 50. Default is 30.
#' @param n_iter Numeric value for the number of iterations. Default is 1000.
#' @param scale_data Logical value indicating whether to scale the data before t-SNE. Default is TRUE.
#' @param color_by A string indicating the column name in colData to use for coloring. Default is NULL.
#' @param shape_by A string indicating the column name in colData to use for shaping. Default is NULL.
#' @return A list containing:  
#' \item{plot}{A \code{ggplot} object showing t-SNE visualization.}  
#' \item{tsne_coords}{A data frame containing t-SNE coordinates.}
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform t-SNE
#' tsne_result <- SE_tSNE(SE, assayname = "log2", color_by = "group")
#' 
#' # View t-SNE plot
#' print(tsne_result$plot)
#' 
#' # Access t-SNE coordinates
#' print(tsne_result$tsne_coords)
#' @export
SE_tSNE <- function(SE, assayname = "log2", nfeatures = 1000, perplexity = 30, 
                     n_iter = 1000, scale_data = TRUE, color_by = NULL, shape_by = NULL) {
    
    exp_data <- assay(SE, assayname)
    
    variances <- apply(exp_data, 1, var, na.rm = TRUE)
    top_features <- names(sort(variances, decreasing = TRUE))[1:min(nfeatures, length(variances))]
    exp_data_subset <- exp_data[top_features, ]
    
    if (scale_data) {
        exp_data_subset <- t(scale(t(exp_data_subset)))
    }
    
    set.seed(123)
    tsne_result <- Rtsne::Rtsne(t(exp_data_subset), 
                                dims = 2, 
                                perplexity = perplexity, 
                                theta = 0.5, 
                                verbose = FALSE, 
                                max_iter = n_iter)
    
    tsne_df <- data.frame(
        tSNE1 = tsne_result$Y[, 1],
        tSNE2 = tsne_result$Y[, 2],
        Sample = colnames(exp_data_subset)
    )
    
    if (!is.null(color_by) && color_by %in% colnames(colData(SE))) {
        tsne_df[[color_by]] <- colData(SE)[[color_by]]
    }
    
    if (!is.null(shape_by) && shape_by %in% colnames(colData(SE))) {
        tsne_df[[shape_by]] <- colData(SE)[[shape_by]]
    }
    
    if (!is.null(color_by) && !is.null(shape_by)) {
        tsne_plot <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = .data[[color_by]], shape = .data[[shape_by]])) +
            geom_point(size = 3, alpha = 0.7) +
            labs(title = "t-SNE Visualization",
                 x = "t-SNE 1",
                 y = "t-SNE 2") +
            theme_minimal() +
            theme(legend.position = "right")
    } else if (!is.null(color_by)) {
        tsne_plot <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = .data[[color_by]])) +
            geom_point(size = 3, alpha = 0.7) +
            labs(title = "t-SNE Visualization",
                 x = "t-SNE 1",
                 y = "t-SNE 2") +
            theme_minimal() +
            theme(legend.position = "right")
    } else {
        tsne_plot <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2)) +
            geom_point(size = 3, alpha = 0.7) +
            labs(title = "t-SNE Visualization",
                 x = "t-SNE 1",
                 y = "t-SNE 2") +
            theme_minimal()
    }
    
    cat("t-SNE analysis completed\n")
    cat("Number of features used:", length(top_features), "\n")
    cat("Perplexity:", perplexity, "\n")
    cat("Number of iterations:", n_iter, "\n")
    
    return(list(
        plot = tsne_plot,
        tsne_coords = tsne_df
    ))
}
