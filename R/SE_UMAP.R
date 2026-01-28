#' @title SE_UMAP: UMAP Dimensionality Reduction and Visualization
#' @description This function performs Uniform Manifold Approximation and Projection (UMAP) dimensionality reduction on a SummarizedExperiment object and creates visualizations.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for UMAP. The default value is \code{"log2"}.
#' @param nfeatures Numeric value indicating the number of top variable features to use. Default is 1000.
#' @param n_neighbors Numeric value for UMAP n_neighbors parameter. Controls the balance between local and global structure. Default is 15.
#' @param min_dist Numeric value for UMAP min_dist parameter. Controls how tightly points are packed. Default is 0.1.
#' @param metric A character string specifying the distance metric. Options include "euclidean", "manhattan", "cosine", "correlation". Default is "euclidean".
#' @param scale_data Logical value indicating whether to scale the data before UMAP. Default is TRUE.
#' @param color_by A string indicating the column name in colData to use for coloring. Default is NULL.
#' @param shape_by A string indicating the column name in colData to use for shaping. Default is NULL.
#' @return A list containing:  
#' \item{plot}{A \code{ggplot} object showing UMAP visualization.}  
#' \item{umap_coords}{A data frame containing UMAP coordinates.}
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform UMAP
#' umap_result <- SE_UMAP(SE, assayname = "log2", color_by = "group")
#' 
#' # View UMAP plot
#' print(umap_result$plot)
#' 
#' # Access UMAP coordinates
#' print(umap_result$umap_coords)
#' @export
SE_UMAP <- function(SE, assayname = "log2", nfeatures = 1000, n_neighbors = 15, 
                     min_dist = 0.1, metric = "euclidean", scale_data = TRUE, 
                     color_by = NULL, shape_by = NULL) {
    
    exp_data <- assay(SE, assayname)
    
    variances <- apply(exp_data, 1, var, na.rm = TRUE)
    top_features <- names(sort(variances, decreasing = TRUE))[1:min(nfeatures, length(variances))]
    exp_data_subset <- exp_data[top_features, ]
    
    if (scale_data) {
        exp_data_subset <- t(scale(t(exp_data_subset)))
    }
    
    set.seed(123)
    umap_result <- uwot::umap(t(exp_data_subset), 
                              n_neighbors = n_neighbors, 
                              min_dist = min_dist, 
                              metric = metric)
    
    umap_df <- data.frame(
        UMAP1 = umap_result[, 1],
        UMAP2 = umap_result[, 2],
        Sample = colnames(exp_data_subset)
    )
    
    if (!is.null(color_by) && color_by %in% colnames(colData(SE))) {
        umap_df[[color_by]] <- colData(SE)[[color_by]]
    }
    
    if (!is.null(shape_by) && shape_by %in% colnames(colData(SE))) {
        umap_df[[shape_by]] <- colData(SE)[[shape_by]]
    }
    
    if (!is.null(color_by) && !is.null(shape_by)) {
        umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = .data[[color_by]], shape = .data[[shape_by]])) +
            geom_point(size = 3, alpha = 0.7) +
            labs(title = "UMAP Visualization",
                 x = "UMAP 1",
                 y = "UMAP 2") +
            theme_minimal() +
            theme(legend.position = "right")
    } else if (!is.null(color_by)) {
        umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = .data[[color_by]])) +
            geom_point(size = 3, alpha = 0.7) +
            labs(title = "UMAP Visualization",
                 x = "UMAP 1",
                 y = "UMAP 2") +
            theme_minimal() +
            theme(legend.position = "right")
    } else {
        umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
            geom_point(size = 3, alpha = 0.7) +
            labs(title = "UMAP Visualization",
                 x = "UMAP 1",
                 y = "UMAP 2") +
            theme_minimal()
    }
    
    cat("UMAP analysis completed\n")
    cat("Number of features used:", length(top_features), "\n")
    cat("n_neighbors:", n_neighbors, "\n")
    cat("min_dist:", min_dist, "\n")
    
    return(list(
        plot = umap_plot,
        umap_coords = umap_df
    ))
}
