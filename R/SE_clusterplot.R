#' @title SE_clusterplot: Visualize Clustering Results
#' @description This function creates visualizations for clustering results from a SummarizedExperiment object. It can create PCA plots, heatmaps, and dendrograms colored by cluster assignments.
#' @param SE A \code{SummarizedExperiment} object containing cluster assignments in colData or rowData.
#' @param assayname A string indicating which assay to use for visualization. The default value is \code{"log2"}.
#' @param plot_type A character string specifying the type of plot. Options include "pca", "heatmap", "dendrogram". Default is "pca".
#' @param cluster_by A character string specifying what was clustered. Options include "samples" or "features". Default is "samples".
#' @param nfeatures Numeric value indicating the number of top variable features to use for PCA or heatmap. Default is 1000.
#' @param scale_data Logical value indicating whether to scale the data before visualization. Default is TRUE.
#' @param cluster_col A string indicating the column name in colData or rowData containing cluster assignments. Default is "cluster".
#' @return A \code{ggplot} object or a list containing plot objects.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform K-means clustering on samples
#' SE_clustered <- SE_kmeans(SE, assayname = "log2", cluster_by = "samples", 
#'                           n_clusters = 3)
#' 
#' # Create PCA plot colored by clusters
#' pca_plot <- SE_clusterplot(SE_clustered, plot_type = "pca", cluster_by = "samples")
#' print(pca_plot)
#' 
#' # Create heatmap colored by clusters
#' heatmap_plot <- SE_clusterplot(SE_clustered, plot_type = "heatmap", cluster_by = "samples")
#' print(heatmap_plot)
#' @export
SE_clusterplot <- function(SE, assayname = "log2", plot_type = "pca", cluster_by = "samples", 
                            nfeatures = 1000, scale_data = TRUE, cluster_col = "cluster") {
    
    exp_data <- assay(SE, assayname)
    
    if (scale_data) {
        if (cluster_by == "samples") {
            exp_data <- t(scale(t(exp_data)))
        } else {
            exp_data <- scale(exp_data)
        }
    }
    
    if (plot_type == "pca") {
        variances <- apply(exp_data, 1, var, na.rm = TRUE)
        top_features <- names(sort(variances, decreasing = TRUE))[1:min(nfeatures, length(variances))]
        exp_data_subset <- exp_data[top_features, ]
        
        pca_result <- prcomp(t(exp_data_subset), scale. = FALSE, center = TRUE)
        
        pca_df <- data.frame(
            PC1 = pca_result$x[, 1],
            PC2 = pca_result$x[, 2],
            PC3 = pca_result$x[, 3],
            Sample = rownames(pca_result$x)
        )
        
        if (cluster_by == "samples" && cluster_col %in% colnames(colData(SE))) {
            pca_df$Cluster <- colData(SE)[[cluster_col]]
        }
        
        cluster_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
            geom_point(size = 3, alpha = 0.7) +
            labs(title = "PCA Plot Colored by Clusters",
                 x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
                 y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)")) +
            theme_minimal() +
            theme(legend.position = "right")
        
    } else if (plot_type == "heatmap") {
        variances <- apply(exp_data, 1, var, na.rm = TRUE)
        top_features <- names(sort(variances, decreasing = TRUE))[1:min(nfeatures, length(variances))]
        exp_data_subset <- exp_data[top_features, ]
        
        if (cluster_by == "samples" && cluster_col %in% colnames(colData(SE))) {
            annotation_col <- data.frame(Cluster = colData(SE)[[cluster_col]])
            rownames(annotation_col) <- colnames(exp_data_subset)
            
            pheatmap(exp_data_subset, 
                     annotation_col = annotation_col,
                     show_rownames = FALSE,
                     cluster_rows = TRUE,
                     cluster_cols = TRUE)
        } else {
            pheatmap(exp_data_subset, 
                     show_rownames = FALSE,
                     cluster_rows = TRUE,
                     cluster_cols = TRUE)
        }
        
        cluster_plot <- NULL
        
    } else if (plot_type == "dendrogram") {
        if (cluster_by == "samples") {
            dist_matrix <- dist(t(exp_data))
        } else {
            dist_matrix <- dist(exp_data)
        }
        
        hc <- hclust(dist_matrix)
        
        if (cluster_by == "samples" && cluster_col %in% colnames(colData(SE))) {
            dend_colors <- as.numeric(colData(SE)[[cluster_col]])
            dendextend::labels_colors(hc) <- dend_colors
        }
        
        cluster_plot <- as.dendrogram(hc)
        
    } else {
        stop("Unknown plot type. Available options: pca, heatmap, dendrogram")
    }
    
    cat("Cluster plot created:", plot_type, "\n")
    
    return(cluster_plot)
}
