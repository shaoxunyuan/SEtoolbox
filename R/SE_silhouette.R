#' @title SE_silhouette: Calculate Silhouette Width for Clustering
#' @description This function calculates silhouette widths for clustering results from a SummarizedExperiment object. Silhouette width measures how well each object lies within its cluster and provides a measure of clustering quality.
#' @param SE A \code{SummarizedExperiment} object containing cluster assignments in colData or rowData.
#' @param assayname A string indicating which assay to use for silhouette calculation. The default value is \code{"log2"}.
#' @param cluster_by A character string specifying what was clustered. Options include "samples" or "features". Default is "samples".
#' @param cluster_col A string indicating the column name in colData or rowData containing cluster assignments. Default is "cluster".
#' @param distance_method A character string specifying the distance metric. Options include "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "correlation". Default is "euclidean".
#' @return A list containing:  
#' \item{silhouette_widths}{A numeric vector of silhouette widths for each sample or feature.}  
#' \item{mean_silhouette}{The mean silhouette width across all samples or features.}  
#' \item{plot}{A \code{ggplot} object showing silhouette widths by cluster.}
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform K-means clustering on samples
#' SE_clustered <- SE_kmeans(SE, assayname = "log2", cluster_by = "samples", 
#'                           n_clusters = 3)
#' 
#' # Calculate silhouette widths
#' silhouette_result <- SE_silhouette(SE_clustered, cluster_by = "samples")
#' 
#' # View mean silhouette width
#' print(silhouette_result$mean_silhouette)
#' 
#' # Plot silhouette widths
#' print(silhouette_result$plot)
#' @export
SE_silhouette <- function(SE, assayname = "log2", cluster_by = "samples", 
                           cluster_col = "cluster", distance_method = "euclidean") {
    
    exp_data <- assay(SE, assayname)
    
    if (cluster_by == "samples") {
        if (!cluster_col %in% colnames(colData(SE))) {
            stop(paste0("Column '", cluster_col, "' not found in colData"))
        }
        clusters <- as.numeric(colData(SE)[[cluster_col]])
        dist_matrix <- dist(t(exp_data), method = distance_method)
    } else if (cluster_by == "features") {
        if (!cluster_col %in% colnames(rowData(SE))) {
            stop(paste0("Column '", cluster_col, "' not found in rowData"))
        }
        clusters <- as.numeric(rowData(SE)[[cluster_col]])
        dist_matrix <- dist(exp_data, method = distance_method)
    } else {
        stop("cluster_by must be either 'samples' or 'features'")
    }
    
    sil <- silhouette(clusters, dist_matrix)
    
    sil_df <- data.frame(
        id = rownames(sil),
        cluster = as.factor(sil[, 1]),
        neighbor = as.factor(sil[, 2]),
        silhouette_width = sil[, 3]
    )
    
    mean_silhouette <- mean(sil[, 3])
    
    sil_plot <- ggplot(sil_df, aes(x = cluster, y = silhouette_width, fill = cluster)) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
        geom_hline(yintercept = mean_silhouette, linetype = "dashed", color = "red") +
        labs(title = "Silhouette Width by Cluster",
             x = "Cluster",
             y = "Silhouette Width",
             fill = "Cluster") +
        theme_minimal() +
        theme(legend.position = "right")
    
    cat("Silhouette analysis completed\n")
    cat("Mean silhouette width:", round(mean_silhouette, 3), "\n")
    cat("Silhouette width range:", round(min(sil[, 3]), 3), "-", round(max(sil[, 3]), 3), "\n")
    
    if (mean_silhouette > 0.7) {
        cat("Strong structure found\n")
    } else if (mean_silhouette > 0.5) {
        cat("Reasonable structure found\n")
    } else if (mean_silhouette > 0.25) {
        cat("Weak structure found\n")
    } else {
        cat("No substantial structure found\n")
    }
    
    return(list(
        silhouette_widths = sil,
        mean_silhouette = mean_silhouette,
        plot = sil_plot
    ))
}
