#' @title SE_hierarchical: Hierarchical Clustering Analysis
#' @description This function performs hierarchical clustering analysis on a SummarizedExperiment object. It can cluster samples or features based on various distance metrics and linkage methods.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for clustering. The default value is \code{"log2"}.
#' @param cluster_by A character string specifying what to cluster. Options include "samples" or "features". Default is "samples".
#' @param distance_method A character string specifying the distance metric. Options include "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "correlation". Default is "euclidean".
#' @param hclust_method A character string specifying the linkage method. Options include "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid". Default is "ward.D2".
#' @param n_clusters Numeric value indicating the number of clusters to cut the dendrogram into. Default is 2.
#' @param scale_data Logical value indicating whether to scale the data before clustering. Default is TRUE.
#' @return A \code{SummarizedExperiment} object with cluster assignments added to colData (if clustering samples) or rowData (if clustering features).
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform hierarchical clustering on samples
#' SE_clustered <- SE_hierarchical(SE, assayname = "log2", cluster_by = "samples", 
#'                                  n_clusters = 3)
#' 
#' # View cluster assignments
#' print(colData(SE_clustered)$cluster)
#' @export
SE_hierarchical <- function(SE, assayname = "log2", cluster_by = "samples", 
                            distance_method = "euclidean", hclust_method = "ward.D2", 
                            n_clusters = 2, scale_data = TRUE) {
    
    exp_data <- assay(SE, assayname)
    
    if (scale_data) {
        if (cluster_by == "samples") {
            exp_data <- t(scale(t(exp_data)))
        } else {
            exp_data <- scale(exp_data)
        }
    }
    
    if (cluster_by == "samples") {
        dist_matrix <- dist(t(exp_data), method = distance_method)
    } else if (cluster_by == "features") {
        dist_matrix <- dist(exp_data, method = distance_method)
    } else {
        stop("cluster_by must be either 'samples' or 'features'")
    }
    
    hc <- hclust(dist_matrix, method = hclust_method)
    
    clusters <- cutree(hc, k = n_clusters)
    
    if (cluster_by == "samples") {
        colData(SE)$cluster <- factor(clusters)
        colData(SE)$cluster_method <- "hierarchical"
        cat("Hierarchical clustering completed on samples\n")
        cat("Number of clusters:", n_clusters, "\n")
        cat("Cluster sizes:", table(clusters), "\n")
    } else {
        rowData(SE)$cluster <- factor(clusters)
        rowData(SE)$cluster_method <- "hierarchical"
        cat("Hierarchical clustering completed on features\n")
        cat("Number of clusters:", n_clusters, "\n")
        cat("Cluster sizes:", table(clusters), "\n")
    }
    
    return(SE)
}
