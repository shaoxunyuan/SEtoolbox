#' @title SE_kmeans: K-means Clustering Analysis
#' @description This function performs K-means clustering analysis on a SummarizedExperiment object. It can cluster samples or features based on the K-means algorithm.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for clustering. The default value is \code{"log2"}.
#' @param cluster_by A character string specifying what to cluster. Options include "samples" or "features". Default is "samples".
#' @param n_clusters Numeric value indicating the number of clusters. Default is 3.
#' @param nstart Numeric value indicating the number of random starts for K-means. Default is 25.
#' @param iter_max Numeric value indicating the maximum number of iterations. Default is 100.
#' @param scale_data Logical value indicating whether to scale the data before clustering. Default is TRUE.
#' @return A \code{SummarizedExperiment} object with cluster assignments added to colData (if clustering samples) or rowData (if clustering features).
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform K-means clustering on samples
#' SE_clustered <- SE_kmeans(SE, assayname = "log2", cluster_by = "samples", 
#'                           n_clusters = 3)
#' 
#' # View cluster assignments
#' print(colData(SE_clustered)$cluster)
#' @export
SE_kmeans <- function(SE, assayname = "log2", cluster_by = "samples", 
                      n_clusters = 3, nstart = 25, iter_max = 100, scale_data = TRUE) {
    
    exp_data <- assay(SE, assayname)
    
    if (scale_data) {
        if (cluster_by == "samples") {
            exp_data <- t(scale(t(exp_data)))
        } else {
            exp_data <- scale(exp_data)
        }
    }
    
    if (cluster_by == "samples") {
        data_for_clustering <- t(exp_data)
    } else if (cluster_by == "features") {
        data_for_clustering <- exp_data
    } else {
        stop("cluster_by must be either 'samples' or 'features'")
    }
    
    set.seed(123)
    km <- kmeans(data_for_clustering, centers = n_clusters, nstart = nstart, iter.max = iter_max)
    
    clusters <- km$cluster
    
    if (cluster_by == "samples") {
        colData(SE)$cluster <- factor(clusters)
        colData(SE)$cluster_method <- "kmeans"
        cat("K-means clustering completed on samples\n")
        cat("Number of clusters:", n_clusters, "\n")
        cat("Cluster sizes:", table(clusters), "\n")
        cat("Within-cluster sum of squares:", km$tot.withinss, "\n")
    } else {
        rowData(SE)$cluster <- factor(clusters)
        rowData(SE)$cluster_method <- "kmeans"
        cat("K-means clustering completed on features\n")
        cat("Number of clusters:", n_clusters, "\n")
        cat("Cluster sizes:", table(clusters), "\n")
        cat("Within-cluster sum of squares:", km$tot.withinss, "\n")
    }
    
    return(SE)
}
