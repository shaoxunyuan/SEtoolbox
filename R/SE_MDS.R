#' @title SE_MDS: Multidimensional Scaling Analysis
#' @description This function performs Multidimensional Scaling (MDS) analysis on a SummarizedExperiment object and creates visualizations. MDS is a technique for visualizing the level of similarity of individual cases of a dataset.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for MDS. The default value is \code{"log2"}.
#' @param nfeatures Numeric value indicating the number of top variable features to use. Default is 1000.
#' @param distance_method A character string specifying the distance metric. Options include "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "correlation". Default is "euclidean".
#' @param scale_data Logical value indicating whether to scale the data before MDS. Default is TRUE.
#' @param color_by A string indicating the column name in colData to use for coloring. Default is NULL.
#' @param shape_by A string indicating the column name in colData to use for shaping. Default is NULL.
#' @return A list containing:  
#' \item{plot}{A \code{ggplot} object showing MDS visualization.}  
#' \item{mds_coords}{A data frame containing MDS coordinates.}  
#' \item{stress}{The stress value indicating goodness of fit.}
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform MDS
#' mds_result <- SE_MDS(SE, assayname = "log2", color_by = "group")
#' 
#' # View MDS plot
#' print(mds_result$plot)
#' 
#' # Access MDS coordinates
#' print(mds_result$mds_coords)
#' 
#' # Check stress value
#' print(mds_result$stress)
#' @export
SE_MDS <- function(SE, assayname = "log2", nfeatures = 1000, distance_method = "euclidean", 
                   scale_data = TRUE, color_by = NULL, shape_by = NULL) {
    
    exp_data <- assay(SE, assayname)
    
    variances <- apply(exp_data, 1, var, na.rm = TRUE)
    top_features <- names(sort(variances, decreasing = TRUE))[1:min(nfeatures, length(variances))]
    exp_data_subset <- exp_data[top_features, ]
    
    if (scale_data) {
        exp_data_subset <- t(scale(t(exp_data_subset)))
    }
    
    dist_matrix <- dist(t(exp_data_subset), method = distance_method)
    
    mds_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)
    
    mds_df <- data.frame(
        MDS1 = mds_result$points[, 1],
        MDS2 = mds_result$points[, 2],
        Sample = colnames(exp_data_subset)
    )
    
    if (!is.null(color_by) && color_by %in% colnames(colData(SE))) {
        mds_df[[color_by]] <- colData(SE)[[color_by]]
    }
    
    if (!is.null(shape_by) && shape_by %in% colnames(colData(SE))) {
        mds_df[[shape_by]] <- colData(SE)[[shape_by]]
    }
    
    if (!is.null(color_by) && !is.null(shape_by)) {
        mds_plot <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = .data[[color_by]], shape = .data[[shape_by]])) +
            geom_point(size = 3, alpha = 0.7) +
            labs(title = "MDS Visualization",
                 x = paste0("MDS 1 (", round(mds_result$GOF[1] * 100, 1), "% variance)"),
                 y = paste0("MDS 2 (", round(mds_result$GOF[2] * 100, 1), "% variance)")) +
            theme_minimal() +
            theme(legend.position = "right")
    } else if (!is.null(color_by)) {
        mds_plot <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = .data[[color_by]])) +
            geom_point(size = 3, alpha = 0.7) +
            labs(title = "MDS Visualization",
                 x = paste0("MDS 1 (", round(mds_result$GOF[1] * 100, 1), "% variance)"),
                 y = paste0("MDS 2 (", round(mds_result$GOF[2] * 100, 1), "% variance)")) +
            theme_minimal() +
            theme(legend.position = "right")
    } else {
        mds_plot <- ggplot(mds_df, aes(x = MDS1, y = MDS2)) +
            geom_point(size = 3, alpha = 0.7) +
            labs(title = "MDS Visualization",
                 x = paste0("MDS 1 (", round(mds_result$GOF[1] * 100, 1), "% variance)"),
                 y = paste0("MDS 2 (", round(mds_result$GOF[2] * 100, 1), "% variance)")) +
            theme_minimal()
    }
    
    stress <- mds_result$eig[mds_result$eig < 0]
    if (length(stress) > 0) {
        stress_value <- sum(abs(stress)) / sum(abs(mds_result$eig))
    } else {
        stress_value <- 0
    }
    
    cat("MDS analysis completed\n")
    cat("Number of features used:", length(top_features), "\n")
    cat("Distance method:", distance_method, "\n")
    cat("Stress value:", round(stress_value, 4), "\n")
    
    return(list(
        plot = mds_plot,
        mds_coords = mds_df,
        stress = stress_value
    ))
}
