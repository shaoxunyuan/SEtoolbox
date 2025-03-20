#' Generate PCA plots  
#'  
#' This function takes a SummarizedExperiment object, computes PCA, and visualizes the results.  
#' @param SE SummarizedExperiment object containing gene expression data.  
#' @param assayname Name of the expression data, default is "TPM".  
#' @param groupname Name of the grouping column, default is "group".  
#' @param outlier_threshold Outlier filtering threshold, default is 2.  
#' @param scale Whether to standardize the data, default is TRUE.  
#' @param feature_of_interesting Vector of specific feature names; if NULL, all features are used, default is NULL.  
#' @param show_caption Logical value indicating whether to display caption, default is TRUE.  
#' @return A list containing two ggplot objects: the original PCA plot and the filtered PCA plot.  
#' @export  
SE_PCAplot = function(SE, assayname = "TPM", groupname = "group", outlier_threshold = 2, scale = TRUE, feature_of_interesting = NULL, show_caption = TRUE){  
  library(ggplot2)  
  library(plyr)  
  library(patchwork)  
  library(cluster)  

  SCvalue = function(data){   
    num_clusters <- length(unique(data$group))  
    kmeans_result <- kmeans(data[, c("PC1", "PC2")], centers = num_clusters, nstart = 25)  
    data$cluster <- as.factor(kmeans_result$cluster)  
    distance_matrix <- dist(data[, c("PC1", "PC2")])  
    silhouette_result <- silhouette(as.numeric(data$cluster), distance_matrix)  
    mean_silhouette_width <- round(mean(silhouette_result[, "sil_width"]), 3)  
    mean_silhouette_width  
  }  

  expdata = assay(SE, assayname)   
  
  # Select specific features if provided  
  if (!is.null(feature_of_interesting)) {  
    expdata <- expdata[rownames(expdata) %in% feature_of_interesting, ]  
  }  

  row_variances <- apply(expdata, 1, var, na.rm = TRUE)  
  constant_rows <- which(row_variances == 0)  

  if (length(constant_rows) > 0) {  
    print(paste0("Delete ", length(constant_rows), " common express features"))  
    expdata <- expdata[-constant_rows, ]  
  } else {  
    print("No common express features to delete.")  
  }  
  num_feature <- nrow(expdata)   

  sample_info = colData(SE)  
  pca_result <- prcomp(as.data.frame(t(expdata)), scale. = scale)   
  pca_data <- as.data.frame(pca_result$x)  

  missing_samples <- setdiff(rownames(pca_data), rownames(sample_info))  
  if (length(missing_samples) > 0) {  
    print(paste("Warning: Missing sample information for:", paste(missing_samples, collapse = ", ")))  
    pca_data <- pca_data[!(rownames(pca_data) %in% missing_samples), ]  
  }  

  pca_data$group <- mapvalues(rownames(pca_data), rownames(sample_info), sample_info[, groupname], warn_missing = FALSE)  
  
  mean_pc1 <- mean(pca_data$PC1)  
  sd_pc1 <- sd(pca_data$PC1)  
  mean_pc2 <- mean(pca_data$PC2)  
  sd_pc2 <- sd(pca_data$PC2)  
  pca_data_filter <- pca_data[(pca_data$PC1 > (mean_pc1 - outlier_threshold * sd_pc1)) & (pca_data$PC1 < (mean_pc1 + outlier_threshold * sd_pc1)) &  
                               (pca_data$PC2 > (mean_pc2 - outlier_threshold * sd_pc2)) & (pca_data$PC2 < (mean_pc2 + outlier_threshold * sd_pc2)),]  

  group_counts <- table(pca_data$group)  
  caption_text1 <- paste("feature:", num_feature, ";", "sample:", paste(names(group_counts), ":", group_counts, collapse = ", "))  

  group_counts_filtered <- table(pca_data_filter$group)  
  caption_text2 <- paste("feature:", num_feature, ";", "sample:", paste(names(group_counts_filtered), ":", group_counts_filtered, collapse = ", "))  

  plot_pca <- function(data, title, caption, show_legend = TRUE) {  
    p <- ggplot(data, aes(x = PC1, y = PC2, color = group)) +  
      geom_point(size = 3) +  
      labs(title = paste0("PCAplot: ", SCvalue(data)), x = "PCA1", y = "PCA2") +  
      theme_minimal() +  
      stat_ellipse(type = "norm", level = 0.95, linetype = 2) +  
      theme(axis.title.x = element_text(size = 12),  
            axis.title.y = element_text(size = 12),  
            axis.text = element_text(size = 12),  
            plot.title = element_text(size = 12),  
            panel.grid.major = element_blank(),  
            panel.grid.minor = element_blank(),  
            strip.text = element_text(size = 12),  
            panel.border = element_rect(color = "gray", fill = NA, size = 1),  
            legend.position = ifelse(show_legend, "right", "none"),  # 将图例位置设置为右侧  
            legend.text = element_text(size = 12),  
            legend.title = element_text(size = 12))  

    if (show_caption) {  
      p <- p + labs(caption = caption)  
    }  
    return(p)  
	}  

  pca_plot1 <- plot_pca(pca_data, "PCAplot", caption_text1, show_legend = FALSE)     # 不显示图例  
  pca_plot2 <- plot_pca(pca_data_filter, "PCAplot Filtered", caption_text2, show_legend = TRUE)  # 显示图例，位置设置为右侧  

  return(plot_grid(pca_plot1, pca_plot2, align = 'hv', ncol = 2))
}