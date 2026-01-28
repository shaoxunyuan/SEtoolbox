#' @title SE_QCreport: Generate Quality Control Report
#' @description This function generates a comprehensive quality control report for a SummarizedExperiment object. It includes various QC metrics and visualizations to assess data quality.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for QC analysis. The default value is \code{"TPM"}.
#' @param output_dir A string indicating the directory to save QC report. Default is NULL (no saving).
#' @param save_plots Logical value indicating whether to save plots. Default is FALSE.
#' @return A list containing:  
#' \item{summary}{A data frame with QC summary statistics.}  
#' \item{plots}{A list of ggplot objects with various QC visualizations.}
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Generate QC report
#' qc_report <- SE_QCreport(SE, assayname = "TPM")
#' 
#' # View QC summary
#' print(qc_report$summary)
#' 
#' # View QC plots
#' print(qc_report$plots$boxplot)
#' print(qc_report$plots$density)
#' print(qc_report$plots$pca)
#' @export
SE_QCreport <- function(SE, assayname = "TPM", output_dir = NULL, save_plots = FALSE) {
    
    exp_data <- assay(SE, assayname)
    
    summary_df <- data.frame(
        Metric = c("Number of features", "Number of samples", "Total reads", 
                   "Mean expression", "Median expression", "SD expression",
                   "Missing values", "Zero values"),
        Value = c(nrow(exp_data), ncol(exp_data), sum(exp_data, na.rm = TRUE),
                  mean(exp_data, na.rm = TRUE), median(exp_data, na.rm = TRUE),
                  sd(exp_data, na.rm = TRUE),
                  sum(is.na(exp_data)), sum(exp_data == 0, na.rm = TRUE)),
        stringsAsFactors = FALSE
    )
    
    sample_sums <- colSums(exp_data, na.rm = TRUE)
    sample_means <- colMeans(exp_data, na.rm = TRUE)
    sample_zeros <- colSums(exp_data == 0, na.rm = TRUE)
    
    sample_qc_df <- data.frame(
        Sample = colnames(exp_data),
        Total = sample_sums,
        Mean = sample_means,
        Zeros = sample_zeros,
        stringsAsFactors = FALSE
    )
    
    boxplot_df <- data.frame(
        Sample = rep(colnames(exp_data), each = nrow(exp_data)),
        Expression = as.vector(exp_data)
    )
    
    boxplot <- ggplot(boxplot_df, aes(x = Sample, y = Expression)) +
        geom_boxplot(fill = "steelblue", alpha = 0.7) +
        labs(title = "Sample Expression Boxplot",
             x = "Sample",
             y = "Expression") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    density_df <- data.frame(
        Sample = rep(colnames(exp_data), each = nrow(exp_data)),
        Expression = as.vector(exp_data)
    )
    
    density_plot <- ggplot(density_df, aes(x = Expression, color = Sample)) +
        geom_density(alpha = 0.5) +
        labs(title = "Sample Expression Density",
             x = "Expression",
             y = "Density") +
        theme_minimal()
    
    variances <- apply(exp_data, 1, var, na.rm = TRUE)
    top_features <- names(sort(variances, decreasing = TRUE))[1:min(1000, length(variances))]
    exp_data_subset <- exp_data[top_features, ]
    
    pca_result <- prcomp(t(exp_data_subset), scale. = TRUE, center = TRUE)
    
    pca_df <- data.frame(
        PC1 = pca_result$x[, 1],
        PC2 = pca_result$x[, 2],
        Sample = rownames(pca_result$x)
    )
    
    pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
        geom_point(size = 3, alpha = 0.7) +
        geom_text(aes(label = Sample), vjust = -0.5, hjust = 0.5, size = 2) +
        labs(title = "PCA of Samples",
             x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
             y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)")) +
        theme_minimal()
    
    plots <- list(
        boxplot = boxplot,
        density = density_plot,
        pca = pca_plot
    )
    
    cat("Quality control report generated\n")
    cat("Number of features:", nrow(exp_data), "\n")
    cat("Number of samples:", ncol(exp_data), "\n")
    
    if (!is.null(output_dir) && save_plots) {
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        
        ggsave(file.path(output_dir, "qc_boxplot.png"), plots$boxplot, width = 10, height = 6)
        ggsave(file.path(output_dir, "qc_density.png"), plots$density, width = 10, height = 6)
        ggsave(file.path(output_dir, "qc_pca.png"), plots$pca, width = 8, height = 8)
        
        write.csv(summary_df, file.path(output_dir, "qc_summary.csv"), row.names = FALSE)
        write.csv(sample_qc_df, file.path(output_dir, "sample_qc.csv"), row.names = FALSE)
        
        cat("QC report saved to:", output_dir, "\n")
    }
    
    return(list(
        summary = summary_df,
        sample_qc = sample_qc_df,
        plots = plots
    ))
}
