#' @title SE_MAplot: Create MA Plot for Differential Expression Results
#' @description This function creates an MA plot to visualize differential expression results from a SummarizedExperiment object. It shows the relationship between average expression and log2 fold change.
#' @param SE A \code{SummarizedExperiment} object containing differential expression results in rowData.
#' @param assayname A string indicating which assay to use for average expression calculation. Default is "log2".
#' @param logFC_col A string indicating the column name in rowData for log2 fold change. Default is "logFC".
#' @param pvalue_col A string indicating the column name in rowData for p-values. Default is "adj.P.Val".
#' @param logFC_threshold Numeric value for log2 fold change threshold to highlight. Default is 1.
#' @param pvalue_threshold Numeric value for p-value threshold to highlight. Default is 0.05.
#' @param title A character string for the plot title. Default is "MA Plot".
#' @param label_genes Numeric value indicating the number of top genes to label. Default is 10.
#' @return A \code{ggplot} object representing the MA plot.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform differential expression analysis
#' SE_limma_result <- SE_limma(SE, assayname = "log2", group_colname = "group")
#' 
#' # Create MA plot
#' ma_plot <- SE_MAplot(SE_limma_result, assayname = "log2", logFC_col = "logFC", 
#'                      pvalue_col = "adj.P.Val")
#' print(ma_plot)
#' @export
SE_MAplot <- function(SE, assayname = "log2", logFC_col = "logFC", pvalue_col = "adj.P.Val", 
                      logFC_threshold = 1, pvalue_threshold = 0.05, 
                      title = "MA Plot", label_genes = 10) {
    
    exp_data <- assay(SE, assayname)
    results <- as.data.frame(rowData(SE))
    
    results$A <- rowMeans(exp_data, na.rm = TRUE)
    results$M <- results[[logFC_col]]
    results$significant <- (abs(results[[logFC_col]]) >= logFC_threshold) & 
                           (results[[pvalue_col]] < pvalue_threshold)
    results$direction <- ifelse(results$M > 0, "Up", "Down")
    results$direction[!results$significant] <- "Not Significant"
    
    ma_plot <- ggplot(results, aes(x = A, y = M)) +
        geom_point(aes(color = direction), alpha = 0.6, size = 1.5) +
        scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Significant" = "grey")) +
        geom_hline(yintercept = 0, linetype = "solid", color = "darkgrey") +
        geom_hline(yintercept = c(-logFC_threshold, logFC_threshold), 
                   linetype = "dashed", color = "darkgrey") +
        labs(title = title,
             x = "Average Expression (A)",
             y = "Log2 Fold Change (M)",
             color = "Significance") +
        theme_minimal() +
        theme(legend.position = "right",
              plot.title = element_text(hjust = 0.5))
    
    if (label_genes > 0 && sum(results$significant) > 0) {
        top_genes <- results[results$significant, ]
        top_genes <- top_genes[order(top_genes[[pvalue_col]]), ]
        top_genes <- top_genes[1:min(label_genes, nrow(top_genes)), ]
        
        ma_plot <- ma_plot + 
            geom_text(data = top_genes, 
                     aes(label = rownames(top_genes)),
                     vjust = -0.5, hjust = 0.5, size = 2, color = "black")
    }
    
    cat("MA plot created\n")
    cat("Up-regulated genes:", sum(results$significant & results$M > 0), "\n")
    cat("Down-regulated genes:", sum(results$significant & results$M < 0), "\n")
    
    return(ma_plot)
}
