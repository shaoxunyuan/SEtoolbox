#' @title SE_volcano: Create Volcano Plot for Differential Expression Results
#' @description This function creates a volcano plot to visualize differential expression results from a SummarizedExperiment object. It shows the relationship between log2 fold change and statistical significance.
#' @param SE A \code{SummarizedExperiment} object containing differential expression results in rowData.
#' @param logFC_col A string indicating the column name in rowData for log2 fold change. Default is "logFC".
#' @param pvalue_col A string indicating the column name in rowData for p-values. Default is "adj.P.Val".
#' @param logFC_threshold Numeric value for log2 fold change threshold to highlight. Default is 1.
#' @param pvalue_threshold Numeric value for p-value threshold to highlight. Default is 0.05.
#' @param title A character string for the plot title. Default is "Volcano Plot".
#' @param label_genes Numeric value indicating the number of top genes to label. Default is 10.
#' @return A \code{ggplot} object representing the volcano plot.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform differential expression analysis
#' SE_limma_result <- SE_limma(SE, assayname = "log2", group_colname = "group")
#' 
#' # Create volcano plot
#' volcano_plot <- SE_volcano(SE_limma_result, logFC_col = "logFC", pvalue_col = "adj.P.Val")
#' print(volcano_plot)
#' @export
SE_volcano <- function(SE, logFC_col = "logFC", pvalue_col = "adj.P.Val", 
                       logFC_threshold = 1, pvalue_threshold = 0.05, 
                       title = "Volcano Plot", label_genes = 10) {
    
    results <- as.data.frame(rowData(SE))
    results$neg_log10_pvalue <- -log10(results[[pvalue_col]])
    results$significant <- (abs(results[[logFC_col]]) >= logFC_threshold) & 
                           (results[[pvalue_col]] < pvalue_threshold)
    results$direction <- ifelse(results[[logFC_col]] > 0, "Up", "Down")
    results$direction[!results$significant] <- "Not Significant"
    
    volcano_plot <- ggplot(results, aes(x = .data[[logFC_col]], y = neg_log10_pvalue)) +
        geom_point(aes(color = direction), alpha = 0.6, size = 1.5) +
        scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Significant" = "grey")) +
        geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), 
                   linetype = "dashed", color = "darkgrey") +
        geom_hline(yintercept = -log10(pvalue_threshold), 
                   linetype = "dashed", color = "darkgrey") +
        labs(title = title,
             x = "Log2 Fold Change",
             y = "-Log10 P-value",
             color = "Significance") +
        theme_minimal() +
        theme(legend.position = "right",
              plot.title = element_text(hjust = 0.5))
    
    if (label_genes > 0 && sum(results$significant) > 0) {
        top_genes <- results[results$significant, ]
        top_genes <- top_genes[order(top_genes[[pvalue_col]]), ]
        top_genes <- top_genes[1:min(label_genes, nrow(top_genes)), ]
        
        volcano_plot <- volcano_plot + 
            geom_text(data = top_genes, 
                     aes(label = rownames(top_genes)),
                     vjust = -0.5, hjust = 0.5, size = 2, color = "black")
    }
    
    cat("Volcano plot created\n")
    cat("Up-regulated genes:", sum(results$significant & results[[logFC_col]] > 0), "\n")
    cat("Down-regulated genes:", sum(results$significant & results[[logFC_col]] < 0), "\n")
    
    return(volcano_plot)
}
