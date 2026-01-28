#' @title SE_limma: Differential Expression Analysis Using limma
#' @description This function performs differential expression analysis using the limma package, which is suitable for both microarray and RNA-seq data. It uses linear models and empirical Bayes methods to identify differentially expressed features.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for analysis. The default value is \code{"log2"}.
#' @param group_colname A string representing the column name in \code{colData} that contains group information. Default is "group".
#' @param contrast A character string specifying the contrast for differential analysis (e.g., "Treatment-Control"). Default is NULL, which will compare the first two groups.
#' @param adjust_method A character string specifying the p-value adjustment method. Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default is "BH".
#' @param pvalue_threshold Numeric value for p-value threshold. Default is 0.05.
#' @param logFC_threshold Numeric value for log2 fold change threshold. Default is 1.
#' @return A \code{SummarizedExperiment} object with differential expression results added to rowData, including columns: logFC, AveExpr, t, P.Value, adj.P.Val, B, and significant (logical).
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform differential expression analysis using limma
#' SE_limma_result <- SE_limma(SE, assayname = "log2", group_colname = "group", 
#'                              contrast = "Treatment-Control")
#' 
#' # Extract significant genes
#' sig_genes <- rowData(SE_limma_result)$significant
#' @export
SE_limma <- function(SE, assayname = "log2", group_colname = "group", contrast = NULL, 
                     adjust_method = "BH", pvalue_threshold = 0.05, logFC_threshold = 1) {
    
    exp_data <- assay(SE, assayname)
    groups <- colData(SE)[[group_colname]]
    
    design <- model.matrix(~0 + groups)
    colnames(design) <- levels(factor(groups))
    
    fit <- lmFit(exp_data, design)
    
    if (is.null(contrast)) {
        contrast <- paste(colnames(design)[2], colnames(design)[1], sep = "-")
    }
    
    contrast_matrix <- makeContrasts(contrasts = contrast, levels = design)
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    results <- topTable(fit2, number = Inf, adjust.method = adjust_method)
    
    results$significant <- (results$adj.P.Val < pvalue_threshold) & (abs(results$logFC) >= logFC_threshold)
    
    rowData(SE)$logFC <- results$logFC[rownames(SE)]
    rowData(SE)$AveExpr <- results$AveExpr[rownames(SE)]
    rowData(SE)$t <- results$t[rownames(SE)]
    rowData(SE)$P.Value <- results$P.Value[rownames(SE)]
    rowData(SE)$adj.P.Val <- results$adj.P.Val[rownames(SE)]
    rowData(SE)$B <- results$B[rownames(SE)]
    rowData(SE)$significant <- results$significant[rownames(SE)]
    
    cat("limma differential expression analysis completed\n")
    cat("Contrast:", contrast, "\n")
    cat("Significant features:", sum(results$significant, na.rm = TRUE), "\n")
    
    return(SE)
}
