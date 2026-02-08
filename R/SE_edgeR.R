#' @title SE_edgeR: Differential Expression Analysis Using edgeR
#' @description This function performs differential expression analysis using the edgeR package, which is specifically designed for RNA-seq count data. It uses negative binomial generalized linear models to identify differentially expressed features.
#' @param SE A \code{SummarizedExperiment} object containing count data.
#' @param assayname A string indicating which assay to use for analysis. The default value is \code{"Counts"}.
#' @param group_colname A string representing the column name in \code{colData} that contains group information. Default is "group".
#' @param contrast A character string specifying the contrast for differential analysis (e.g., "Treatment-Control"). Default is NULL, which will compare the first two groups.
#' @param adjust_method A character string specifying the p-value adjustment method. Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default is "BH".
#' @param pvalue_threshold Numeric value for p-value threshold. Default is 0.05.
#' @param logFC_threshold Numeric value for log2 fold change threshold. Default is 1.
#' @param normalize_method A character string specifying the normalization method. Options include "TMM", "RLE", "upperquartile", "none". Default is "TMM".
#' @return A \code{SummarizedExperiment} object with differential expression results added to rowData, including columns: logFC, logCPM, PValue, FDR, and significant (logical).
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform differential expression analysis using edgeR
#' SE_edgeR_result <- SE_edgeR(SE, assayname = "Counts", group_colname = "group", 
#'                             contrast = "Treatment-Control")
#' 
#' # Extract significant genes
#' sig_genes <- rowData(SE_edgeR_result)$significant
#' @export
SE_edgeR <- function(SE, assayname = "Counts", group_colname = "group", contrast = NULL, 
                     adjust_method = "BH", pvalue_threshold = 0.05, logFC_threshold = 1, 
                     normalize_method = "TMM") {
    
    exp_data <- assay(SE, assayname)
    groups <- colData(SE)[[group_colname]]
    
    dge <- DGEList(counts = exp_data, group = groups)
    
    dge <- calcNormFactors(dge, method = normalize_method)
    
    design <- model.matrix(~0 + groups)
    colnames(design) <- levels(factor(groups))
    
    dge <- estimateDisp(dge, design)
    
    fit <- glmFit(dge, design)
    
    if (is.null(contrast)) {
        contrast <- paste(colnames(design)[2], colnames(design)[1], sep = "-")
    }
    
    contrast_matrix <- makeContrasts(contrasts = contrast, levels = design)
    lrt <- glmLRT(fit, contrast = contrast_matrix)
    
    results <- topTags(lrt, n = Inf, adjust.method = adjust_method)
    results_df <- as.data.frame(results)
    
    results_df$significant <- (results_df$FDR < pvalue_threshold) & (abs(results_df$logFC) >= logFC_threshold)
    
    # 格式化结果表格
    results_df <- format_result_table(results_df, pvalue_cols = c("PValue", "FDR"))
    
    rowData(SE)$logFC <- results_df$logFC[rownames(SE)]
    rowData(SE)$logCPM <- results_df$logCPM[rownames(SE)]
    rowData(SE)$PValue <- results_df$PValue[rownames(SE)]
    rowData(SE)$FDR <- results_df$FDR[rownames(SE)]
    rowData(SE)$significant <- results_df$significant[rownames(SE)]
    
    cat("edgeR differential expression analysis completed\n")
    cat("Contrast:", contrast, "\n")
    cat("Normalization method:", normalize_method, "\n")
    cat("Significant features:", sum(results_df$significant, na.rm = TRUE), "\n")
    
    return(SE)
}
