#' @title SE_GO: GO Enrichment Analysis
#' @description This function performs Gene Ontology (GO) enrichment analysis for differentially expressed genes from a SummarizedExperiment object. It identifies over-represented GO terms in Biological Process, Molecular Function, and Cellular Component categories.
#' @param SE A \code{SummarizedExperiment} object containing differential expression results.
#' @param gene_col A string indicating the column name in rowData containing gene identifiers. Default is "gene_id".
#' @param pvalue_threshold Numeric value for p-value threshold to select significant genes. Default is 0.05.
#' @param logFC_threshold Numeric value for log2 fold change threshold to select significant genes. Default is 1.
#' @param ontology A character string specifying the GO ontology to analyze. Options include "BP", "MF", "CC", or "ALL". Default is "ALL".
#' @param p_adjust_method A character string specifying the p-value adjustment method. Default is "BH".
#' @param qvalue_threshold Numeric value for adjusted p-value threshold. Default is 0.05.
#' @param universe A character vector of gene IDs representing the background universe. If NULL, all genes in SE will be used.
#' @return A data frame containing GO enrichment results with columns: GO_ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, and count.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform differential expression analysis
#' SE_limma_result <- SE_limma(SE, assayname = "log2", group_colname = "group")
#' 
#' # Perform GO enrichment analysis
#' go_results <- SE_GO(SE_limma_result, ontology = "BP")
#' print(go_results)
#' @export
SE_GO <- function(SE, gene_col = "gene_id", pvalue_threshold = 0.05, logFC_threshold = 1, 
                  ontology = "ALL", p_adjust_method = "BH", qvalue_threshold = 0.05, universe = NULL) {
    
    if (!"significant" %in% colnames(rowData(SE))) {
        stop("SE object must contain a 'significant' column in rowData. Run differential expression analysis first.")
    }
    
    sig_genes <- rownames(SE)[rowData(SE)$significant & 
                              !is.na(rowData(SE)$significant)]
    
    if (length(sig_genes) == 0) {
        stop("No significant genes found. Check pvalue_threshold and logFC_threshold parameters.")
    }
    
    if (is.null(universe)) {
        universe <- rownames(SE)
    }
    
    cat("GO enrichment analysis initiated\n")
    cat("Number of significant genes:", length(sig_genes), "\n")
    cat("Number of genes in universe:", length(universe), "\n")
    
    ego <- enrichGO(gene = sig_genes,
                    universe = universe,
                    OrgDb = org.Hs.eg.db,
                    keyType = "SYMBOL",
                    ont = ontology,
                    pAdjustMethod = p_adjust_method,
                    pvalueCutoff = qvalue_threshold,
                    qvalueCutoff = qvalue_threshold,
                    readable = TRUE)
    
    if (is.null(ego) || nrow(ego) == 0) {
        cat("No significant GO terms found\n")
        return(data.frame())
    }
    
    go_results <- as.data.frame(ego)
    
    cat("GO enrichment analysis completed\n")
    cat("Number of enriched GO terms:", nrow(go_results), "\n")
    
    return(go_results)
}
