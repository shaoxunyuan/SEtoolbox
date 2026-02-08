#' @title SE_KEGG: KEGG Pathway Enrichment Analysis
#' @description This function performs KEGG pathway enrichment analysis for differentially expressed genes from a SummarizedExperiment object. It identifies over-represented KEGG pathways.
#' @param SE A \code{SummarizedExperiment} object containing differential expression results.
#' @param gene_col A string indicating the column name in rowData containing gene identifiers. Default is "gene_id".
#' @param pvalue_threshold Numeric value for p-value threshold to select significant genes. Default is 0.05.
#' @param logFC_threshold Numeric value for log2 fold change threshold to select significant genes. Default is 1.
#' @param organism A character string specifying the organism code (e.g., "hsa" for human, "mmu" for mouse). Default is "hsa".
#' @param p_adjust_method A character string specifying the p-value adjustment method. Default is "BH".
#' @param qvalue_threshold Numeric value for adjusted p-value threshold. Default is 0.05.
#' @param universe A character vector of gene IDs representing the background universe. If NULL, all genes in SE will be used.
#' @return A data frame containing KEGG enrichment results with columns: ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, and count.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform differential expression analysis
#' SE_limma_result <- SE_limma(SE, assayname = "log2", group_colname = "group")
#' 
#' # Perform KEGG enrichment analysis
#' kegg_results <- SE_KEGG(SE_limma_result, organism = "hsa")
#' print(kegg_results)
#' @export
SE_KEGG <- function(SE, gene_col = "gene_id", pvalue_threshold = 0.05, logFC_threshold = 1, 
                     organism = "hsa", p_adjust_method = "BH", qvalue_threshold = 0.05, universe = NULL) {
    
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
    
    cat("KEGG pathway enrichment analysis initiated\n")
    cat("Number of significant genes:", length(sig_genes), "\n")
    cat("Number of genes in universe:", length(universe), "\n")
    
    gene_entrez <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", 
                        OrgDb = org.Hs.eg.db)
    
    if (nrow(gene_entrez) == 0) {
        stop("No gene IDs could be converted to Entrez IDs")
    }
    
    universe_entrez <- bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", 
                            OrgDb = org.Hs.eg.db)
    
    ekegg <- enrichKEGG(gene = gene_entrez$ENTREZID,
                        universe = universe_entrez$ENTREZID,
                        organism = organism,
                        pAdjustMethod = p_adjust_method,
                        pvalueCutoff = qvalue_threshold,
                        qvalueCutoff = qvalue_threshold)
    
    if (is.null(ekegg) || nrow(ekegg) == 0) {
        cat("No significant KEGG pathways found\n")
        return(data.frame())
    }
    
    kegg_results <- as.data.frame(ekegg)
    
    # 格式化结果表格
    kegg_results <- format_result_table(kegg_results, pvalue_cols = c("pvalue", "p.adjust", "qvalue"))
    
    cat("KEGG pathway enrichment analysis completed\n")
    cat("Number of enriched KEGG pathways:", nrow(kegg_results), "\n")
    
    return(kegg_results)
}
