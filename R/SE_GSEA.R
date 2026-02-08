#' @title SE_GSEA: Gene Set Enrichment Analysis
#' @description This function performs Gene Set Enrichment Analysis (GSEA) using ranked gene lists from a SummarizedExperiment object. It identifies pathways or gene sets that are coordinately up- or down-regulated.
#' @param SE A \code{SummarizedExperiment} object containing expression data.
#' @param assayname A string indicating which assay to use for ranking. The default value is \code{"log2"}.
#' @param group_colname A string representing the column name in \code{colData} that contains group information. Default is "group".
#' @param gene_sets A list of gene sets for enrichment analysis. Each element should be a character vector of gene names.
#' @param ranking_metric A character string specifying the metric for ranking genes. Options include "logFC", "t", "signal2noise". Default is "logFC".
#' @param nperm Numeric value indicating the number of permutations for significance testing. Default is 1000.
#' @param min_size Numeric value indicating the minimum size of gene sets to consider. Default is 10.
#' @param max_size Numeric value indicating the maximum size of gene sets to consider. Default is 500.
#' @return A data frame containing GSEA results with columns: pathway, NES, NOM p-value, FDR q-value, and leading_edge_genes.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Define gene sets (example)
#' gene_sets <- list(
#'     Pathway1 = c("Gene1", "Gene2", "Gene3"),
#'     Pathway2 = c("Gene4", "Gene5", "Gene6")
#' )
#' 
#' # Perform GSEA
#' gsea_results <- SE_GSEA(SE, assayname = "log2", group_colname = "group", 
#'                          gene_sets = gene_sets)
#' print(gsea_results)
#' @export
SE_GSEA <- function(SE, assayname = "log2", group_colname = "group", gene_sets, 
                    ranking_metric = "logFC", nperm = 1000, min_size = 10, max_size = 500) {
    
    exp_data <- assay(SE, assayname)
    groups <- colData(SE)[[group_colname]]
    
    group_levels <- unique(groups)
    if (length(group_levels) != 2) {
        stop("GSEA requires exactly two groups for comparison")
    }
    
    group1 <- groups == group_levels[1]
    group2 <- groups == group_levels[2]
    
    if (ranking_metric == "logFC") {
        mean1 <- rowMeans(exp_data[, group1, drop = FALSE], na.rm = TRUE)
        mean2 <- rowMeans(exp_data[, group2, drop = FALSE], na.rm = TRUE)
        ranking <- mean2 - mean1
    } else if (ranking_metric == "t") {
        mean1 <- rowMeans(exp_data[, group1, drop = FALSE], na.rm = TRUE)
        mean2 <- rowMeans(exp_data[, group2, drop = FALSE], na.rm = TRUE)
        sd1 <- apply(exp_data[, group1, drop = FALSE], 1, sd, na.rm = TRUE)
        sd2 <- apply(exp_data[, group2, drop = FALSE], 1, sd, na.rm = TRUE)
        n1 <- sum(group1)
        n2 <- sum(group2)
        pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
        ranking <- (mean2 - mean1) / (pooled_sd * sqrt(1/n1 + 1/n2))
    } else if (ranking_metric == "signal2noise") {
        mean1 <- rowMeans(exp_data[, group1, drop = FALSE], na.rm = TRUE)
        mean2 <- rowMeans(exp_data[, group2, drop = FALSE], na.rm = TRUE)
        sd1 <- apply(exp_data[, group1, drop = FALSE], 1, sd, na.rm = TRUE)
        sd2 <- apply(exp_data[, group2, drop = FALSE], 1, sd, na.rm = TRUE)
        ranking <- (mean2 - mean1) / ((sd1 + sd2) / 2)
    } else {
        stop("Unknown ranking metric. Available options: logFC, t, signal2noise")
    }
    
    ranking <- sort(ranking, decreasing = TRUE)
    
    gsea_results <- data.frame(
        pathway = character(),
        NES = numeric(),
        NOM_pvalue = numeric(),
        FDR_qvalue = numeric(),
        leading_edge_genes = character(),
        stringsAsFactors = FALSE
    )
    
    for (pathway_name in names(gene_sets)) {
        genes_in_pathway <- gene_sets[[pathway_name]]
        genes_in_data <- intersect(genes_in_pathway, names(ranking))
        
        if (length(genes_in_data) < min_size || length(genes_in_data) > max_size) {
            next
        }
        
        pathway_ranks <- ranking[genes_in_data]
        
        running_es <- cumsum(pathway_ranks)
        max_es <- max(running_es)
        min_es <- min(running_es)
        
        es <- ifelse(abs(max_es) > abs(min_es), max_es, min_es)
        
        null_es <- replicate(nperm, {
            random_genes <- sample(names(ranking), length(genes_in_data))
            random_ranks <- ranking[random_genes]
            random_running_es <- cumsum(random_ranks)
            random_max_es <- max(random_running_es)
            random_min_es <- min(random_running_es)
            ifelse(abs(random_max_es) > abs(random_min_es), random_max_es, random_min_es)
        })
        
        if (es > 0) {
            nom_pvalue <- sum(null_es >= es) / nperm
        } else {
            nom_pvalue <- sum(null_es <= es) / nperm
        }
        
        nes <- es / mean(abs(null_es))
        
        leading_edge_idx <- which.max(abs(running_es))
        leading_edge_genes <- paste(genes_in_data[1:leading_edge_idx], collapse = ", ")
        
        gsea_results <- rbind(gsea_results, data.frame(
            pathway = pathway_name,
            NES = nes,
            NOM_pvalue = nom_pvalue,
            FDR_qvalue = p.adjust(gsea_results$NOM_pvalue, method = "BH"),
            leading_edge_genes = leading_edge_genes,
            stringsAsFactors = FALSE
        ))
    }
    
    gsea_results$FDR_qvalue <- p.adjust(gsea_results$NOM_pvalue, method = "BH")
    gsea_results <- gsea_results[order(gsea_results$NES, decreasing = TRUE), ]
    rownames(gsea_results) <- NULL
    
    # 格式化结果表格
    gsea_results <- format_result_table(gsea_results, pvalue_cols = c("NOM_pvalue", "FDR_qvalue"))
    
    cat("GSEA analysis completed\n")
    cat("Number of gene sets analyzed:", nrow(gsea_results), "\n")
    
    return(gsea_results)
}
