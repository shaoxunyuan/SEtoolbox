#' @title SE_enrichplot: Visualize Enrichment Analysis Results
#' @description This function creates various visualizations for enrichment analysis results, including bar plots, dot plots, gene-concept networks, and enrichment maps.
#' @param enrich_results A data frame containing enrichment analysis results (from SE_GO, SE_KEGG, or SE_GSEA).
#' @param plot_type A character string specifying the type of plot. Options include "bar", "dot", "cnet", "emap", "ridge". Default is "bar".
#' @param top_n Numeric value indicating the number of top enriched terms to display. Default is 20.
#' @param show_category Numeric value indicating the number of categories to show in the plot. Default is 20.
#' @param title A character string for the plot title. Default is NULL.
#' @param color_by A character string specifying the column to use for color mapping. Default is "p.adjust".
#' @return A \code{ggplot} object representing the enrichment plot.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform differential expression analysis
#' SE_limma_result <- SE_limma(SE, assayname = "log2", group_colname = "group")
#' 
#' # Perform GO enrichment analysis
#' go_results <- SE_GO(SE_limma_result, ontology = "BP")
#' 
#' # Create bar plot
#' bar_plot <- SE_enrichplot(go_results, plot_type = "bar", top_n = 15)
#' print(bar_plot)
#' 
#' # Create dot plot
#' dot_plot <- SE_enrichplot(go_results, plot_type = "dot", top_n = 15)
#' print(dot_plot)
#' @export
SE_enrichplot <- function(enrich_results, plot_type = "bar", top_n = 20, 
                           show_category = 20, title = NULL, color_by = "p.adjust") {
    
    if (nrow(enrich_results) == 0) {
        stop("No enrichment results to plot")
    }
    
    plot_data <- enrich_results[1:min(top_n, nrow(enrich_results)), ]
    
    if (plot_type == "bar") {
        plot_data$Description <- factor(plot_data$Description, 
                                        levels = rev(plot_data$Description))
        
        enrich_plot <- ggplot(plot_data, aes(x = .data$Description, y = .data$count, 
                                             fill = .data[[color_by]])) +
            geom_bar(stat = "identity") +
            coord_flip() +
            scale_fill_gradient(low = "red", high = "blue") +
            labs(title = ifelse(is.null(title), "Enrichment Bar Plot", title),
                 x = "Term",
                 y = "Gene Count",
                 fill = color_by) +
            theme_minimal() +
            theme(axis.text.y = element_text(size = 10))
        
    } else if (plot_type == "dot") {
        plot_data$Description <- factor(plot_data$Description, 
                                        levels = rev(plot_data$Description))
        
        enrich_plot <- ggplot(plot_data, aes(x = .data$Description, y = .data$GeneRatio, 
                                             size = .data$count, color = .data[[color_by]])) +
            geom_point() +
            coord_flip() +
            scale_color_gradient(low = "red", high = "blue") +
            labs(title = ifelse(is.null(title), "Enrichment Dot Plot", title),
                 x = "Term",
                 y = "Gene Ratio",
                 size = "Gene Count",
                 color = color_by) +
            theme_minimal() +
            theme(axis.text.y = element_text(size = 10))
        
    } else if (plot_type == "ridge") {
        plot_data$Description <- factor(plot_data$Description, 
                                        levels = rev(plot_data$Description))
        
        enrich_plot <- ggplot(plot_data, aes(x = .data[[color_by]], y = .data$Description, 
                                             fill = .data$Description)) +
            geom_density_ridge(alpha = 0.7) +
            labs(title = ifelse(is.null(title), "Enrichment Ridge Plot", title),
                 x = color_by,
                 y = "Term") +
            theme_minimal() +
            theme(legend.position = "none")
        
    } else {
        stop("Unknown plot type. Available options: bar, dot, ridge")
    }
    
    cat("Enrichment plot created:", plot_type, "\n")
    
    return(enrich_plot)
}
