#' @title SE_networkplot: Create Network Plot from Correlation Matrix
#' @description This function creates a network plot from a correlation matrix or adjacency matrix. It visualizes relationships between features or samples as a network graph.
#' @param cor_matrix A correlation matrix or adjacency matrix.
#' @param threshold Numeric value indicating the minimum absolute correlation to include in the network. Default is 0.5.
#' @param max_nodes Numeric value indicating the maximum number of nodes to display in the network. Default is 50.
#' @param layout A character string specifying the layout algorithm. Options include "fr", "kk", "grid", "circle", "random". Default is "fr".
#' @param node_size Numeric value for node size. Default is 10.
#' @param edge_width Numeric value for edge width. Default is 1.
#' @param label_nodes Logical value indicating whether to label nodes. Default is FALSE.
#' @return A \code{ggplot} object representing the network.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Compute feature-feature correlations
#' cor_matrix <- SE_correlation(SE, assayname = "log2", method = "features")
#' 
#' # Create network plot
#' network_plot <- SE_networkplot(cor_matrix, threshold = 0.7)
#' print(network_plot)
#' @export
SE_networkplot <- function(cor_matrix, threshold = 0.5, max_nodes = 50, 
                          layout = "fr", node_size = 10, edge_width = 1, 
                          label_nodes = FALSE) {
    
    if (threshold < 0 || threshold > 1) {
        stop("threshold must be between 0 and 1")
    }
    
    adj_matrix <- cor_matrix
    adj_matrix[abs(adj_matrix) < threshold] <- 0
    diag(adj_matrix) <- 0
    
    node_degrees <- rowSums(abs(adj_matrix))
    top_nodes <- names(sort(node_degrees, decreasing = TRUE))[1:min(max_nodes, length(node_degrees))]
    
    adj_matrix_subset <- adj_matrix[top_nodes, top_nodes]
    
    g <- igraph::graph_from_adjacency_matrix(adj_matrix_subset, mode = "undirected", 
                                             weighted = TRUE, diag = FALSE)
    
    if (layout == "fr") {
        layout_coords <- igraph::layout_with_fr(g)
    } else if (layout == "kk") {
        layout_coords <- igraph::layout_with_kk(g)
    } else if (layout == "grid") {
        layout_coords <- igraph::layout_on_grid(g)
    } else if (layout == "circle") {
        layout_coords <- igraph::layout_in_circle(g)
    } else if (layout == "random") {
        layout_coords <- igraph::layout_randomly(g)
    } else {
        layout_coords <- igraph::layout_with_fr(g)
    }
    
    nodes_df <- data.frame(
        id = igraph::V(g)$name,
        x = layout_coords[, 1],
        y = layout_coords[, 2],
        degree = igraph::degree(g)
    )
    
    edges_df <- igraph::as_data_frame(g, what = "edges")
    colnames(edges_df) <- c("from", "to", "weight")
    
    network_plot <- ggplot() +
        geom_segment(data = edges_df, aes(x = .data$x[match(.data$from, .data$id)], 
                                      y = .data$y[match(.data$from, .data$id)],
                                      xend = .data$x[match(.data$to, .data$id)],
                                      yend = .data$y[match(.data$to, .data$id)],
                                      alpha = abs(.data$weight)),
                     color = "grey50", size = edge_width) +
        geom_point(data = nodes_df, aes(x = x, y = y, size = degree), 
                   color = "steelblue", alpha = 0.7) +
        scale_size_continuous(range = c(node_size * 0.5, node_size * 1.5)) +
        labs(title = "Network Plot",
             x = "", y = "") +
        theme_void() +
        theme(legend.position = "none")
    
    if (label_nodes) {
        network_plot <- network_plot + 
            geom_text(data = nodes_df, aes(x = x, y = y, label = id), 
                     size = 2, vjust = -0.5)
    }
    
    cat("Network plot created\n")
    cat("Number of nodes:", nrow(nodes_df), "\n")
    cat("Number of edges:", nrow(edges_df), "\n")
    cat("Correlation threshold:", threshold, "\n")
    
    return(network_plot)
}
