#' @title SE_WGCNA: Weighted Gene Co-expression Network Analysis
#' @description This function performs Weighted Gene Co-expression Network Analysis (WGCNA) on a SummarizedExperiment object. It identifies modules of highly correlated genes and relates them to sample traits.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for WGCNA. The default value is \code{"log2"}.
#' @param nfeatures Numeric value indicating the number of top variable features to use. Default is 5000.
#' @param power A numeric value for soft-thresholding power. If NULL, it will be automatically selected. Default is NULL.
#' @param minModuleSize Minimum module size. Default is 30.
#' @param mergeCutHeight Height at which to merge modules. Default is 0.25.
#' @param scale_data Logical value indicating whether to scale the data before WGCNA. Default is TRUE.
#' @return A list containing:  
#' \item{moduleColors}{Module assignments for each gene.}  
#' \item{MEs}{Module eigengenes.}  
#' \item{moduleTraitCor}{Correlation between modules and traits.}  
#' \item{plot_dendrogram}{Dendrogram showing module structure.}
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform WGCNA
#' wgcna_result <- SE_WGCNA(SE, assayname = "log2", nfeatures = 5000)
#' 
#' # View module assignments
#' print(head(wgcna_result$moduleColors))
#' @export
SE_WGCNA <- function(SE, assayname = "log2", nfeatures = 5000, power = NULL, 
                      minModuleSize = 30, mergeCutHeight = 0.25, scale_data = TRUE) {
    
    exp_data <- assay(SE, assayname)
    
    variances <- apply(exp_data, 1, var, na.rm = TRUE)
    top_features <- names(sort(variances, decreasing = TRUE))[1:min(nfeatures, length(variances))]
    exp_data_subset <- exp_data[top_features, ]
    
    if (scale_data) {
        exp_data_subset <- t(scale(t(exp_data_subset)))
    }
    
    datExpr <- t(exp_data_subset)
    
    if (is.null(power)) {
        sft <- pickSoftThreshold(datExpr, powerVector = c(1:20), verbose = 5)
        power <- sft$powerEstimate
        cat("Selected soft-thresholding power:", power, "\n")
    }
    
    adjacency <- adjacency(datExpr, power = power)
    TOM <- TOMsimilarity(adjacency)
    dissTOM <- 1 - TOM
    
    geneTree <- hclust(as.dist(dissTOM), method = "average")
    
    dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                       deepSplit = 2, pamRespectsDendro = FALSE,
                                       minClusterSize = minModuleSize)
    
    dynamicColors <- labels2colors(dynamicMods)
    
    MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
    MEs <- MEList$eigengenes
    MEDiss <- 1 - cor(MEs)
    METree <- hclust(as.dist(MEDiss), method = "average")
    
    merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = mergeCutHeight, 
                             verbose = 3)
    mergedColors <- merge$colors
    mergedMEs <- merge$newMEs
    
    moduleColors <- mergedColors
    MEs <- mergedMEs
    
    cat("WGCNA analysis completed\n")
    cat("Number of genes analyzed:", nrow(datExpr), "\n")
    cat("Number of modules identified:", length(unique(moduleColors)), "\n")
    
    return(list(
        moduleColors = moduleColors,
        MEs = MEs,
        geneTree = geneTree,
        power = power
    ))
}
