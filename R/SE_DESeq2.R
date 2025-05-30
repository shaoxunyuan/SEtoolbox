#' @title Perform differential expression analysis using DESeq2  
#'  
#' @description  
#' This function conducts differential expression analysis on a given SummarizedExperiment object  
#' using the DESeq2 package. It takes count data from the specified assay and performs  
#' two-group comparisons based on the grouping factor provided in the colData.  
#'  
#' @importFrom DESeq2 DESeqDataSetFromMatrix  
#' @importFrom DESeq2 DESeq  
#' @importFrom DESeq2 results  
#' @import SummarizedExperiment  
#' @param SE A SummarizedExperiment object containing count data.  
#' @param assayname The name of the assay to use for the analysis. Default is "Count".  
#' @param group_colname The name of the column in colData(SE) that contains the factor for grouping samples. Default is "group".  
#' @param set_compare A character vector specifying the groups to compare. It should contain valid group names present in colData.  
#' @return A SummarizedExperiment object with DESeq2 results stored in its metadata.  
#' @examples  
#'  
#' # Load example SummarizedExperiment object  
#' SE = loadSE()   
#'  
#' # Perform differential expression analysis  
#' results_se <- SE_DEseq2(SE = SE, assayname = "Count", group_colname = "group", set_compare = c("hc", "patient"))  
#'  
#' @export  
SE_DEseq2 <- function(SE, assayname = "Count", group_colname = "group", set_compare = c("hc", "patient")) {  

  # Validate that SE is of type SummarizedExperiment  
  if (!is(SE, "SummarizedExperiment")) {  
    stop("The input object must be of SummarizedExperiment type.")  
  }  
  
  # Validate that the specified assay exists  
  if (!assayname %in% assayNames(SE)) {  
    stop("The specified assay does not exist: ", assayname)  
  }  
  
  # Validate that the grouping column exists  
  if (!group_colname %in% colnames(colData(SE))) {  
    stop("The grouping column does not exist in colData: ", group_colname)  
  }  
  
  # Check that all groups in set_compare are present in colData  
  if (!all(set_compare %in% unique(colData(SE)[[group_colname]]))) {  
    stop("One or more groups in 'set_compare' are not present in the specified grouping column ", group_colname, ".")  
  }  

  # Extract count matrix and force to integer (DESeq2 requires integers)  
  countData <- assay(SE, assayname)
  if (any(countData < 0)) {  
    min_value <- abs(min(countData[countData < 0])) + 1  # 获取负值的绝对值加1  
    message("Negative values detected, applying transformation... Adding minimum value: ", min_value)  
    countData <- countData + min_value  # 添加最小值以消除负值  
  }  
  if (!all(countData == floor(countData))) {  
    message("Non-integer values detected, applying ceiling transformation...")  
    countData <- ceiling(countData) + 1  
  }  

  # Extract sample information and dynamically construct the design formula  
  sample_info <- colData(SE)  
  design_formula <- reformulate(group_colname)  # Dynamically adapt to the grouping column name  
    
  # --------------------------  
  # DESeq2 differential analysis  
  # --------------------------  
  # Create DESeq2 object  
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = sample_info, design = design_formula)  
  
  # Run the analysis and suppress redundant information  
  suppressMessages({dds <- DESeq(dds, quiet = TRUE)})  
  
  # --------------------------  
  # Results processing  
  # --------------------------  
  # Get grouping information and validate  
  num_groups <- length(set_compare)  
  
  # Store all differential analysis results in a list  
  all_results <- list()  
  # Perform all possible pairwise comparisons  
  for (i in num_groups:2) {  
    for (j in (i - 1):1) {  
      contrast <- c(group_colname, set_compare[j], set_compare[i])  
       
      if (contrast[2] %in% levels(dds[[group_colname]]) && contrast[3] %in% levels(dds[[group_colname]])) {  
        res <- as.data.frame(results(dds, contrast = contrast))
		res$gene = rownames(res);rownames(res) = NULL
        na_rows <- which(is.na(res$pvalue))  
        res[na_rows, c("log2FoldChange", "lfcSE", "stat", "pvalue", "padj")] <- 0  
        res[na_rows, "pvalue"] <- 1  
        res[na_rows, "padj"] <- 1  

        res <- res[order(res$pvalue, decreasing = FALSE),]  
        comparison_name <- paste0(set_compare[j], "_vs_", set_compare[i])  
        all_results[[comparison_name]] <- res  
      } else {  
        warning(paste("Compare group", contrast[3], " vs ", contrast[2], " not exist"))  
      }  
    }  
  }  

  # Add differential analysis results to the metadata of the new object  
  metadata(SE)$DEresults <- all_results  

  # Return the new object (original SE remains unchanged)  
  return(SE)  
}