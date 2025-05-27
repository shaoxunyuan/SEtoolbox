#' @title Identify Highly Variable Genes  
#'  
#' @description  
#' This function identifies highly variable genes from a SummarizedExperiment object  
#' using different methods such as coefficient of variation squared (cv2), variance stabilizing transformation (vst), or mean-variance trend (mvp).  
#'  
#' @param SE A SummarizedExperiment object containing gene expression data.  
#' @param assayname Character string specifying the assay name in the SummarizedExperiment object. Default is "TPM".  
#' @param method Character string specifying the method to use for identifying highly variable genes.   
#' Options are "cv2", "vst", and "mvp". Default is "cv2".  
#' @param ntop Integer specifying the number of top highly variable genes to return. Default is 2000.  
#' @param span Numeric specifying the span for loess fitting in the mean-variance trend method. Default is 0.3.  
#' @param ... Additional arguments passed to the functions.  
#'  
#' @return A character vector of highly variable gene names.  
#' @export  
#'  
#' @examples  
#' # Example usage:  
#' hvg <- SE_HVG(my_se, assayname = "TPM", method = "cv2", ntop = 2000)  
SE_HVG <- function(SE, assayname = "TPM",   
                   method = c("cv2", "vst", "mvp"),   
                   ntop = 2000,   
                   span = 0.3,  
                   ...) {  
  # Check for required packages  
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {  
    stop("The SummarizedExperiment package is required. Please install it.")  
  }  
  
  # Validate input  
  if (!inherits(SE, "SummarizedExperiment")) {  
    stop("Input must be a SummarizedExperiment object.")  
  }  
  
  if (!(assayname %in% assayNames(SE))) {  
    stop(paste("Assay name '", assayname, "' does not exist in SE.", sep = ""))  
  }  
  
  method <- match.arg(method)  
  
  # Get expression matrix  
  expr <- assay(SE, assayname)  
  
  # Calculate highly variable genes  
  switch(method,  
    "cv2" = {  
      # Calculate squared coefficient of variation (CV²)  
      mean_expr <- rowMeans(expr)  
      var_expr <- apply(expr, 1, var)  
      cv2 <- var_expr / mean_expr^2  
      
      # Filter low-expressed genes  
      expr_filtered <- expr[mean_expr > 0.1, ]  
      mean_filtered <- rowMeans(expr_filtered)  
      var_filtered <- apply(expr_filtered, 1, var)  
      cv2_filtered <- var_filtered / mean_filtered^2  
      
      # Order by CV² and select top ntop genes  
      hvg_idx <- order(cv2_filtered, decreasing = TRUE)[1:min(ntop, length(cv2_filtered))]  
      hvg_genes <- rownames(expr_filtered)[hvg_idx]  
    },  
    
    "vst" = {  
      # Use variance stabilizing transformation method (requires DESeq2 package)  
      if (!requireNamespace("DESeq2", quietly = TRUE)) {  
        stop("The vst method requires the DESeq2 package. Please install it.")  
      }  
      
      # Create DESeqDataSet  
      dds <- DESeq2::DESeqDataSetFromMatrix(  
        countData = round(expr),  
        colData = colData(SE),  
        design = ~ 1  # Not using grouping information  
      )  
      
      # Perform VST transformation  
      vsd <- DESeq2::vst(dds, blind = TRUE)  
      
      # Calculate gene variance and select highly variable genes  
      gene_var <- rowVars(assay(vsd))  
      hvg_idx <- order(gene_var, decreasing = TRUE)[1:min(ntop, length(gene_var))]  
      hvg_genes <- rownames(vsd)[hvg_idx]  
    },  
    
    "mvp" = {  
      # Use mean-variance trend method  
      mean_expr <- rowMeans(expr)  
      var_expr <- apply(expr, 1, var)  
      
      # Fit mean-variance trend  
      df <- data.frame(mean = log10(mean_expr), variance = log10(var_expr))  
      df <- df[complete.cases(df), ]  
      
      loess_fit <- loess(variance ~ mean, data = df, span = span)  
      expected_var <- predict(loess_fit, newdata = df)  
      residual_var <- df$variance - expected_var  
      
      # Select genes with largest residual variance  
      hvg_idx <- order(residual_var, decreasing = TRUE)[1:min(ntop, length(residual_var))]  
      hvg_genes <- rownames(expr)[df$mean[hvg_idx]]  
    }  
  )  
  
  # Return names of highly variable genes  
  return(hvg_genes)  
}