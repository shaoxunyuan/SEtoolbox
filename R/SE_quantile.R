#' @title SE_quantile: Quantile Normalization for SummarizedExperiment Object 
#' @description This function performs quantile normalization on expression data in a SummarizedExperiment object. 
#' @param se A \code{SummarizedExperiment} object containing gene expression data. 
#' @param assayname A string indicating which assay to normalize. The default value is \code{"TPM"}. 
#' @param reference A numeric vector of the same length as nrow(se) used as reference distribution. If NULL (default), compute reference from data. 
#' @return A \code{SummarizedExperiment} object with the quantile-normalized data stored in a new assay named after the original assay with "_quantile" suffix. 
#' @examples 
#' # Load example SummarizedExperiment object 
#' # SE <- loadSE() 
#' # Perform quantile normalization on TPM assay 
#' # SE_qn <- SE_quantile(SE, assayname = "TPM") 
#' # Using external reference 
#' # ref <- readRDS("reference_rank_mean.rds") 
#' # SE_qn_ref <- SE_quantile(SE, "TPM", reference = ref) 
#' @export 
SE_quantile <- function(se, assayname = "TPM", reference = NULL) { 
  mat <- as.matrix(assay(se, assayname)) 
  n_features <- nrow(mat) 
  
  if (!is.null(reference)) { 
    if (!is.numeric(reference)) stop("reference must be a numeric vector") 
    if (length(reference) != n_features) stop("reference length must equal nrow(se)") 
    rank_mean <- as.vector(reference) 
  } 
  
  if (ncol(mat) == 1) { 
    qn_mat <- mat 
    if (is.null(reference)) rank_mean <- as.vector(sort(mat[, 1])) 
    assays(se)[[paste0(assayname, "_quantile")]] <- qn_mat 
    metadata(se)$qn_rank_mean <- rank_mean 
    return(se) 
  } 
  
  if (is.null(reference)) { 
    sorted_list <- lapply(seq_len(ncol(mat)), function(i) sort(mat[, i])) 
    sorted_mat <- do.call(cbind, sorted_list) 
    rank_mean <- as.vector(rowMeans(sorted_mat, na.rm = TRUE)) 
  } 
  
  rank_mat <- apply(mat, 2, function(x) rank(x, ties.method = "first")) 
  qn_mat <- mat 
  for (i in seq_len(ncol(mat))) { 
    qn_mat[, i] <- rank_mean[rank_mat[, i]] 
  } 
  
  assays(se)[[paste0(assayname, "_quantile")]] <- qn_mat 
  metadata(se)$qn_rank_mean <- rank_mean 
  se 
}