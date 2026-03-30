#' @title SE_quantile: Quantile Normalization for SummarizedExperiment Object
#' @description This function performs quantile normalization on expression data in a SummarizedExperiment object.
#' @param se A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to normalize. The default value is \code{"TPM"}.
#' @return A \code{SummarizedExperiment} object with the quantile-normalized data stored in a new assay named after the original assay with "_quantile" suffix.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform quantile normalization on TPM assay
#' SE_qn <- SE_quantile(SE, assayname = "TPM")
#' 
#' # Check the new assay
#' assayNames(SE_qn)
#' @export
SE_quantile <- function(se, assayname = "TPM") {
  mat <- assay(se, assayname)
  mat_log <- log2(mat + 1)
  
  # Handle case with single column
  if (ncol(mat_log) == 1) {
    qn_mat <- mat_log
    rank_mean <- sort(mat_log[, 1])
  } else {
    # Ensure sorted_mat is always a matrix
    sorted_list <- lapply(seq_len(ncol(mat_log)), function(i) sort(mat_log[, i]))
    sorted_mat <- do.call(cbind, sorted_list)
    rank_mean <- rowMeans(sorted_mat, na.rm = TRUE)
    rank_mat <- apply(mat_log, 2, function(x) rank(x, ties.method = "first"))
    qn_mat <- mat_log
    for (i in 1:ncol(mat_log)) {
      qn_mat[, i] <- rank_mean[rank_mat[, i]]
    }
  }
  
  new_name <- paste0(assayname, "_quantile")
  assays(se)[[new_name]] <- qn_mat
  metadata(se)$qn_rank_mean <- rank_mean
  se
}