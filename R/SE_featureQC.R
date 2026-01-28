#' @title SE_featureQC: Feature Quality Assessment
#' @description This function assesses the quality of features (genes/proteins) in a SummarizedExperiment object. It calculates various quality metrics and identifies low-quality features.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for QC analysis. The default value is \code{"TPM"}.
#' @param min_detectratio Numeric value for minimum detection ratio. Default is 0.1.
#' @param min_variance Numeric value for minimum variance. Default is 0.
#' @param min_cv Numeric value for minimum coefficient of variation. Default is 0.
#' @param max_missing Numeric value for maximum proportion of missing values allowed. Default is 0.5.
#' @return A \code{SummarizedExperiment} object with feature QC metrics added to rowData and a 'low_quality' column indicating low-quality features.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Assess feature quality
#' SE_qc <- SE_featureQC(SE, assayname = "TPM")
#' 
#' # View QC metrics
#' print(rowData(SE_qc))
#' 
#' # Identify low-quality features
#' low_quality <- which(rowData(SE_qc)$low_quality)
#' print("Low-quality features:")
#' print(rownames(SE_qc)[low_quality])
#' @export
SE_featureQC <- function(SE, assayname = "TPM", min_detectratio = 0.1, 
                          min_variance = 0, min_cv = 0, max_missing = 0.5) {
    
    exp_data <- assay(SE, assayname)
    
    feature_sums <- rowSums(exp_data, na.rm = TRUE)
    feature_means <- rowMeans(exp_data, na.rm = TRUE)
    feature_sds <- apply(exp_data, 1, sd, na.rm = TRUE)
    feature_medians <- apply(exp_data, 1, median, na.rm = TRUE)
    feature_missing <- rowMeans(is.na(exp_data))
    feature_zeros <- rowMeans(exp_data == 0, na.rm = TRUE)
    feature_detectratio <- 1 - feature_zeros
    
    feature_cv <- feature_sds / feature_means
    feature_variance <- apply(exp_data, 1, var, na.rm = TRUE)
    
    rowData(SE)$total_expression <- feature_sums
    rowData(SE)$mean_expression <- feature_means
    rowData(SE)$sd_expression <- feature_sds
    rowData(SE)$median_expression <- feature_medians
    rowData(SE)$missing_ratio <- feature_missing
    rowData(SE)$zero_ratio <- feature_zeros
    rowData(SE)$detectratio <- feature_detectratio
    rowData(SE)$cv <- feature_cv
    rowData(SE)$variance <- feature_variance
    
    rowData(SE)$low_detectratio <- feature_detectratio < min_detectratio
    rowData(SE)$low_variance <- feature_variance < min_variance
    rowData(SE)$low_cv <- feature_cv < min_cv
    rowData(SE)$high_missing <- feature_missing > max_missing
    
    rowData(SE)$low_quality <- rowData(SE)$low_detectratio | 
                                  rowData(SE)$low_variance | 
                                  rowData(SE)$low_cv | 
                                  rowData(SE)$high_missing
    
    n_low_detectratio <- sum(rowData(SE)$low_detectratio)
    n_low_variance <- sum(rowData(SE)$low_variance)
    n_low_cv <- sum(rowData(SE)$low_cv)
    n_high_missing <- sum(rowData(SE)$high_missing)
    n_low_quality <- sum(rowData(SE)$low_quality)
    
    cat("Feature quality assessment completed\n")
    cat("Total features:", nrow(exp_data), "\n")
    cat("Low detection ratio features:", n_low_detectratio, "\n")
    cat("Low variance features:", n_low_variance, "\n")
    cat("Low CV features:", n_low_cv, "\n")
    cat("High missing features:", n_high_missing, "\n")
    cat("Total low-quality features:", n_low_quality, "\n")
    cat("Proportion of low-quality features:", round(n_low_quality / nrow(exp_data) * 100, 2), "%\n")
    
    return(SE)
}
