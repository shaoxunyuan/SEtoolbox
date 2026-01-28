#' @title SE_sampleQC: Sample Quality Assessment
#' @description This function assesses the quality of samples in a SummarizedExperiment object. It calculates various quality metrics and identifies potential outlier samples.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for QC analysis. The default value is \code{"TPM"}.
#' @param outlier_threshold Numeric value for outlier detection in standard deviations. Default is 2.
#' @param min_expression Numeric value for minimum expression threshold. Default is 0.
#' @param max_missing Numeric value for maximum proportion of missing values allowed. Default is 0.5.
#' @return A \code{SummarizedExperiment} object with sample QC metrics added to colData and an 'outlier' column indicating outlier samples.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Assess sample quality
#' SE_qc <- SE_sampleQC(SE, assayname = "TPM")
#' 
#' # View QC metrics
#' print(colData(SE_qc))
#' 
#' # Identify outlier samples
#' outliers <- which(colData(SE_qc)$outlier)
#' print("Outlier samples:")
#' print(colnames(SE_qc)[outliers])
#' @export
SE_sampleQC <- function(SE, assayname = "TPM", outlier_threshold = 2, 
                       min_expression = 0, max_missing = 0.5) {
    
    exp_data <- assay(SE, assayname)
    
    sample_sums <- colSums(exp_data, na.rm = TRUE)
    sample_means <- colMeans(exp_data, na.rm = TRUE)
    sample_sds <- apply(exp_data, 2, sd, na.rm = TRUE)
    sample_medians <- apply(exp_data, 2, median, na.rm = TRUE)
    sample_missing <- colMeans(is.na(exp_data))
    sample_zeros <- colMeans(exp_data == 0, na.rm = TRUE)
    
    colData(SE)$total_expression <- sample_sums
    colData(SE)$mean_expression <- sample_means
    colData(SE)$sd_expression <- sample_sds
    colData(SE)$median_expression <- sample_medians
    colData(SE)$missing_ratio <- sample_missing
    colData(SE)$zero_ratio <- sample_zeros
    
    z_sums <- scale(sample_sums)
    z_means <- scale(sample_means)
    z_sds <- scale(sample_sds)
    
    outlier_sums <- abs(z_sums) > outlier_threshold
    outlier_means <- abs(z_means) > outlier_threshold
    outlier_sds <- abs(z_sds) > outlier_threshold
    
    colData(SE)$outlier <- outlier_sums | outlier_means | outlier_sds
    colData(SE)$low_expression <- sample_means < min_expression
    colData(SE)$high_missing <- sample_missing > max_missing
    
    n_outliers <- sum(colData(SE)$outlier)
    n_low_exp <- sum(colData(SE)$low_expression)
    n_high_missing <- sum(colData(SE)$high_missing)
    
    cat("Sample quality assessment completed\n")
    cat("Total samples:", ncol(exp_data), "\n")
    cat("Outlier samples:", n_outliers, "\n")
    cat("Low expression samples:", n_low_exp, "\n")
    cat("High missing samples:", n_high_missing, "\n")
    
    return(SE)
}
