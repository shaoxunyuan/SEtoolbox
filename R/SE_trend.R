#' @title SE_trend: Trend Analysis
#' @description This function performs trend analysis on a SummarizedExperiment object. It identifies increasing, decreasing, or stable trends across features or samples.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for trend analysis. The default value is \code{"TPM"}.
#' @param trend_col A string indicating the column name in colData or rowData containing the variable for trend analysis. Default is "time".
#' @param by A character string specifying what to analyze trends for. Options include "features", "samples". Default is "features".
#' @param method A character string specifying the trend detection method. Options include "linear", "mann_kendall", "spearman". Default is "linear".
#' @param pvalue_threshold Numeric value for p-value threshold to consider a trend significant. Default is 0.05.
#' @return A \code{SummarizedExperiment} object with trend results added to rowData (if by = "features") or colData (if by = "samples").
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Analyze trends across features
#' SE_trend <- SE_trend(SE, assayname = "TPM", trend_col = "time", by = "features")
#' 
#' # View trend results
#' print(rowData(SE_trend)$trend_direction)
#' print(rowData(SE_trend)$trend_significant)
#' @export
SE_trend <- function(SE, assayname = "TPM", trend_col = "time", by = "features", 
                   method = "linear", pvalue_threshold = 0.05) {
    
    exp_data <- assay(SE, assayname)
    
    if (by == "features") {
        if (!trend_col %in% colnames(colData(SE))) {
            stop(paste0("Column '", trend_col, "' not found in colData"))
        }
        
        trend_values <- colData(SE)[[trend_col]]
        
        slopes <- numeric(nrow(exp_data))
        pvalues <- numeric(nrow(exp_data))
        
        for (i in 1:nrow(exp_data)) {
            if (method == "linear") {
                model <- lm(exp_data[i, ] ~ trend_values)
                slopes[i] <- coef(model)[2]
                pvalues[i] <- summary(model)$coefficients[2, 4]
                
            } else if (method == "mann_kendall") {
                mk_result <- Kendall::MannKendall(exp_data[i, ])
                slopes[i] <- mk_result$tau
                pvalues[i] <- mk_result$sl
                
            } else if (method == "spearman") {
                cor_result <- cor.test(exp_data[i, ], trend_values, method = "spearman")
                slopes[i] <- cor_result$estimate
                pvalues[i] <- cor_result$p.value
                
            } else {
                stop("Unknown method. Available options: linear, mann_kendall, spearman")
            }
        }
        
        rowData(SE)$trend_slope <- slopes
        rowData(SE)$trend_pvalue <- pvalues
        rowData(SE)$trend_significant <- pvalues < pvalue_threshold
        rowData(SE)$trend_direction <- ifelse(slopes > 0, "increasing", 
                                              ifelse(slopes < 0, "decreasing", "stable"))
        
        cat("Trend analysis completed for features\n")
        cat("Method:", method, "\n")
        cat("Significant trends:", sum(pvalues < pvalue_threshold), "\n")
        cat("Increasing trends:", sum(slopes > 0 & pvalues < pvalue_threshold), "\n")
        cat("Decreasing trends:", sum(slopes < 0 & pvalues < pvalue_threshold), "\n")
        
    } else if (by == "samples") {
        if (!trend_col %in% colnames(rowData(SE))) {
            stop(paste0("Column '", trend_col, "' not found in rowData"))
        }
        
        trend_values <- rowData(SE)[[trend_col]]
        
        slopes <- numeric(ncol(exp_data))
        pvalues <- numeric(ncol(exp_data))
        
        for (i in 1:ncol(exp_data)) {
            if (method == "linear") {
                model <- lm(exp_data[, i] ~ trend_values)
                slopes[i] <- coef(model)[2]
                pvalues[i] <- summary(model)$coefficients[2, 4]
                
            } else if (method == "mann_kendall") {
                mk_result <- Kendall::MannKendall(exp_data[, i])
                slopes[i] <- mk_result$tau
                pvalues[i] <- mk_result$sl
                
            } else if (method == "spearman") {
                cor_result <- cor.test(exp_data[, i], trend_values, method = "spearman")
                slopes[i] <- cor_result$estimate
                pvalues[i] <- cor_result$p.value
                
            } else {
                stop("Unknown method. Available options: linear, mann_kendall, spearman")
            }
        }
        
        colData(SE)$trend_slope <- slopes
        colData(SE)$trend_pvalue <- pvalues
        colData(SE)$trend_significant <- pvalues < pvalue_threshold
        colData(SE)$trend_direction <- ifelse(slopes > 0, "increasing", 
                                              ifelse(slopes < 0, "decreasing", "stable"))
        
        cat("Trend analysis completed for samples\n")
        cat("Method:", method, "\n")
        cat("Significant trends:", sum(pvalues < pvalue_threshold), "\n")
        cat("Increasing trends:", sum(slopes > 0 & pvalues < pvalue_threshold), "\n")
        cat("Decreasing trends:", sum(slopes < 0 & pvalues < pvalue_threshold), "\n")
        
    } else {
        stop("by must be either 'features' or 'samples'")
    }
    
    return(SE)
}
