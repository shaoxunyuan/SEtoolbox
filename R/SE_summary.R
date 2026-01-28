#' @title SE_summary: Generate Summary Statistics for SummarizedExperiment Object
#' @description This function generates comprehensive summary statistics for a SummarizedExperiment object, including expression statistics, missing value analysis, and feature/sample information.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to summarize. The default value is \code{"TPM"}.
#' @param group_colname A string representing the column name in \code{colData} that contains group information. Default is "group".
#' @return A list containing:  
#' \item{overall}{Overall summary statistics.}  
#' \item{sample_stats}{Sample-level statistics.}  
#' \item{feature_stats}{Feature-level statistics.}  
#' \item{group_stats}{Group-level statistics (if group_colname is provided).}
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Generate summary statistics
#' summary_result <- SE_summary(SE, assayname = "TPM", group_colname = "group")
#' 
#' # View overall summary
#' print(summary_result$overall))
#' 
#' # View sample statistics
#' print(summary_result$sample_stats)
#' 
#' # View group statistics
#' print(summary_result$group_stats)
#' @export
SE_summary <- function(SE, assayname = "TPM", group_colname = "group") {
    
    exp_data <- assay(SE, assayname)
    
    overall_stats <- list(
        n_features = nrow(exp_data),
        n_samples = ncol(exp_data),
        total_expression = sum(exp_data, na.rm = TRUE),
        mean_expression = mean(exp_data, na.rm = TRUE),
        median_expression = median(exp_data, na.rm = TRUE),
        sd_expression = sd(exp_data, na.rm = TRUE),
        min_expression = min(exp_data, na.rm = TRUE),
        max_expression = max(exp_data, na.rm = TRUE),
        total_missing = sum(is.na(exp_data)),
        total_zeros = sum(exp_data == 0, na.rm = TRUE)
    )
    
    sample_stats <- data.frame(
        Sample = colnames(exp_data),
        Total = colSums(exp_data, na.rm = TRUE),
        Mean = colMeans(exp_data, na.rm = TRUE),
        Median = apply(exp_data, 2, median, na.rm = TRUE),
        SD = apply(exp_data, 2, sd, na.rm = TRUE),
        Min = apply(exp_data, 2, min, na.rm = TRUE),
        Max = apply(exp_data, 2, max, na.rm = TRUE),
        Missing = colSums(is.na(exp_data)),
        Zeros = colSums(exp_data == 0, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
    
    feature_stats <- data.frame(
        Feature = rownames(exp_data),
        Total = rowSums(exp_data, na.rm = TRUE),
        Mean = rowMeans(exp_data, na.rm = TRUE),
        Median = apply(exp_data, 1, median, na.rm = TRUE),
        SD = apply(exp_data, 1, sd, na.rm = TRUE),
        Min = apply(exp_data, 1, min, na.rm = TRUE),
        Max = apply(exp_data, 1, max, na.rm = TRUE),
        Missing = rowSums(is.na(exp_data)),
        Zeros = rowSums(exp_data == 0, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
    
    group_stats <- NULL
    if (!is.null(group_colname) && group_colname %in% colnames(colData(SE))) {
        groups <- colData(SE)[[group_colname]]
        group_levels <- unique(groups)
        
        group_stats <- data.frame(
            Group = character(),
            N_Samples = numeric(),
            Mean_Expression = numeric(),
            Median_Expression = numeric(),
            SD_Expression = numeric(),
            stringsAsFactors = FALSE
        )
        
        for (g in group_levels) {
            group_samples <- groups == g
            group_data <- exp_data[, group_samples, drop = FALSE]
            
            group_stats <- rbind(group_stats, data.frame(
                Group = g,
                N_Samples = sum(group_samples),
                Mean_Expression = mean(group_data, na.rm = TRUE),
                Median_Expression = median(group_data, na.rm = TRUE),
                SD_Expression = sd(group_data, na.rm = TRUE),
                stringsAsFactors = FALSE
            ))
        }
    }
    
    cat("Summary statistics generated\n")
    cat("Number of features:", overall_stats$n_features, "\n")
    cat("Number of samples:", overall_stats$n_samples, "\n")
    cat("Total expression:", round(overall_stats$total_expression, 2), "\n")
    
    return(list(
        overall = overall_stats,
        sample_stats = sample_stats,
        feature_stats = feature_stats,
        group_stats = group_stats
    ))
}
