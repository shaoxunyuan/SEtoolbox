#' @title Circular RNA Differential Expression Test Using Beta - Binomial Model
#'
#' @description This function conducts differential expression analysis for circular RNAs (circRNAs).
#' It compares the expression ratios of circRNAs among different groups by utilizing a beta - binomial model.
#' The function also addresses overdispersion in count data and provides adjusted p - values for multiple testing.
#'
#' @param SEcirc A SummarizedExperiment object that holds circular RNA count data.
#' @param SElinear A SummarizedExperiment object that holds linear RNA count data.
#' @param assayname Character string specifying the name of the assay containing count data. Defaults to "Count".
#' @param group_colname Character string indicating the name of the column in \code{colData(SECirc)}
#' that provides group - related information. Defaults to "group".
#' @param select_group Character vector of groups to include, eg c("control","case"); NULL includes all groups (default: NULL).
#'
#' @return A data frame with circRNA names as row names. It contains columns for:
#'         - Mean expression ratios in each group.
#'         - Raw p - values obtained from likelihood ratio tests.
#'         - Adjusted p - values using the Benjamini - Hochberg method.
#'         - The data frame is sorted by raw p - values in ascending order.
#'
#' @details
#' The function first extracts common features and samples between circular and linear RNA data.
#' Then, it calculates the ratio of circular RNA expression to the total (circular + linear) RNA expression
#' for each circRNA. After that, it fits null and alternative beta - binomial models to test for significant
#' differences in these ratios between groups. Finally, it adjusts the p - values for multiple testing
#' using the Benjamini - Hochberg method and sorts the results by raw p - values.
#'
#' @importFrom stats p.adjust
#' @importFrom SummarizedExperiment assay colData
#' @import aod
#'
#' @examples
#' \dontrun{
#' # Example usage
#' results <- SE_circTest(SEcirc = SEcirc, SElinear = SElinear, 
#'                        group_colname = "group", select_group = c("control", "case"))
#' }
#'
#' @export
SE_circTest <- function(SEcirc, SElinear, assayname = "Count", group_colname = "group", select_group = NULL) {
  
  # Package dependency checks
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("The 'SummarizedExperiment' package is required but not installed.")
  }
  if (!requireNamespace("aod", quietly = TRUE)) {
    stop("The 'aod' package is required but not installed.")
  }
  
  # Object type checks
  if (!is(SEcirc, "SummarizedExperiment")) {
    stop("SEcirc must be a SummarizedExperiment object")
  }
  if (!is(SElinear, "SummarizedExperiment")) {
    stop("SElinear must be a SummarizedExperiment object")
  }
  
  # Parameter validity checks
  if (!is.character(assayname) || length(assayname) != 1) {
    stop("assayname must be a single character string")
  }
  if (!is.character(group_colname) || length(group_colname) != 1) {
    stop("group_colname must be a single character string")
  }
  
  # Feature and sample intersection
  common_features <- intersect(rownames(SEcirc), rownames(SElinear))
  common_samples <- intersect(colnames(SEcirc), colnames(SElinear))
  message(sprintf("Processing %d features for %d samples!", length(common_features), length(common_samples)))
  
  # Subset data to common features and samples
  SEcirc <- SEcirc[rownames(SEcirc) %in% common_features, colnames(SEcirc) %in% common_samples]
  SElinear <- SElinear[rownames(SElinear) %in% common_features, colnames(SElinear) %in% common_samples]
  
  # Extract count data and group information
  circ_counts <- ceiling(as.data.frame(assay(SEcirc, assayname)))
  linear_counts <- ceiling(as.data.frame(assay(SElinear, assayname)))
  sample_metadata <- as.data.frame(colData(SEcirc))
  sample_metadata <- sample_metadata[rownames(sample_metadata) %in% common_samples, ]
  
  # Handle group selection
  if (!is.null(select_group)) {
    if (!all(select_group %in% unique(sample_metadata[[group_colname]]))) {
      unknown_groups <- select_group[!select_group %in% unique(sample_metadata[[group_colname]])]
      stop(sprintf("Unknown groups in select_group: %s", paste(unknown_groups, collapse = ", ")))
    }
    message(sprintf("Using selected groups: %s", paste(select_group, collapse = ", ")))
    sample_metadata <- sample_metadata[sample_metadata[[group_colname]] %in% select_group, ]
  } else {
    message("Including all groups in analysis")
  }
  
  group <- factor(sample_metadata[[group_colname]])
  
  # Validate data dimensions
  if (nrow(circ_counts) != nrow(linear_counts) || ncol(circ_counts) != ncol(linear_counts)) {
    stop("circ_counts and linear_counts must have matching dimensions")
  }
  if (length(group) != ncol(circ_counts)) {
    stop("Group vector length must match number of samples")
  }
  
  # Initialize results dataframe
  unique_groups <- unique(group)
  result_cols <- c(paste0("RatioMean_", unique_groups), "p.value", "p.adj")
  results <- data.frame(matrix(nrow = length(common_features), ncol = length(result_cols)))
  colnames(results) <- result_cols
  rownames(results) <- common_features
  
  # Process each circRNA feature
  progress_counter <- 0
  total_features <- length(common_features)
  
  for (feature in common_features) {
    progress_counter <- progress_counter + 1
    
    # Calculate total counts and handle zeros
    circ <- as.numeric(circ_counts[feature, ])
    linear <- as.numeric(linear_counts[feature, ])
    total <- circ + linear
    total[total == 0] <- 1  # Avoid division by zero
    
    # Prepare data for beta - binomial model
    model_data <- data.frame(
      circ = circ,
      total = total,
      group = group
    )
    
    # Calculate mean ratios per group
    for (g in unique_groups) {
      group_indices <- model_data$group == g
      results[feature, paste0("RatioMean_", g)] <- mean(circ[group_indices] / total[group_indices], na.rm = TRUE)
    }
    
    # Fit beta - binomial models with error handling
    p_value <- NA
    tryCatch({
      null_model <- aod::betabin(cbind(circ, total - circ) ~ 1, ~1, data = model_data)
      alt_model <- aod::betabin(cbind(circ, total - circ) ~ group, ~1, data = model_data, control = list(maxit = 100))
      lr_test <- aod::anova(null_model, alt_model)
      # 这里获取p值的方式更明确一些
      p_value <- lr_test@anova.table[["Pr(>Chisq)"]][2]
    }, error = function(e) {
      message(sprintf("Model fitting failed for feature %s: %s", feature, e$message))
    })
    
    results[feature, "p.value"] <- p_value
    
    # Progress update
    if (progress_counter %% 20 == 0) {
      message(sprintf("%d/%d features processed", progress_counter, total_features))
    }
  }
  
  # Adjust p - values
  results$p.adj <- p.adjust(results$p.value, method = "BH", n = sum(!is.na(results$p.value)))
  
  # Sort results by raw p - values in ascending order
  results <- results[order(results$p.value, decreasing = FALSE), ]
  
  return(results)
}