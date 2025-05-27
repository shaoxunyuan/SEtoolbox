#' @title Circular RNA Differential Expression Test Using Beta-Binomial Model
#'
#' @description This function conducts differential expression analysis for circular RNAs (circRNAs)
#' by comparing expression ratios among different groups using a beta-binomial model.
#'
#' @param SECirc A SummarizedExperiment object containing circRNA count data.
#' @param SELinear A SummarizedExperiment object containing linear RNA count data.
#' @param assayname Character string specifying the assay name with count data (default: "Count").
#' @param group_colname Character string for the group column in colData(SECirc) (default: "group").
#' @param select_group Character vector of groups to include; NULL includes all groups (default: NULL).
#'
#' @return Data frame with circRNA names as row names, containing:
#'         - Mean expression ratios per group
#'         - Raw p-values from likelihood ratio tests
#'         - FDR-adjusted p-values (Benjamini-Hochberg method)
#'
#' @importFrom stats anova p.adjust
#' @importFrom SummarizedExperiment assay colData
#' @import aod
#'
#' @examples
#' \dontrun{
#' results <- SE_circTest(SECirc = SEcirc, SELinear = SElinear, 
#'                       group_colname = "treatment_group",select_group=NULL)
#' }
#'
#' @export
SE_circTest <- function(SECirc, SELinear, assayname = "Count", group_colname = "group", select_group = NULL) {
  
  # Package dependency checks
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("The 'SummarizedExperiment' package is required but not installed.")
  }
  if (!requireNamespace("aod", quietly = TRUE)) {
    stop("The 'aod' package is required but not installed.")
  }
  
  # Object type checks
  if (!is(SECirc, "SummarizedExperiment")) {
    stop("SECirc must be a SummarizedExperiment object")
  }
  if (!is(SELinear, "SummarizedExperiment")) {
    stop("SELinear must be a SummarizedExperiment object")
  }
  
  # Parameter validity checks
  if (!is.character(assayname) || length(assayname) != 1) {
    stop("assayname must be a single character string")
  }
  if (!is.character(group_colname) || length(group_colname) != 1) {
    stop("group_colname must be a single character string")
  }
  
  # Feature and sample intersection
  common_features <- intersect(rownames(SECirc), rownames(SELinear))
  common_samples <- intersect(colnames(SECirc), colnames(SELinear))
  message(sprintf("Processing %d features for %d samples!", length(common_features), length(common_samples)))
  
  # Subset data to common features and samples
  SECirc <- SECirc[rownames(SECirc) %in% common_features, colnames(SECirc) %in% common_samples]
  SELinear <- SELinear[rownames(SELinear) %in% common_features, colnames(SELinear) %in% common_samples]
  
  # Extract count data and group information
  circ_counts <- as.data.frame(assay(SECirc, assayname))
  linear_counts <- as.data.frame(assay(SELinear, assayname))
  sample_metadata <- as.data.frame(colData(SECirc))
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
    
    # Prepare data for beta-binomial model
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
    
    # Fit beta-binomial models with error handling
    p_value <- NA
    tryCatch({
      null_model <- aod::betabin(cbind(circ, total - circ) ~ 1, ~1, data = model_data)
      alt_model <- aod::betabin(cbind(circ, total - circ) ~ group, ~1, data = model_data, control = list(maxit = 100))
      lr_test <- anova(null_model, alt_model)
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
  
  # Adjust p-values
  results$p.adj <- p.adjust(results$p.value, method = "BH", n = sum(!is.na(results$p.value)))
  
  return(results)
}