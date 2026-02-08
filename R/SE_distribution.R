#' @title SE_distribution: Distribution Visualization for SummarizedExperiment
#'
#' @description
#' Visualizes gene-level and sample-level expression distributions from a
#' \code{SummarizedExperiment} object. Generates density plots and boxplots
#' with optional grouping and significance letters (a/b/c) for group comparisons.
#'
#' @param SE A \code{SummarizedExperiment} object containing expression data.
#' @param assayname A string indicating the assay to use. Default \code{"TPM"}.
#' @param group_colname A string specifying the column name in \code{colData}
#'   for grouping. Default \code{NULL} (no grouping).
#' @param transform A string specifying the transformation method.
#'   Options: \code{"raw"}, \code{"log2"}, \code{"zscore"}, \code{"scale"}.
#'   Default \code{"raw"}.
#' @param drop_zero Logical. Whether to drop zero expression values.
#'   Default \code{TRUE}.
#' @param add_significance Logical. Whether to add Tukey HSD significance
#'   letters (a/b/c) to boxplots when grouping is provided. Default \code{TRUE}.
#'
#' @return A list containing four ggplot objects:
#'   \item{gene_density}{Density plot of gene-level expression}
#'   \item{gene_box_scatter}{Boxplot with points of gene-level expression}
#'   \item{sample_density}{Density plot of sample-level total expression}
#'   \item{sample_box_scatter}{Boxplot with points of sample-level total expression}
#'
#' @examples
#' # Create example data
#' library(SummarizedExperiment)
#' data_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(data_matrix) <- paste0("Gene", 1:100)
#' colnames(data_matrix) <- paste0("Sample", 1:10)
#' sample_info <- DataFrame(group = rep(c("A", "B"), each = 5))
#' SE <- SummarizedExperiment(assays = list(TPM = data_matrix), colData = sample_info)
#'
#' # Generate distribution plots
#' plots <- SE_distribution(SE, group_colname = "group", transform = "log2")
#' print(plots$gene_density)
#' print(plots$gene_box_scatter)
#'
#' @importFrom ggplot2 ggplot aes geom_density geom_boxplot geom_point theme_bw
#'   theme labs element_blank element_text scale_color_manual scale_fill_manual
#'   geom_text
#' @importFrom SummarizedExperiment assay colData
#' @importFrom reshape2 melt
#' @importFrom plyr mapvalues
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr group_by summarise left_join n_distinct
#' @importFrom multcompView multcompLetters
#' @importFrom stats aov TukeyHSD
#' @export
SE_distribution <- function(SE,
                             assayname = "TPM",
                             group_colname = NULL,
                             transform = c("raw", "log2", "zscore", "scale"),
                             drop_zero = TRUE,
                             add_significance = TRUE) {

  transform <- match.arg(transform)

  expmat <- as.matrix(SummarizedExperiment::assay(SE, assayname))
  sample_info <- as.data.frame(SummarizedExperiment::colData(SE))

  if (transform == "log2") expmat <- log2(expmat + 1)
  if (transform == "zscore") expmat <- t(scale(t(expmat)))
  if (transform == "scale") expmat <- scale(expmat)

  ## ================= gene-level long =================
  expdf <- as.data.frame(expmat)
  expdf <- tibble::rownames_to_column(expdf, var = "feature")
  exp_long <- reshape2::melt(expdf,
                             id.vars = "feature",
                             variable.name = "sample",
                             value.name = "express")

  if (drop_zero)
    exp_long <- exp_long[!is.na(exp_long$express) & exp_long$express != 0, ]

  has_group <- !(is.null(group_colname) || is.na(group_colname))

  if (has_group) {
    exp_long$group <- plyr::mapvalues(
      exp_long$sample,
      sample_info$BioSample,
      sample_info[[group_colname]],
      warn_missing = FALSE
    )
    exp_long$group <- as.factor(exp_long$group)
  }

  ## ================= sample-level summary =================
  sample_stat <- data.frame(
    sample = colnames(expmat),
    total_expr = colSums(expmat, na.rm = TRUE),
    detected_gene = colSums(expmat > 0)
  )

  if (has_group) {
    sample_stat$group <- plyr::mapvalues(
      sample_stat$sample,
      sample_info$BioSample,
      sample_info[[group_colname]],
      warn_missing = FALSE
    )
    sample_stat$group <- as.factor(sample_stat$group)
  }

  colors <- RColorBrewer::brewer.pal(8, "Set2")
  title_text <- paste0("Assay: ", assayname, " | Transform: ", transform)
  alpha_val <- 0.3

  ## ================= gene-level plots =================
  if (!has_group) {

    gene_density <- ggplot(exp_long, aes(x = express)) +
      geom_density(fill = colors[1], alpha = alpha_val, linewidth = 0.5) +
      theme_bw() +
      labs(title = paste0(title_text, " | Gene-level"),
           x = paste0("Expression (", transform, ")"),
           y = "density")

    gene_box_scatter <- ggplot(exp_long, aes(x = 1, y = express)) +
      geom_boxplot(fill = colors[1], alpha = alpha_val, outlier.shape = NA) +
      geom_point(color = colors[1], alpha = alpha_val, size = 0.6) +
      theme_bw() +
      labs(title = paste0(title_text, " | Gene-level"),
           x = "",
           y = paste0("Expression (", transform, ")")) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())

  } else {

    gene_density <- ggplot(exp_long, aes(x = express, color = group, fill = group)) +
      geom_density(alpha = alpha_val, linewidth = 0.5) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      theme_bw() +
      labs(title = paste0(title_text, " | Gene-level"),
           x = paste0("Expression (", transform, ")"),
           y = "density") +
      theme(legend.title = element_blank())

    # Create gene-level boxplot with significance letters
    gene_box_scatter <- ggplot(exp_long, aes(x = group, y = express, color = group)) +
      geom_boxplot(fill = NA, linewidth = 0.5, outlier.shape = NA) +
      geom_point(alpha = alpha_val, size = 0.6) +
      scale_color_manual(values = colors) +
      theme_bw() +
      labs(title = paste0(title_text, " | Gene-level"),
           x = group_colname,
           y = paste0("Expression (", transform, ")")) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))

    # Add significance letters if requested
    if (add_significance && n_distinct(exp_long$group) >= 2) {
      # Calculate significance letters for gene-level data
      gene_letters <- calculate_significance_letters(exp_long, "express", "group")
      if (nrow(gene_letters) > 0) {
        gene_box_scatter <- gene_box_scatter +
          geom_text(data = gene_letters,
                    aes(x = group, y = y_position, label = label),
                    inherit.aes = FALSE, size = 3, color = "black")
      }
    }
  }

  ## ================= sample-level plots =================
  if (!has_group) {

    sample_density <- ggplot(sample_stat, aes(x = total_expr)) +
      geom_density(fill = colors[2], alpha = alpha_val, linewidth = 0.5) +
      theme_bw() +
      labs(title = paste0(title_text, " | Sample-level"),
           x = "Total expression per sample",
           y = "density")

    sample_box_scatter <- ggplot(sample_stat, aes(x = 1, y = total_expr)) +
      geom_boxplot(fill = colors[2], alpha = alpha_val, outlier.shape = NA) +
      geom_point(color = colors[2], alpha = alpha_val, size = 2) +
      theme_bw() +
      labs(title = paste0(title_text, " | Sample-level"),
           x = "",
           y = "Total expression per sample") +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())

  } else {

    sample_density <- ggplot(sample_stat, aes(x = total_expr, color = group, fill = group)) +
      geom_density(alpha = alpha_val, linewidth = 0.5) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      theme_bw() +
      labs(title = paste0(title_text, " | Sample-level"),
           x = "Total expression per sample",
           y = "density") +
      theme(legend.title = element_blank())

    # Create sample-level boxplot with significance letters
    sample_box_scatter <- ggplot(sample_stat, aes(x = group, y = total_expr, color = group)) +
      geom_boxplot(fill = NA, linewidth = 0.5, outlier.shape = NA) +
      geom_point(alpha = alpha_val, size = 2) +
      scale_color_manual(values = colors) +
      theme_bw() +
      labs(title = paste0(title_text, " | Sample-level"),
           x = group_colname,
           y = "Total expression per sample") +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))

    # Add significance letters if requested
    if (add_significance && n_distinct(sample_stat$group) >= 2) {
      # Calculate significance letters for sample-level data
      sample_letters <- calculate_significance_letters(sample_stat, "total_expr", "group")
      if (nrow(sample_letters) > 0) {
        sample_box_scatter <- sample_box_scatter +
          geom_text(data = sample_letters,
                    aes(x = group, y = y_position, label = label),
                    inherit.aes = FALSE, size = 3, color = "black")
      }
    }
  }

  return(list(
    gene_density = gene_density,
    gene_box_scatter = gene_box_scatter,
    sample_density = sample_density,
    sample_box_scatter = sample_box_scatter
  ))
}

#' Calculate Significance Letters for Boxplots
#'
#' Internal helper function to calculate Tukey HSD significance letters
#' for group comparisons in boxplots.
#'
#' @param data A data frame containing the data
#' @param value_col Character. Name of the value column
#' @param group_col Character. Name of the grouping column
#' @param letter_ymin_frac Numeric. Letter y position as fraction of y range
#'   above each group max. Default \code{0.05}.
#'
#' @return A data frame with columns: group, label, y_position
#'
#' @importFrom dplyr group_by summarise
#' @importFrom multcompView multcompLetters
#' @importFrom stats aov TukeyHSD
calculate_significance_letters <- function(data, value_col, group_col, letter_ymin_frac = 0.05) {
  # Check if there are at least 2 groups
  n_groups <- dplyr::n_distinct(data[[group_col]])
  if (n_groups < 2) {
    return(data.frame(group = character(), label = character(), y_position = numeric()))
  }

  # Check if any group has less than 2 observations
  group_counts <- data %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop")

  if (any(group_counts$count < 2)) {
    return(data.frame(group = character(), label = character(), y_position = numeric()))
  }

  # Check for variation in values
  if (sd(data[[value_col]], na.rm = TRUE) == 0) {
    return(data.frame(group = character(), label = character(), y_position = numeric()))
  }

  # Perform ANOVA and Tukey HSD
  formula_str <- stats::as.formula(paste0(value_col, " ~ ", group_col))
  fit <- stats::aov(formula_str, data = data)
  tukey <- stats::TukeyHSD(fit)

  # Extract p-values
  tukey_key <- group_col
  if (!(tukey_key %in% names(tukey))) {
    tukey_key <- names(tukey)[1]
  }

  pmat <- tukey[[tukey_key]][, "p adj", drop = TRUE]
  names(pmat) <- rownames(tukey[[tukey_key]])

  # Calculate significance letters
  letters_vec <- multcompView::multcompLetters(pmat)$Letters

  # Create letter data frame
  letter_df <- data.frame(
    group = names(letters_vec),
    label = as.character(letters_vec),
    stringsAsFactors = FALSE
  )
  names(letter_df)[1] <- group_col

  # Calculate y positions
  y_pos <- data %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(y_max = max(.data[[value_col]], na.rm = TRUE), .groups = "drop")

  y_range <- diff(range(data[[value_col]], na.rm = TRUE))
  if (is.na(y_range) || y_range == 0) y_range <- 1

  y_pos$y_position <- y_pos$y_max + letter_ymin_frac * y_range
  names(y_pos)[1] <- group_col

  # Merge letter and position data
  letter_df <- dplyr::left_join(letter_df, y_pos, by = group_col)

  return(letter_df)
}
