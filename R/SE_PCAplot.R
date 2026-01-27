#' @title Generate PCA plots
#'
#' @export
SE_PCAplot <- function(SE,assayname = "TPM",groupname = "group",outlier_threshold = 2,scale = TRUE,feature_of_interesting = NULL,show_legend = FALSE) {

  stopifnot(inherits(SE, "SummarizedExperiment"))

  ## ---------- 数据准备 ----------
  expdata <- SummarizedExperiment::assay(SE, assayname)
  sample_info <- as.data.frame(SummarizedExperiment::colData(SE))

  if (!is.null(feature_of_interesting)) {
    expdata <- expdata[rownames(expdata) %in% feature_of_interesting, , drop = FALSE]
  }

  row_var <- apply(expdata, 1, var, na.rm = TRUE)
  expdata <- expdata[row_var > 0, , drop = FALSE]
  num_feature <- nrow(expdata)

  if (num_feature < 2)
    stop("Not enough variable features for PCA.")

  ## ---------- PCA ----------
  pca_result <- prcomp(t(expdata), scale. = scale)
  pca_data <- as.data.frame(pca_result$x)

  ## ---------- group ----------
  if (is.null(groupname) || !(groupname %in% colnames(sample_info))) {
    pca_data$group <- "n/a"
  } else {
    pca_data$group <- plyr::mapvalues(
      rownames(pca_data),
      rownames(sample_info),
      sample_info[[groupname]],
      warn_missing = FALSE
    )
  }
  pca_data$group <- as.factor(pca_data$group)

  ## ---------- PCA variance ----------
  pca_var <- pca_result$sdev^2
  pca_var_percent <- round(100 * pca_var / sum(pca_var), 2)

  ## ---------- outlier filtering ----------
  mean_pc1 <- mean(pca_data$PC1)
  sd_pc1 <- sd(pca_data$PC1)
  mean_pc2 <- mean(pca_data$PC2)
  sd_pc2 <- sd(pca_data$PC2)

  keep_idx <- (pca_data$PC1 > mean_pc1 - outlier_threshold * sd_pc1) &
              (pca_data$PC1 < mean_pc1 + outlier_threshold * sd_pc1) &
              (pca_data$PC2 > mean_pc2 - outlier_threshold * sd_pc2) &
              (pca_data$PC2 < mean_pc2 + outlier_threshold * sd_pc2)

  pca_data_filter <- pca_data[keep_idx, , drop = FALSE]

  ## ---------- 写回 SE ----------
  sample_info$outlier <- ifelse(
    rownames(sample_info) %in% rownames(pca_data_filter),
    "keep", "delete"
  )
  SummarizedExperiment::colData(SE) <- S4Vectors::DataFrame(sample_info)

  ## ---------- silhouette ----------
  SCvalue <- function(df) {
    if (length(unique(df$group)) < 2) return(NA)
    k <- length(unique(df$group))
    km <- kmeans(df[, c("PC1", "PC2")], centers = k, nstart = 20)
    sil <- cluster::silhouette(km$cluster, dist(df[, c("PC1", "PC2")]))
    round(mean(sil[, "sil_width"]), 3)
  }

  sc_raw <- SCvalue(pca_data)
  sc_filter <- SCvalue(pca_data_filter)

  ## ---------- 绘图函数 ----------
  plot_pca <- function(df, title, show_legend) {
    ggplot(df, aes(PC1, PC2, color = group)) +
      geom_point(size = 3, alpha = 0.8) +
      stat_ellipse(type = "norm", level = 0.95, linetype = 2) +
      theme_bw() +
      labs(
        title = title,
        x = paste0("PC1 (", pca_var_percent[1], "%)"),
        y = paste0("PC2 (", pca_var_percent[2], "%)")
      ) +
      theme(
        legend.position = ifelse(show_legend, "right", "none"),
        legend.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 13, hjust = 0.5)
      )
  }

  pca_plot1 <- plot_pca(
    pca_data,
    paste0("PCA (all) | features=", num_feature,
           " | SC=", sc_raw),
    show_legend
  )

  pca_plot2 <- plot_pca(
    pca_data_filter,
    paste0("PCA (filtered) | kept=", nrow(pca_data_filter),
           " | SC=", sc_filter),
    show_legend
  )

  plot <- cowplot::plot_grid(
    pca_plot1,
    pca_plot2,
    nrow = 1,
    align = "hv",
    labels = c("A", "B")
  )

  ## ---------- 删除样本表 ----------
  sample_outlier <- data.frame(sample_info[sample_info$outlier == "delete", ])

  return(list(
    SE = SE,
    plot = plot,
    sample_outlier = sample_outlier,
    silhouette_raw = sc_raw,
    silhouette_filtered = sc_filter
  ))
}
