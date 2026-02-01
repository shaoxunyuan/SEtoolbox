#' @title SE_heatmap: Heatmap for SummarizedExperiment with Annotations
#'
#' @description
#' Draws a heatmap for specified genes from a \code{SummarizedExperiment} object,
#' with optional sample annotations (numeric or categorical), row/column clustering,
#' and normalization (log, scale, or none). Uses \code{ComplexHeatmap}.
#'
#' @param SE A \code{SummarizedExperiment} object containing expression data and
#'   sample metadata in \code{colData}.
#' @param genes_of_interest Character vector of gene (row) names to plot. Must
#'   exist in \code{rownames(SE)}; missing genes trigger a warning and are dropped.
#' @param assayname Assay name to use. Default is \code{"TPM"}.
#' @param select_classcol Optional character vector of \code{colData} column names
#'   for sample annotations. \code{NULL} for no annotation.
#' @param normalization One of \code{"log"}, \code{"scale"}, or \code{"none"}.
#'   \code{"log"} = log2(x+1); \code{"scale"} = z-score per row (constant row → 0).
#'   Default is \code{"log"}.
#' @param cluster_rows Whether to cluster rows. Default is \code{TRUE}.
#' @param cluster_cols Whether to cluster columns. Default is \code{TRUE}.
#' @param show_rownames Whether to show row names. Default is \code{TRUE}.
#' @param show_colnames Whether to show column names. Default is \code{TRUE}.
#' @param use_raster Whether to rasterize the heatmap body. Default is \code{TRUE}.
#' @param raster_quality Raster quality (1 or 2). Default is \code{2}.
#'
#' @return
#' The result of \code{draw(ht, ...)} (the heatmap is plotted; return value is
#' typically used for further layout tweaks).
#'
#' @examples
#' # Build a minimal SummarizedExperiment
#' set.seed(1)
#' expr <- matrix(rnorm(200, mean = 5, sd = 2), nrow = 20, ncol = 10)
#' rownames(expr) <- paste0("Gene", 1:20)
#' colnames(expr) <- paste0("Sample", 1:10)
#' sample_info <- DataFrame(group = rep(c("A", "B"), each = 5))
#' SE <- SummarizedExperiment(assays = list(TPM = expr), colData = sample_info)
#'
#' # Heatmap for selected genes with group annotation, log normalization
#' SE_heatmap(SE,
#'   genes_of_interest = c("Gene1", "Gene5", "Gene10"),
#'   select_classcol = "group",
#'   normalization = "log"
#' )
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar unit
#' @export
SE_heatmap <- function(
  SE,
  genes_of_interest,
  assayname = "TPM",
  select_classcol = NULL,
  normalization = c("log", "scale", "none"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  use_raster = TRUE,
  raster_quality = 2
) {
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
  })

  if (!inherits(SE, "SummarizedExperiment")) {
    stop("Input SE must be a SummarizedExperiment object.")
  }

  normalization <- match.arg(normalization)

  # 特征（基因）：部分不在 SE 中则提示并取交集，全部不在则报错退出
  all_genes <- rownames(SE)
  genes_use <- intersect(genes_of_interest, all_genes)
  missing_genes <- setdiff(genes_of_interest, all_genes)
  if (length(missing_genes) > 0) {
    warning("The following features are not found in SE and are dropped: ",
            paste(missing_genes, collapse = ", "))
  }
  if (length(genes_use) == 0) {
    stop("None of the features in genes_of_interest are present in SE. Please check genes_of_interest.")
  }

  # 样本分组列（select_classcol）：用户指定，必须全部存在于 colData，否则报错退出（不取交集）
  if (!is.null(select_classcol)) {
    coldata_nms <- colnames(colData(SE))
    invalid_cols <- setdiff(select_classcol, coldata_nms)
    if (length(invalid_cols) > 0) {
      stop("The following sample annotation columns are not in colData(SE): ",
           paste(invalid_cols, collapse = ", "),
           ". Please check select_classcol. Available columns: ",
           paste(coldata_nms, collapse = ", "))
    }
  }

  # 用整数索引子集 SE，避免用字符子集时不存在的行导致 "index out of bounds"
  row_idx <- match(genes_use, rownames(SE))
  row_idx <- row_idx[!is.na(row_idx)]
  expmat <- as.matrix(assay(SE[row_idx, ], assayname))
  rownames(expmat) <- rownames(SE)[row_idx]

  if (normalization == "log") {
    expmat <- log2(expmat + 1)
  }

  if (normalization == "scale") {
    expmat <- t(apply(expmat, 1, function(x) {
      if (sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
      as.numeric(scale(x))
    }))
    colnames(expmat) <- colnames(assay(SE))
  }

  expmat <- expmat[rowSums(is.na(expmat)) == 0, , drop = FALSE]
  if (nrow(expmat) == 0) {
    stop("After removing rows with NA, no rows remain for the heatmap.")
  }

  ha <- NULL
  if (!is.null(select_classcol)) {
    sampleinfo <- as.data.frame(colData(SE))
    ann_df <- sampleinfo[, select_classcol, drop = FALSE]
    ann_df <- ann_df[colnames(expmat), , drop = FALSE]

    ann_colors <- list()
    for (nm in colnames(ann_df)) {
      x <- ann_df[[nm]]
      x_no_na <- x[!is.na(x)]

      if (is.numeric(x_no_na)) {
        uq <- unique(x_no_na)
        if (length(uq) == 1) {
          cols <- c("#4a90a4")
          names(cols) <- as.character(uq)
          ann_colors[[nm]] <- cols
        } else {
          qs <- quantile(x_no_na, c(0.05, 0.5, 0.95))
          if (length(unique(qs)) < 3) {
            qs <- range(x_no_na)
            ann_colors[[nm]] <- colorRamp2(qs, c("#26456e", "#c94b3a"))
          } else {
            ann_colors[[nm]] <- colorRamp2(
              qs,
              c("#26456e", "#faf8f5", "#c94b3a")
            )
          }
        }
      } else {
        x_no_na <- as.character(x_no_na)
        lev <- unique(x_no_na)
        if (length(lev) == 1) {
          cols <- c("#3d8b82")
          names(cols) <- lev
          ann_colors[[nm]] <- cols
        } else {
          # 柔和、易区分的多色方案（Paul Tol 风格）
          pal <- c("#4477aa", "#66ccee", "#228833", "#ccbb44", "#ee6677",
                   "#aa3377", "#bbbbbb")
          cols <- colorRampPalette(pal)(length(lev))
          names(cols) <- lev
          ann_colors[[nm]] <- cols
        }
      }
    }

    ha <- HeatmapAnnotation(
      df = ann_df,
      col = ann_colors,
      annotation_name_side = "left",
      annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
      gap = unit(1.5, "mm")
    )
  }

  # 热图本体：蓝–米白–红 发散配色，更柔和
  col_fun <- colorRamp2(
    c(quantile(expmat, 0.05), 0, quantile(expmat, 0.95)),
    c("#26456e", "#faf8f5", "#c94b3a")
  )

  heatmap_name <- paste0(assayname, "_", normalization)

  ht <- Heatmap(
    expmat,
    name = heatmap_name,
    col = col_fun,
    top_annotation = ha,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_names = show_rownames,
    show_column_names = show_colnames,
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 45,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      legend_height = unit(4, "cm")
    ),
    border = TRUE,
    use_raster = use_raster,
    raster_quality = raster_quality,
    row_title = "Genes",
    column_title = "Samples",
    row_title_gp = gpar(fontsize = 11, fontface = "bold"),
    column_title_gp = gpar(fontsize = 11, fontface = "bold")
  )

  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legend = TRUE
  )
}
