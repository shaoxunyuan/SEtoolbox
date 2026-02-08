#' @title SE_heatmap: Heatmap for SummarizedExperiment with Annotations
#'
#' @description
#' 基于 \code{SummarizedExperiment} 绘制热图：对指定基因（特征）在样本上的表达矩阵做热图展示，
#' 支持样本注释（数值型或分类型）、行列聚类、多种归一化方式，以及可选的格子内数值标签。
#' 热图配色为改良发散色：表达 0 为浅灰，低表达为蓝，高表达为橙红。依赖 \code{ComplexHeatmap}、
#' \code{circlize} 和 \code{RColorBrewer}。
#'
#' @param SE \code{SummarizedExperiment} 对象。必须包含表达矩阵（通过 \code{assay(SE, assayname)} 获取）
#'   以及样本元数据 \code{colData(SE)}；\code{rownames(SE)} 为基因/特征名，\code{colnames(SE)} 为样本名。
#' @param genes_of_interest 字符向量，指定要在热图中展示的基因（行）名称。仅会绘制在 \code{rownames(SE)}
#'   中存在的特征；不存在的会触发 \code{warning} 并被忽略。若全部不存在则 \code{stop}。
#' @param assayname 字符标量，指定使用的 assay 名称，从中提取表达矩阵。默认为 \code{"TPM"}。
#' @param select_classcol 可选。字符向量，指定 \code{colData(SE)} 中用于样本注释的列名（可多列）。
#'   这些列会作为热图顶部的注释条；数值型自动用连续色，分类型使用 \code{RColorBrewer} 的 Set3 调色板。
#'   若为 \code{NULL}，不绘制样本注释。指定的列必须全部存在于 \code{colData(SE)}，否则报错退出。
#' @param normalization 字符，归一化方式。\code{"log"}：log2(x+1)；\code{"scale"}：按行 z-score（常数列置 0）；
#'   \code{"none"}：不变换。使用 \code{match.arg}，默认 \code{"log"}。
#' @param cluster_rows 逻辑值，是否对行（基因）做层次聚类。默认为 \code{TRUE}。
#' @param cluster_cols 逻辑值，是否对列（样本）做层次聚类。默认为 \code{TRUE}。为 \code{FALSE} 时样本顺序
#'   与 \code{assay(SE)} 的列顺序一致。
#' @param show_rownames 逻辑值，是否显示行名（基因名）。默认为 \code{TRUE}。
#' @param show_colnames 逻辑值，是否显示列名（样本名）。默认为 \code{TRUE}。
#' @param use_raster 逻辑值，是否对热图体进行栅格化以减小内存与加快绘制。默认为 \code{TRUE}。
#'   当 \code{show_cell_value = TRUE} 时会被自动设为 \code{FALSE}，以保证格子内数字可见。
#' @param raster_quality 整数，栅格化时的分辨率倍数（如 1 或 2）。值越大图越清晰、导出文件越大。
#'   仅在 \code{use_raster = TRUE} 且未显示格子数字时生效。默认为 \code{2}。
#' @param show_cell_value 逻辑值，是否在每个格子内显示数值（保留两位小数）。默认为 \code{FALSE}。
#'   为 \code{TRUE} 时会关闭栅格化并在格子内绘制数值。
#'
#' @return
#' 返回 \code{draw(ht, ...)} 的结果；热图会直接绘制到当前设备，返回值可用于后续布局微调。
#'
#' @details
#' 特征（基因）若部分不在 \code{SE} 中会提示并取交集；样本注释列必须全部存在于 \code{colData(SE)}，
#' 否则报错。表达矩阵在归一化后会剔除全行含 \code{NA} 的行再绘图。热图颜色：0 对应浅灰 \code{#F7F7F7}，
#' 低值蓝 \code{#3B6FB6}，高值橙红 \code{#D84A3A}。
#'
#' @examples
#' set.seed(1)
#' expr <- matrix(rnorm(200, mean = 5, sd = 2), nrow = 20, ncol = 10)
#' rownames(expr) <- paste0("Gene", 1:20)
#' colnames(expr) <- paste0("Sample", 1:10)
#' sample_info <- DataFrame(group = rep(c("A", "B"), each = 5))
#' SE <- SummarizedExperiment(assays = list(TPM = expr), colData = sample_info)
#'
#' SE_heatmap(SE,
#'   genes_of_interest = c("Gene1", "Gene5", "Gene10"),
#'   select_classcol = "group",
#'   normalization = "log"
#' )
#'
#' SE_heatmap(SE, genes_of_interest = c("Gene1", "Gene2"),
#'   show_cell_value = TRUE, cluster_cols = FALSE)
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar unit
#' @importFrom RColorBrewer brewer.pal
#' @export
SE_heatmap <- function(
  SE,
  genes_of_interest,
  assayname = "TPM",
  select_classcol = NULL,
  normalization = c("log", "scale", "zscore", "none"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  use_raster = TRUE,
  raster_quality = 2,
  show_cell_value = FALSE
) {
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
    library(RColorBrewer)
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

        if (normalization == "zscore") {
          expmat <- t(scale(t(expmat)))
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
          cols <- brewer.pal(3, "Set3")[1]
          names(cols) <- as.character(uq)
          ann_colors[[nm]] <- cols
        } else {
          qs <- quantile(x_no_na, c(0.05, 0.5, 0.95))
          rdbu <- brewer.pal(11, "RdBu")
          if (length(unique(qs)) < 3) {
            qs <- range(x_no_na)
            ann_colors[[nm]] <- colorRamp2(qs, c(rdbu[1], rdbu[11]))
          } else {
            ann_colors[[nm]] <- colorRamp2(qs, c(rdbu[1], rdbu[6], rdbu[11]))
          }
        }
      } else {
        x_no_na <- as.character(x_no_na)
        lev <- unique(x_no_na)
        set3 <- brewer.pal(12, "Set3")
        if (length(lev) == 1) {
          cols <- set3[1]
          names(cols) <- lev
          ann_colors[[nm]] <- cols
        } else {
          cols <- colorRampPalette(set3)(length(lev))
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

  # 热图本体：表达 0 为浅灰，低→蓝，高→橙红（改良发散配色）
  p_lo <- min(quantile(expmat, 0.05, na.rm = TRUE), 0)
  p_hi <- max(quantile(expmat, 0.95, na.rm = TRUE), 0)
  if (p_lo >= p_hi) {
    p_lo <- min(expmat, na.rm = TRUE)
    p_hi <- max(expmat, na.rm = TRUE)
    if (p_lo >= p_hi) p_hi <- p_lo + 1
  }
  if (p_lo < 0 && p_hi > 0) {
    col_fun <- colorRamp2(
      c(p_lo, 0, p_hi),
      c("#3B6FB6", "#F7F7F7", "#D84A3A")
    )
  } else if (p_hi <= 0) {
    col_fun <- colorRamp2(
      c(p_lo, 0),
      c("#3B6FB6", "#F7F7F7")
    )
  } else {
    col_fun <- colorRamp2(
      c(0, p_hi),
      c("#F7F7F7", "#D84A3A")
    )
  }

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
    cell_fun = if (show_cell_value) {
      function(j, i, x, y, w, h, fill) {
        grid.text(sprintf("%.2f", expmat[i, j]), x, y, gp = gpar(fontsize = 8))
      }
    } else NULL,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      legend_height = unit(4, "cm")
    ),
    border = TRUE,
    use_raster = !show_cell_value && use_raster,
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
