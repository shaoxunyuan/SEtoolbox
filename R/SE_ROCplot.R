#' @title Calculate AUC and Plot ROC for Group Comparisons
#'
#' @description
#' 基于 SummarizedExperiment 计算指定特征（基因）在不同样本分组间的 AUC，并绘制 ROC 曲线。
#' 支持自动两两比较所有分组，或指定特定分组对进行比较。
#'
#' @param SE SummarizedExperiment 对象。
#' @param assayname 使用的 assay 名称，默认 "TPM"。
#' @param features 字符向量，指定要分析的基因/特征名称，如 c("Gene1", "Gene2")。
#' @param groupcolname 分组列名，指定 colData(SE) 中用于分组的列，默认 "group"。
#' @param setcompare 可选。分组比较列表。如果不指定（NULL），则自动计算所有分组的两两组合；
#'   如果指定，如 list(c("HC", "LTB"), c("HC", "ATB"))，则只比较指定的分组对。
#'
#' @return 返回列表，包含 AUC 结果表格和 ROC 曲线对象列表。
#'
#' @importFrom pROC roc plot.roc ci.auc
#' @importFrom grDevices rainbow
#' @importFrom stats glm binomial predict as.formula
#' @export
SE_ROCplot <- function(
  SE,
  assayname = "TPM",
  features,
  groupcolname = "group",
  setcompare = NULL
) {
  # 参数检查
  if (!is.character(features) || length(features) == 0) {
    stop("features must be a non-empty character vector.")
  }
  if (!groupcolname %in% colnames(colData(SE))) {
    stop(paste("Column", groupcolname, "not found in colData(SE)."))
  }

  # 数据处理
  expdata <- as.data.frame(t(assay(SE, assayname)))
  expdata <- log2(expdata + 1)
  sample_info <- as.data.frame(colData(SE))

  # 添加分组信息
  expdata$group <- sample_info[, groupcolname]

  # 筛选存在的特征
  available_features <- intersect(features, names(expdata))
  if (length(available_features) == 0) {
    stop("None of the specified features are found in the assay.")
  }
  if (length(available_features) < length(features)) {
    warning("The following features are not found and ignored: ",
            paste(setdiff(features, available_features), collapse = ", "))
  }

  # 构建逻辑回归公式: group ~ feature1 + feature2 + ...
  formula_str <- paste("group ~", paste(available_features, collapse = " + "))
  formula_obj <- as.formula(formula_str)

  # 获取所有分组
  all_groups <- unique(expdata$group)

  # 确定要比较的分组对
  if (is.null(setcompare)) {
    # 自动生成所有两两组合
    if (length(all_groups) < 2) {
      stop("Need at least 2 groups for comparison.")
    }
    compare_pairs <- list()
    for (i in 1:(length(all_groups) - 1)) {
      for (j in (i + 1):length(all_groups)) {
        compare_pairs[[length(compare_pairs) + 1]] <- c(all_groups[i], all_groups[j])
      }
    }
  } else {
    # 使用用户指定的分组对
    if (!is.list(setcompare)) {
      stop("setcompare must be a list of character vectors, e.g., list(c('A', 'B'), c('A', 'C'))")
    }
    compare_pairs <- setcompare
    # 验证分组是否存在
    for (pair in compare_pairs) {
      if (!all(pair %in% all_groups)) {
        stop(paste("Groups", paste(pair, collapse = ", "), "not all found in", groupcolname))
      }
    }
  }

  # 为每个分组对计算 AUC
  results_list <- list()
  roc_objects <- list()

  for (pair in compare_pairs) {
    # 筛选两个分组的数据
    subset_data <- expdata[expdata$group %in% pair, ]

    # 将分组转换为二元变量 (0/1)
    subset_data$group <- ifelse(subset_data$group == pair[2], 1, 0)

    # 拟合逻辑回归模型
    model <- glm(formula_obj, data = subset_data, family = binomial())

    # 获取预测概率
    subset_data$score <- predict(model, type = "response")

    # 计算 ROC
    roc_curve <- roc(
      response = subset_data$group,
      predictor = subset_data$score,
      quiet = TRUE
    )

    comparison_name <- paste(pair, collapse = " vs ")

    # 获取 AUC 的 95% 置信区间
    auc_ci <- ci.auc(roc_curve)

    results_list[[comparison_name]] <- data.frame(
      comparison = comparison_name,
      group1 = pair[1],
      group2 = pair[2],
      auc = as.numeric(roc_curve$auc),
      auc_ci_lower = auc_ci[1],
      auc_ci_upper = auc_ci[3],
      stringsAsFactors = FALSE
    )
    roc_objects[[comparison_name]] <- roc_curve
  }

  # 合并结果
  auc_results <- do.call(rbind, results_list)
  rownames(auc_results) <- NULL

  # 绘制 ROC 曲线（所有比较放在一张图上）
  colors <- rainbow(length(roc_objects))

  for (i in seq_along(roc_objects)) {
    if (i == 1) {
      plot.roc(
        roc_objects[[i]],
        main = "ROC Curves",
        col = colors[i],
        lwd = 2,
        grid = TRUE
      )
    } else {
      plot.roc(roc_objects[[i]], add = TRUE, col = colors[i], lwd = 2)
    }
  }

  # 添加图例
  legend(
    "bottomright",
    legend = names(roc_objects),
    col = colors,
    lwd = 2,
    cex = 0.7,
    title = "Comparison"
  )

  return(list(
    AUC_results = auc_results,
    ROC_objects = roc_objects
  ))
}
