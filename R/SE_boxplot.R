#' @title SE_boxplot: Boxplot with Scatter and Significance Letters for SummarizedExperiment.
#'
#' @description
#' This function generates boxplots with overlaid points (no jitter) for specified
#' genes within a \code{SummarizedExperiment} object, and adds Tukey HSD a/b/c
#' significance letters. Optionally performs data normalization.
#'
#' @param SE A \code{SummarizedExperiment} object that contains gene expression data.
#' @param feature_of_interest A character vector specifying the gene names to be plotted.
#' Defaults to \code{c("AAGAB", "ABCA13", "ABCC4", "ABHD2")}.
#' @param assayname A string indicating the assay from which to extract the data.
#' The default value is \code{"TPM"}.
#' @param group_colname A string representing the column name in \code{colData} that
#' holds group information. Defaults to \code{"group"}.
#' @param normalization A string specifying the normalization method.
#' Available options are \code{"none"}, \code{"scale"}, or \code{"log"}. Defaults to \code{"none"}.
#' 
#' @return
#' A \code{ggplot} object: boxplot + scatter (no jitter) with significance letters.
#'
#' @examples
#' # Create a dummy SummarizedExperiment object
#' data_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(data_matrix) <- paste0("Gene", 1:100)
#' colnames(data_matrix) <- paste0("Sample", 1:10)
#' sample_info <- DataFrame(group = rep(c("A", "B"), each = 5))
#' SE <- SummarizedExperiment(assays = list(TPM = data_matrix), colData = sample_info)
#' 
#' # Call the SE_boxplot function
#' plot <- SE_boxplot(SE, feature_of_interest = c("Gene1", "Gene2"), 
#'                    group_colname = "group", normalization = "log")
#' print(plot)
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point geom_text theme_minimal
#'   theme element_rect element_blank position_identity
#' @import ggforce
#' @import multcompView
#' @importFrom multcompView multcompLetters
#' @import tidyverse
#' @import rstatix
#' @export
SE_boxplot <- function(SE,   
                       feature_of_interest = c("AAGAB", "ABCA13", "ABCC4", "ABHD2"),   
                       assayname = "TPM",   
                       group_colname = "group",   
                       normalization = "none") {  

    if (!inherits(SE, "SummarizedExperiment")) {  
        stop("Input SE must be a SummarizedExperiment object.")  
    }  

    expdata = as.data.frame(assay(SE, assayname))  
    expdata_subset <- expdata[rownames(expdata) %in% feature_of_interest, , drop = FALSE]    

    missing_genes <- setdiff(feature_of_interest, rownames(expdata_subset))  
    if (length(missing_genes) > 0) {  
        warning(paste("The following genes are not found in the dataset:",   
                      paste(missing_genes, collapse = ", ")))  
        feature_of_interest <- intersect(feature_of_interest, rownames(expdata_subset))  
        expdata_subset <- expdata_subset[feature_of_interest, , drop = FALSE]  

        # 检查是否仍然有有效基因  
        if (length(feature_of_interest) == 0) {  
            stop("None of the input genes are found in the dataset.")  # 弹出消息并停止执行  
        }  
    }  

    if (normalization == "scale") {  
        samplename = colnames(expdata_subset)  
        expdata_subset <- as.data.frame(t(apply(expdata_subset, 1, scale)))  
        colnames(expdata_subset) = samplename  
    } else if (normalization == "log") {  
        expdata_subset <- log2(expdata_subset + 1)   
    }  

    sample_info <- colData(SE)  
    sample_info[] <- lapply(sample_info, function(x) {  
        if (inherits(x, "integer64")) {  
            return(as.integer(x))  # integer64 to numeric  
        } else {  
            return(x)  
        }  
    })  
    sample_info = as.data.frame(sample_info)  

    exp_data_long <- as.data.frame(expdata_subset) %>%   
                     rownames_to_column(var = "feature") %>%   
                     pivot_longer(cols = -feature, names_to = "sample", values_to = "express")  

    if (!is.null(group_colname) && !(group_colname %in% colnames(sample_info))) {  
        stop("Provided group_colname does not exist in colData of the SummarizedExperiment.")  
    }  

    exp_data_long <- exp_data_long %>%  
                     left_join(sample_info %>% rownames_to_column(var = "sample"), by = "sample") %>%  
                     mutate(group = .data[[group_colname]])  

    MakeSigMarker <- function(onedata) {
        if (n_distinct(onedata$group) < 2) {
            warning("Only one group found. Skipping analysis.")
            return(data.frame(group = unique(onedata$group), label = NA, 
                              y_position = NA, feature = unique(onedata$feature)))
        }
        
        # 检查表达值是否有变异
        if (sd(onedata$express, na.rm = TRUE) == 0) {
            warning("No variation in expression values. Skipping analysis.")
            return(data.frame(group = unique(onedata$group), label = NA, 
                              y_position = NA, feature = unique(onedata$feature)))
        }
        
        # 检查每组样本数是否至少为 2
        group_counts <- onedata %>%
                        group_by(group) %>%
                        summarise(count = n(), .groups = "drop")
        
        if (any(group_counts$count < 2)) {
            warning("Some groups have insufficient sample size (< 2). Skipping analysis.")
            return(data.frame(group = unique(onedata$group), label = NA, 
                              y_position = NA, feature = unique(onedata$feature)))
        }
        
        formula_str <- as.formula("express ~ group")
        anova_result <- aov(formula_str, data = onedata)
         
        # Get pvalue（模型项名为 "group"，与 group_colname 无关）
        tukey_result <- TukeyHSD(anova_result)
        tukey_group <- tukey_result[["group"]]
        p_values <- tukey_group[, "p adj"]
        names(p_values) <- rownames(tukey_group)

        letters <- multcompView::multcompLetters(p_values)$Letters
        letter_df <- data.frame(group = names(letters), label = as.character(letters), stringsAsFactors = FALSE)

        y_range <- diff(range(onedata$express, na.rm = TRUE))
        if (is.na(y_range) || y_range == 0) y_range <- 1
        max_vals <- onedata %>%
                    group_by(group) %>%
                    summarise(y_position = max(express, na.rm = TRUE) + 0.05 * y_range, .groups = "drop")   
        
        results <- left_join(max_vals, letter_df, by = "group")  
        results$feature <- unique(onedata$feature)  
        
        return(results)  
    }  
    
    dfout_summary <- data.frame()

    for (feature in unique(exp_data_long$feature)) {
        onedata <- exp_data_long[exp_data_long$feature == feature, , drop = FALSE]
        out <- tryCatch(MakeSigMarker(onedata), error = function(e) {
            warning("计算显著性标记时出错 (feature: ", feature, "): ", conditionMessage(e))
            data.frame(group = unique(onedata$group), label = NA_character_, y_position = NA_real_, feature = feature, stringsAsFactors = FALSE)
        })
        dfout_summary <- rbind(dfout_summary, out)
    }
    dfout_summary <- dfout_summary[!is.na(dfout_summary$label), , drop = FALSE]

    # 无任何显著性标记时提示
    features_all <- unique(exp_data_long$feature)
    if (nrow(dfout_summary) == 0) {
        warning(
            "未生成任何显著性标记 (a/b/c)。可能原因：\n",
            "  - 每个基因仅有一个分组（需至少 2 组才能做 Tukey HSD）；\n",
            "  - group_colname 在 colData 中对应的分组水平数 < 2；\n",
            "  - 某组样本量过少或全为 NA。\n",
            "请检查 group_colname = \"", group_colname, "\" 及各组的样本量。"
        )
    } else {
        features_with_letters <- unique(dfout_summary$feature)
        missing <- setdiff(features_all, features_with_letters)
        if (length(missing) > 0) {
            warning(
                "以下基因/特征未显示显著性标记: ", paste(missing, collapse = ", "), "。\n",
                "可能原因：该基因仅有一个分组，或该基因下分组数 < 2。"
            )
        }
    }

    plot <- ggplot(exp_data_long, aes(x = group, y = express, fill = group)) +
            geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA) +
            geom_point(position = position_identity(), alpha = 0.5, size = 1.2) +
            geom_text(data = dfout_summary, aes(x = group, y = y_position, label = label), size = 5, inherit.aes = FALSE) +
            facet_row(~ feature, scales = "free", space = "free") +
            theme_minimal() +
            theme(
                legend.position = "none",
                strip.text = element_text(size = 12, colour = "black"),
                panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
                axis.text.x = element_text(size = 12, colour = "black"),
                axis.text.y = element_text(size = 12, colour = "black"),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size = 12, colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()
            )

    # 执行差异分析
    if (n_distinct(exp_data_long$group) >= 2) {
        # 计算每个feature在每个组中的样本量
        group_sample_sizes <- exp_data_long %>% 
            group_by(feature, group) %>% 
            summarise(sample_size = n(), .groups = "drop") %>% 
            pivot_wider(names_from = group, values_from = sample_size, 
                       names_prefix = "n_", values_fill = 0)
        
        # 过滤掉表达值没有变异的 feature
        features_with_variation <- exp_data_long %>%
            group_by(feature) %>%
            summarise(has_variation = sd(express, na.rm = TRUE) > 0, .groups = "drop") %>%
            filter(has_variation) %>%
            pull(feature)
        
        # 只对有变异的 feature 进行差异分析
        if (length(features_with_variation) > 0) {
            exp_data_long_filtered <- exp_data_long %>%
                filter(feature %in% features_with_variation)
            
            # 使用 BH 校正的差异分析
            anova_res_BH <- exp_data_long_filtered %>% 
                dplyr::group_by(feature) %>% 
                anova_test(express ~ group) %>% 
                adjust_pvalue(method = "BH") %>% 
                add_significance()
            
            pairwise_res_BH <- exp_data_long_filtered %>% 
                group_by(feature) %>% 
                tukey_hsd(express ~ group) %>% 
                select(feature, group1, group2, p.adj, p.adj.signif)
            
            # 重命名 pairwise_res 中的列，避免与 anova_res 中的列冲突
            pairwise_res_BH <- pairwise_res_BH %>%
                rename(p.adj.tukey = p.adj, p.adj.signif.tukey = p.adj.signif)
            
            # 将 grouped_anova_test 对象转换为普通数据框
            anova_res_BH_df <- as.data.frame(anova_res_BH)
            
            diff_results_BH <- anova_res_BH_df %>% 
                left_join(pairwise_res_BH, by = "feature") %>% 
                left_join(group_sample_sizes, by = "feature") %>%
                # 添加 group1 和 group2 对应的样本量列
                mutate(
                    group1_n = apply(., 1, function(x) {
                        group1_val <- x["group1"]
                        if (!is.na(group1_val)) {
                            return(as.numeric(x[paste0("n_", group1_val)]))
                        } else {
                            return(NA)
                        }
                    }),
                    group2_n = apply(., 1, function(x) {
                        group2_val <- x["group2"]
                        if (!is.na(group2_val)) {
                            return(as.numeric(x[paste0("n_", group2_val)]))
                        } else {
                            return(NA)
                        }
                    }),
                    # 添加校正方法标识符
                    adjustment_method = "BH"
                ) %>%
                # 调整列顺序：feature, 样本量列, 校正方法, ANOVA 结果, Tukey HSD 结果（包含组样本量）
                select(feature, starts_with("n_"), adjustment_method, DFn, DFd, F, p, p.adj, p.adj.signif, 
                       group1, group1_n, group2, group2_n, p.adj.tukey, p.adj.signif.tukey)
            
            # 使用 FDR 校正的差异分析
            anova_res_FDR <- exp_data_long_filtered %>% 
                dplyr::group_by(feature) %>% 
                anova_test(express ~ group) %>% 
                adjust_pvalue(method = "fdr") %>% 
                add_significance()
            
            # 将 grouped_anova_test 对象转换为普通数据框
            anova_res_FDR_df <- as.data.frame(anova_res_FDR)
            
            pairwise_res_FDR <- exp_data_long_filtered %>% 
                group_by(feature) %>% 
                tukey_hsd(express ~ group) %>% 
                select(feature, group1, group2, p.adj, p.adj.signif) %>%
                rename(p.adj.tukey = p.adj, p.adj.signif.tukey = p.adj.signif)
            
            diff_results_FDR <- anova_res_FDR_df %>% 
                left_join(pairwise_res_FDR, by = "feature") %>% 
                left_join(group_sample_sizes, by = "feature") %>%
                # 添加 group1 和 group2 对应的样本量列
                mutate(
                    group1_n = apply(., 1, function(x) {
                        group1_val <- x["group1"]
                        if (!is.na(group1_val)) {
                            return(as.numeric(x[paste0("n_", group1_val)]))
                        } else {
                            return(NA)
                        }
                    }),
                    group2_n = apply(., 1, function(x) {
                        group2_val <- x["group2"]
                        if (!is.na(group2_val)) {
                            return(as.numeric(x[paste0("n_", group2_val)]))
                        } else {
                            return(NA)
                        }
                    }),
                    # 添加校正方法标识符
                    adjustment_method = "FDR"
                ) %>%
                # 调整列顺序：feature, 样本量列, 校正方法, ANOVA 结果, Tukey HSD 结果（包含组样本量）
                select(feature, starts_with("n_"), adjustment_method, DFn, DFd, F, p, p.adj, p.adj.signif, 
                       group1, group1_n, group2, group2_n, p.adj.tukey, p.adj.signif.tukey)
            
            # 合并 BH 和 FDR 结果到一个表格
            diff_results <- bind_rows(diff_results_BH, diff_results_FDR)
        } else {
            # 所有 feature 都没有变异时，返回空的差异分析结果
            diff_results <- data.frame()
            warning("No variation in expression values for all features. Skipping differential expression analysis.")
        }
    } else {
        # 只有一个组时，返回空的差异分析结果
        diff_results <- data.frame()
    }
    
    # 准备表达矩阵
    expdata <- exp_data_long %>% 
        distinct(feature, sample, .keep_all = TRUE) %>% 
        pivot_wider(id_cols = c(feature, group), 
                   names_from = sample, 
                   values_from = express, 
                   values_fill = list(express = 0)) %>% 
        select(feature, group, everything())
    
    # 准备样本矩阵
    sampledata <- sample_info
    
    # 返回结果列表
    return(list(
        plot = plot, 
        diff_results = diff_results, 
        expdata = expdata,
        sampledata = sampledata
    ))
}
