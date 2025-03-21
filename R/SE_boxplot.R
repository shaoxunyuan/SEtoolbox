#' Generates a boxplot for specified features in a SummarizedExperiment object.
#'
#' @param SE A SummarizedExperiment object.
#' @param feature_of_interest Character vector of gene identifiers.
#' @param assayname The assay name in the SummarizedExperiment object.
#' @param group_col Column name in colData(SE) for grouping.
#' @param normalization Normalization method ("scale", "log", or "none").
#' @export
SE_boxplot <- function(SE, feature_of_interest = c("AAGAB", "ABCA13", "ABCC4", "ABHD2"), assayname = "TPM", group_col = NULL, normalization = "none") {  
    library(SummarizedExperiment)  
    library(reshape2)  
    library(ggplot2)  
    # 检查输入有效性  
    if (!inherits(SE, "SummarizedExperiment")) {  
        stop("Input SE must be a SummarizedExperiment object.")  
    }  
    
    # 提取表达数据  
    exp_data_subset <- assay(SE, assayname)[feature_of_interest, , drop = FALSE]  
    
    # 检查基因是否存在于数据中  
    missing_genes <- setdiff(feature_of_interest, rownames(exp_data_subset))  
    if (length(missing_genes) > 0) {  
        warning(paste("The following genes are not found in the dataset:", paste(missing_genes, collapse = ", ")))  
        feature_of_interest <- intersect(feature_of_interest, rownames(exp_data_subset))  
        exp_data_subset <- exp_data_subset[feature_of_interest, , drop = FALSE]  
    }  

    # 标准化表达数据  
    if (normalization == "scale") {  
        exp_data_subset <- t(apply(exp_data_subset, 1, scale))  # 对每个基因进行 z-score 标准化  
    } else if (normalization == "log") {  
        exp_data_subset <- log2(exp_data_subset)  # 对数据进行 log1p 转换  
    }  # 如果 sel_mode 是 "none"，则保留原始数据  

    # 如果提供了分组且不是 NULL，则获取分组信息  
    if (!is.null(group_col) && group_col %in% names(colData(SE))) {  
        group_info <- colData(SE)[[group_col]]  
        exp_data_melted <- melt(as.matrix(exp_data_subset))  
        exp_data_melted$group <- rep(group_info, each = nrow(exp_data_subset))  
        
        # 绘制 boxplot（带分组）  
        ggplot(exp_data_melted, aes(x=Var1, y=value, fill=group)) +   
              geom_boxplot(outlier.shape = NA) +  
              geom_point(position = position_dodge(width = 0.75), size = 1, color = "black") +  
              labs(title="", x="Feature", y=assayname) +   
              theme(
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                panel.background = element_blank(),  
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),   
                panel.border = element_rect(color = "gray", fill = NA, size = 1)
          ) +
                scale_fill_discrete(name = group_col) +  
                facet_wrap(~ Var1,scales  = "free_x",nrow = 1)  
    } else {  
        # 不进行分组，直接绘制  
        exp_data_melted <- melt(as.matrix(exp_data_subset))       
        # 绘制 boxplot（不带分组）  
        ggplot(exp_data_melted, aes(x=Var1, y=value)) +   
            geom_boxplot(outlier.shape = NA) +  # 省略离群点的显示  
            geom_point(position = position_dodge(width = 0.75), size = 1, color = "black") +
            labs(title="", x="Feature", y=assayname) +   
            theme(
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                panel.background = element_blank(),  
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),   
                panel.border = element_rect(color = "gray", fill = NA, size = 1)
          ) +
                scale_fill_discrete(name = group_col) +  
                facet_wrap(~ Var1,scales  = "free_x",nrow = 1)  
    }  
}  
