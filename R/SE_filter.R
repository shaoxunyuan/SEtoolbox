#' @title SE_filter: Filter Features in SummarizedExperiment Object
#' @description This function filters features (genes/proteins) in a SummarizedExperiment object based on various criteria including expression level, detection ratio, variance, and coefficient of variation. It provides flexible filtering options to remove low-quality features before downstream analysis.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for filtering. The default value is \code{"TPM"}.
#' @param min_expr Numeric value. Minimum expression threshold. Features with mean expression below this value will be removed. Default is 0.
#' @param min_detectratio Numeric value between 0 and 1. Minimum detection ratio threshold. Features with detection ratio below this value will be removed. Default is 0.
#' @param min_variance Numeric value. Minimum variance threshold. Features with variance below this value will be removed. Default is 0.
#' @param min_cv Numeric value. Minimum coefficient of variation (CV) threshold. Features with CV below this value will be removed. Default is 0.
#' @param group_colname A string representing the column name in \code{colData} that contains group information. This is optional; if provided, filtering will be performed within each group separately.
#' @return A filtered \code{SummarizedExperiment} object with features that pass all specified filtering criteria.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Filter features with minimum expression of 1 and minimum detection ratio of 0.5
#' SE_filtered <- SE_filter(SE, assayname = "TPM", min_expr = 1, min_detectratio = 0.5)
#' 
#' # Filter features with minimum CV of 0.1
#' SE_filtered <- SE_filter(SE, assayname = "TPM", min_cv = 0.1)
#' 
#' # Filter features within each group
#' SE_filtered <- SE_filter(SE, assayname = "TPM", min_expr = 1, group_colname = "group")
#' @export
SE_filter <- function(SE, assayname = "TPM", min_expr = 0, min_detectratio = 0, 
                      min_variance = 0, min_cv = 0, group_colname = NULL) {
    
    exp_data <- assay(SE, assayname)
    
    keep_features <- rep(TRUE, nrow(exp_data))
    
    if (min_expr > 0) {
        if (!is.null(group_colname)) {
            groups <- colData(SE)[[group_colname]]
            expr_by_group <- sapply(unique(groups), function(g) {
                rowMeans(exp_data[, groups == g, drop = FALSE], na.rm = TRUE)
            })
            keep_expr <- apply(expr_by_group, 1, function(x) any(x >= min_expr, na.rm = TRUE))
        } else {
            keep_expr <- rowMeans(exp_data, na.rm = TRUE) >= min_expr
        }
        keep_features <- keep_features & keep_expr
        cat(paste0("After expression filter (min_expr = ", min_expr, "): ", sum(keep_features), " features remaining\n"))
    }
    
    if (min_detectratio > 0) {
        if ("detectratio" %in% colnames(rowData(SE))) {
            keep_detect <- rowData(SE)$detectratio >= min_detectratio
        } else {
            detect_ratio <- rowMeans(exp_data > 0, na.rm = TRUE)
            keep_detect <- detect_ratio >= min_detectratio
        }
        keep_features <- keep_features & keep_detect
        cat(paste0("After detection ratio filter (min_detectratio = ", min_detectratio, "): ", sum(keep_features), " features remaining\n"))
    }
    
    if (min_variance > 0) {
        keep_var <- apply(exp_data, 1, var, na.rm = TRUE) >= min_variance
        keep_features <- keep_features & keep_var
        cat(paste0("After variance filter (min_variance = ", min_variance, "): ", sum(keep_features), " features remaining\n"))
    }
    
    if (min_cv > 0) {
        row_means <- rowMeans(exp_data, na.rm = TRUE)
        row_sds <- apply(exp_data, 1, sd, na.rm = TRUE)
        cv <- row_sds / row_means
        keep_cv <- cv >= min_cv
        keep_features <- keep_features & keep_cv
        cat(paste0("After CV filter (min_cv = ", min_cv, "): ", sum(keep_features), " features remaining\n"))
    }
    
    SE_filtered <- SE[keep_features, ]
    cat(paste0("Final filtered SE: ", nrow(SE_filtered), " features, ", ncol(SE_filtered), " samples\n"))
    
    return(SE_filtered)
}
