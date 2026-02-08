# 工具函数：统一格式化输出表格
# 规则：
# 1. 所有数值保留两位小数
# 2. p值使用科学计数法，保留两位有效数字

#' Format numeric values to scientific notation for p-values
#' @param x Numeric vector of p-values
#' @return Formatted character vector
format_pvalue <- function(x) {
  if (is.null(x) || all(is.na(x))) return(x)
  sapply(x, function(val) {
    if (is.na(val)) return(NA_character_)
    # 所有p值统一使用科学计数法，保留两位有效数字
    format(val, scientific = TRUE, digits = 2)
  })
}

#' Format result table with unified numeric formatting
#' 
#' @param df A data frame to format
#' @param pvalue_cols Character vector of column names that contain p-values
#' @return Formatted data frame
#' @keywords internal
format_result_table <- function(df, pvalue_cols = NULL) {
  if (!is.data.frame(df) || nrow(df) == 0) {
    return(df)
  }
  
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  # 自动检测p值列
  if (is.null(pvalue_cols)) {
    pvalue_pattern <- grep("^(p|P)[_.]?(value|val|adj|adjust)|^(p|P)$|pvalue|pval|padj|p_adj|pvalue|adj\\.P\\.Val|P.Value|PValue",
                           colnames(df), ignore.case = TRUE, value = TRUE)
    pvalue_cols <- pvalue_pattern
  }
  
  # 格式化数值列
  for (col in colnames(df)) {
    if (is.numeric(df[[col]])) {
      if (col %in% pvalue_cols) {
        # p值使用科学计数法或保留3位小数
        df[[col]] <- format_pvalue(df[[col]])
      } else {
        # 其他数值保留两位小数
        df[[col]] <- round(df[[col]], 2)
      }
    }
  }
  
  return(df)
}

#' Format numeric columns in a data frame
#' 
#' @param df A data frame
#' @param digits Number of digits to round to (default: 2)
#' @return Data frame with formatted numeric columns
#' @keywords internal
format_numeric_cols <- function(df, digits = 2) {
  if (!is.data.frame(df)) return(df)
  
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  num_cols <- sapply(df, is.numeric)
  
  if (any(num_cols)) {
    df[num_cols] <- lapply(df[num_cols], function(x) round(x, digits))
  }
  
  return(df)
}
