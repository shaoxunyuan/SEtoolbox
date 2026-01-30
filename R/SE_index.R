#' @title SE_index: Compute Weighted Gene Index and Group Comparisons
#'
#' @description
#' This function computes a sample-level index (e.g. TB-Index) as a weighted
#' sum of z-scored gene expression from a \code{SummarizedExperiment} object.
#' When \code{group_cols} is provided, it generates boxplots by group and
#' performs pairwise ROC (AUC) and Wilcoxon tests between groups.
#'
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param DEfeature A named numeric vector: names are feature (e.g. gene) IDs,
#'   values are weights (e.g. \code{-1} or \code{1}). Only features present
#'   in \code{SE} are used.
#' @param assayname A string indicating which assay to use. Default is \code{"TPM"}.
#' @param group_cols A character vector of column names in \code{colData(SE)}
#'   used for grouping. If \code{NULL}, only the index table is returned
#'   (no plots or group statistics).
#'
#' @return A list with:
#' \item{index_df}{A data frame with columns \code{Sample}, \code{Index}, and
#'   one column per feature (log2 expression).}
#' \item{plot}{When \code{group_cols} is not \code{NULL}, the last boxplot
#'   produced (one per group column). Otherwise \code{NULL}.}
#' \item{data}{When \code{group_cols} is not \code{NULL}, a data frame of
#'   group comparison statistics (N, mean, SD, AUC, CI, P_value, P_adj).
#'   Otherwise an empty data frame with the same columns.}
#'
#' @examples
#' # Build a named numeric vector of gene weights (e.g. from differential analysis)
#' DEfeatures <- c(
#'   SLC39A8 = -1, IQGAP2 = 1, DDX60L = -1, QRICH1 = 1,
#'   PPP1R7 = -1, SOD2 = 1, PLEK = -1, CCT3 = 1
#' )
#'
#' # Load example SE (or use your own filtered SE)
#' SE <- loadSE()
#' # Compute index and group comparisons; use group column name in colData(SE)
#' x <- SE_index(SE = SE, assayname = "TPM", group_cols = "group", DEfeature = DEfeatures)
#' head(x$index_df)   # sample-level index
#' x$plot             # boxplot by group (if group_cols provided)
#' x$data             # AUC and Wilcoxon by group pair
#'
#' @importFrom pROC roc
#' @importFrom dplyr group_by summarise n rename filter left_join bind_rows
#'   mutate ungroup
#' @importFrom tibble remove_rownames rownames_to_column
#' @export
SE_index <- function(SE, DEfeature, assayname = "TPM", group_cols = NULL) {
    suppressPackageStartupMessages({
        library(pROC)
        library(ggplot2)
        library(dplyr)
        library(tibble)
    })

    features <- names(DEfeature)
    exp_raw <- as.matrix(assay(SE, assayname))
    common_f <- intersect(features, rownames(exp_raw))
    if (length(common_f) == 0) {
        stop("SE 对象中未找到 DEfeature 中的任何 ID，请检查 assay 与特征名。")
    }

    exp_log <- log2(exp_raw[common_f, , drop = FALSE] + 1)

    z_scores <- t(apply(exp_log, 1, function(x) {
        if (all(x == x[1])) return(rep(0, length(x)))
        as.numeric(scale(x))
    }))
    colnames(z_scores) <- colnames(exp_log)

    weights <- DEfeature[common_f]
    index_vec <- as.numeric(t(z_scores) %*% weights[rownames(z_scores)])

    output_df <- data.frame(
        Sample = colnames(exp_log),
        Index = index_vec,
        t(exp_log),
        stringsAsFactors = FALSE
    ) %>% tibble::remove_rownames()

    res_boxplots <- list()
    audit_table <- list()
    plot <- NULL

    if (!is.null(group_cols)) {
        meta <- as.data.frame(colData(SE)) %>% tibble::rownames_to_column("Sample")
        audit_df <- dplyr::left_join(
            output_df[, c("Sample", "Index")],
            meta,
            by = "Sample"
        )

        for (g_col in group_cols) {
            if (!(g_col %in% colnames(audit_df))) next

            sub_df <- audit_df %>% dplyr::filter(!is.na(.data[[g_col]]))
            if (nrow(sub_df) < 2) next
            sub_df[[g_col]] <- factor(sub_df[[g_col]])
            grps <- levels(sub_df[[g_col]])
            if (length(grps) < 2) next

            plot <- R_boxplot(
                data = sub_df,
                class_col = g_col,
                feature_cols = "Index",
                ylab = "TB-Index Score"
            )

            stat_df <- sub_df %>%
                dplyr::group_by(.data[[g_col]]) %>%
                dplyr::summarise(
                    n = dplyr::n(),
                    mean = mean(.data$Index, na.rm = TRUE),
                    sd = sd(.data$Index, na.rm = TRUE),
                    .groups = "drop"
                ) %>%
                dplyr::rename(group = 1)

            comb_pairs <- combn(grps, 2, simplify = FALSE)

            for (pair in comb_pairs) {
                g1 <- pair[1]
                g2 <- pair[2]

                pair_df <- sub_df %>% dplyr::filter(.data[[g_col]] %in% c(g1, g2))
                pair_df$Response <- factor(pair_df[[g_col]], levels = c(g1, g2))

                roc_res <- try(
                    pROC::roc(
                        response = pair_df$Response,
                        predictor = pair_df$Index,
                        levels = c(g1, g2),
                        direction = "auto",
                        quiet = TRUE,
                        ci = TRUE
                    ),
                    silent = TRUE
                )
                if (inherits(roc_res, "try-error")) next

                auc_val <- as.numeric(roc_res$auc)
                ci_low <- as.numeric(roc_res$ci[1])
                ci_high <- as.numeric(roc_res$ci[3])

                if (!is.na(auc_val) && auc_val < 0.5) {
                    auc_val <- 1 - auc_val
                    ci_low_new <- 1 - ci_high
                    ci_high_new <- 1 - ci_low
                    ci_low <- ci_low_new
                    ci_high <- ci_high_new
                }

                p_val <- try(stats::wilcox.test(Index ~ Response, data = pair_df)$p.value, silent = TRUE)
                if (inherits(p_val, "try-error")) p_val <- NA_real_

                s1 <- stat_df %>% dplyr::filter(.data$group == g1)
                s2 <- stat_df %>% dplyr::filter(.data$group == g2)

                audit_table[[length(audit_table) + 1]] <- data.frame(
                    GroupCol = g_col,
                    Group1 = g1,
                    Group2 = g2,
                    N1 = s1$n,
                    N2 = s2$n,
                    Mean1 = s1$mean,
                    Mean2 = s2$mean,
                    SD1 = s1$sd,
                    SD2 = s2$sd,
                    AUC = auc_val,
                    CI_low = ci_low,
                    CI_high = ci_high,
                    P_value = as.numeric(p_val),
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    data_out <- dplyr::bind_rows(audit_table)

    if (nrow(data_out) > 0) {
        data_out <- data_out %>%
            dplyr::group_by(GroupCol) %>%
            dplyr::mutate(P_adj = p.adjust(P_value, method = "BH")) %>%
            dplyr::ungroup()

        num_cols <- sapply(data_out, is.numeric)
        data_out[num_cols] <- lapply(data_out[num_cols], function(x) round(x, 2))
    } else {
        data_out <- data.frame(
            GroupCol = character(),
            Group1 = character(),
            Group2 = character(),
            N1 = integer(),
            N2 = integer(),
            Mean1 = numeric(),
            Mean2 = numeric(),
            SD1 = numeric(),
            SD2 = numeric(),
            AUC = numeric(),
            CI_low = numeric(),
            CI_high = numeric(),
            P_value = numeric(),
            P_adj = numeric(),
            stringsAsFactors = FALSE
        )
    }

    return(list(
        index_df = output_df,
        plot = plot,
        data = data_out
    ))
}
