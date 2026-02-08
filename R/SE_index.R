#' @title SE_index: Compute Weighted Gene Index and Group Comparisons
#'
#' @description
#' This function computes a sample-level index (e.g. TB-Index) as a weighted
#' sum of z-scored gene expression from a \code{SummarizedExperiment} object.
#' When \code{group_cols} is provided, it generates boxplots by group,
#' performs pairwise ROC (AUC) and Wilcoxon tests between groups,
#' and outputs ROC curves similar to \code{R_boxplot} style.
#'
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param DEfeature A named numeric vector: names are feature (e.g. gene) IDs,
#'   values are weights (e.g. \code{-1} or \code{1}). Only features present
#'   in \code{SE} are used.
#' @param assayname A string indicating which assay to use. Default is \code{"TPM"}.
#' @param group_cols A character vector of column names in \code{colData(SE)}
#'   used for grouping. If \code{NULL}, only the index table is returned
#'   (no plots or group statistics).
#' @param setcompare Optional. A list of group pairs to compare, e.g.
#'   \code{list(c("HC", "LTB"), c("HC", "ATB"))}. If \code{NULL}, all pairwise
#'   comparisons are performed.
#'
#' @return A list with:
#' \item{index_df}{A data frame with columns \code{Sample}, \code{Index}, and
#'   one column per feature (log2 expression).}
#' \item{boxplot}{When \code{group_cols} is not \code{NULL}, the last boxplot
#'   produced (one per group column). Otherwise \code{NULL}.}
#' \item{roc_plot}{ROC curve plot for group comparisons.}
#' \item{auc_df}{AUC results data frame with CI for each group comparison.}
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
#' x$boxplot          # boxplot by group (if group_cols provided)
#' x$roc_plot         # ROC curves
#' x$auc_df           # AUC results table
#' x$data             # detailed statistics by group pair
#'
#' @importFrom pROC roc plot.roc ci.auc
#' @importFrom grDevices rainbow
#' @importFrom dplyr group_by summarise n rename filter left_join bind_rows
#'   mutate ungroup
#' @importFrom tibble remove_rownames rownames_to_column
#' @export
SE_index <- function(SE, DEfeature, assayname = "TPM", group_cols = NULL, setcompare = NULL) {
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
    boxplot <- NULL
    roc_plot <- NULL
    auc_df <- data.frame()

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

            boxplot <- R_boxplot(
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

            # Determine comparison pairs
            if (is.null(setcompare)) {
                # Generate all pairwise combinations
                comb_pairs <- combn(grps, 2, simplify = FALSE)
            } else {
                # Use user-specified pairs
                comb_pairs <- setcompare
                # Validate groups exist
                for (pair in comb_pairs) {
                    if (!all(pair %in% grps)) {
                        warning(paste("Groups", paste(pair, collapse = ", "), "not all found in", g_col))
                    }
                }
                comb_pairs <- Filter(function(p) all(p %in% grps), comb_pairs)
            }

            # Store ROC objects for plotting
            roc_objects <- list()
            auc_results <- list()

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
                        quiet = TRUE
                    ),
                    silent = TRUE
                )
                if (inherits(roc_res, "try-error")) next

                # Get AUC and CI
                auc_val <- as.numeric(roc_res$auc)
                auc_ci <- pROC::ci.auc(roc_res)
                ci_low <- as.numeric(auc_ci[1])
                ci_high <- as.numeric(auc_ci[3])

                # Store ROC object
                comparison_name <- paste(g1, "vs", g2)
                roc_objects[[comparison_name]] <- roc_res

                # Store AUC result
                auc_results[[comparison_name]] <- data.frame(
                    comparison = comparison_name,
                    group1 = g1,
                    group2 = g2,
                    auc = auc_val,
                    auc_ci_lower = ci_low,
                    auc_ci_upper = ci_high,
                    stringsAsFactors = FALSE
                )

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

            # Plot ROC curves (similar to R_boxplot style)
            if (length(roc_objects) > 0) {
                colors <- rainbow(length(roc_objects))

                for (i in seq_along(roc_objects)) {
                    if (i == 1) {
                        pROC::plot.roc(
                            roc_objects[[i]],
                            main = paste("ROC Curves -", g_col),
                            col = colors[i],
                            lwd = 2,
                            grid = TRUE
                        )
                    } else {
                        pROC::plot.roc(roc_objects[[i]], add = TRUE, col = colors[i], lwd = 2)
                    }
                }

                # Add legend
                legend(
                    "bottomright",
                    legend = names(roc_objects),
                    col = colors,
                    lwd = 2,
                    cex = 0.7,
                    title = "Comparison"
                )

                # Store results
                roc_plot <- recordPlot()
                auc_df <- do.call(rbind, auc_results)
                rownames(auc_df) <- NULL
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
        boxplot = boxplot,
        roc_plot = roc_plot,
        auc_df = auc_df,
        data = data_out
    ))
}
