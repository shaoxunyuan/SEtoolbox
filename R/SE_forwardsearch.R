#' @title Forward Feature Selection for AUC Optimization
#'
#' @description
#' Performs forward feature selection to identify the optimal subset of features
#' that maximizes AUC for binary classification. Starting from an empty set,
#' features are iteratively added based on which addition yields the highest AUC gain.
#'
#' @param SE A \code{SummarizedExperiment} object containing expression data.
#' @param assayname A string indicating which assay to use. Default is "TPM".
#' @param group_colname A string specifying the column name in \code{colData} for grouping.
#' @param feature_of_interest A character vector of feature names to consider for selection.
#' @param max_features Maximum number of features to select. Default is 10.
#' @param min_auc_gain Minimum AUC gain threshold to continue selection. Default is 0.005.
#'
#' @return A list containing:
#' \item{resultdf}{Data frame with selection results at each step}
#' \item{plot}{ggplot object showing AUC trend}
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom plyr mapvalues
#' @importFrom dplyr bind_rows
#' @importFrom pROC roc auc ci.auc
#' @importFrom stats glm binomial predict
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_text scale_x_continuous scale_y_continuous coord_cartesian labs theme_bw theme element_text element_blank
#' @export
SE_forwardsearch <- function(SE, assayname = "TPM", group_colname = "group", 
                              feature_of_interest, max_features = 10, min_auc_gain = 0.005) {
    
    message("[1/6] Loading libraries")
    suppressPackageStartupMessages({library(plyr);library(dplyr);library(pROC);library(ggplot2)})

    message("[2/6] Preparing expression matrix and group labels")
    expdata <- as.data.frame(t(assay(SE, assayname)))
    expdata <- log2(expdata + 1)
    sample_info <- as.data.frame(colData(SE))
    expdata$group <- mapvalues(rownames(expdata), rownames(sample_info), 
                                sample_info[[group_colname]], warn_missing = FALSE)
    expdata <- expdata[, c(feature_of_interest, "group")]
    expdata$group <- factor(expdata$group)
    stopifnot(nlevels(expdata$group) == 2)

    feature_names <- setdiff(colnames(expdata), "group")
    expr_mat <- as.matrix(expdata[, feature_names])
    y <- expdata$group

    message("[3/6] Initializing forward selection parameters")
    selected_genes <- character(0)
    remaining_genes <- feature_names
    best_auc_prev <- -Inf
    step <- 1
    forward_df <- data.frame()

    message("[4/6] Running forward selection")
    repeat {
        best_gene <- NA
        best_auc <- -Inf
        best_ci <- c(NA, NA, NA)
        
        for (g in remaining_genes) {
            genes_try <- c(selected_genes, g)
            df <- as.data.frame(expr_mat[, genes_try, drop = FALSE])
            df$group <- y
            fit <- glm(group ~ ., data = df, family = binomial)
            prob <- predict(fit, type = "response")
            roc_obj <- roc(y, prob, quiet = TRUE)
            auc_val <- as.numeric(auc(roc_obj))
            ci_val <- as.numeric(ci.auc(roc_obj))
            if (auc_val > best_auc) {
                best_auc <- auc_val
                best_gene <- g
                best_ci <- ci_val
            }
        }
        
        auc_gain <- best_auc - best_auc_prev
        if (auc_gain < min_auc_gain) {
            message("  Stop: AUC gain below threshold")
            break
        }
        
        selected_genes <- c(selected_genes, best_gene)
        remaining_genes <- setdiff(remaining_genes, best_gene)
        best_auc_prev <- best_auc
        
        forward_df <- rbind(forward_df, data.frame(
            Step = step,
            GeneSet = paste(selected_genes, collapse = ","),
            AUC = round(best_auc, 4),
            CI_low = round(best_ci[1], 4),
            CI_mid = round(best_ci[2], 4),
            CI_high = round(best_ci[3], 4),
            N_feature = length(selected_genes),
            stringsAsFactors = FALSE
        ))
        
        step <- step + 1
        if (step > max_features || length(remaining_genes) == 0) break
    }

    message("[5/6] Post-processing results")
    if (nrow(forward_df) > 0) {
        forward_df$AUC_gain <- c(forward_df$AUC[1], diff(forward_df$AUC))
    }

    message("[6/6] Generating AUC trend plot")
    plot <- ggplot(forward_df, aes(x = N_feature, y = AUC)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        geom_text(aes(label = sprintf("%.4f", AUC)), vjust = -0.8, size = 3) +
        scale_x_continuous(breaks = forward_df$N_feature) +
        scale_y_continuous(breaks = seq(0, 1, by = 0.1), 
                           labels = sprintf("%.2f", seq(0, 1, by = 0.1))) +
        coord_cartesian(ylim = c(0, 1)) +
        labs(title = "Forward Selection: AUC vs Number of Features",
             x = "Number of Selected Features",
             y = "AUC",
             subtitle = paste("Stopping threshold (min AUC gain) =", min_auc_gain)) +
        theme_bw() +
        theme(plot.title = element_text(face = "bold"),
              panel.grid.minor = element_blank())

    message("Forward selection finished")
    return(list(resultdf = forward_df, plot = plot))
}
