#' @title Backward Feature Elimination for AUC Optimization
#'
#' @description
#' Performs backward feature elimination to identify the optimal subset of features
#' that maintains high AUC for binary classification. Starting from all features,
#' features are iteratively removed based on which removal causes the least AUC loss.
#'
#' @param SE A \code{SummarizedExperiment} object containing expression data.
#' @param assayname A string indicating which assay to use. Default is "TPM".
#' @param group_colname A string specifying the column name in \code{colData} for grouping.
#' @param feature_of_interest A character vector of feature names to consider for elimination.
#' @param min_features Minimum number of features to retain. Default is 1.
#' @param max_auc_loss Maximum acceptable AUC loss threshold. Default is 0.005.
#'
#' @return A list containing:
#' \item{resultdf}{Data frame with elimination results at each step}
#' \item{plot}{ggplot object showing AUC trend}
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom plyr mapvalues
#' @importFrom dplyr bind_rows
#' @importFrom pROC roc auc ci.auc
#' @importFrom stats glm binomial predict
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_text scale_x_reverse scale_y_continuous coord_cartesian labs theme_bw theme element_text element_blank
#' @export
SE_backwardsearch <- function(SE, assayname = "TPM", group_colname = "group", 
                               feature_of_interest, min_features = 1, max_auc_loss = 0.005) {
    
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

    message("[3/6] Initializing backward selection parameters")
    current_genes <- feature_names
    step <- 1
    backward_df <- data.frame()

    message("[4/6] Evaluating full model AUC")
    df_full <- as.data.frame(expr_mat[, current_genes, drop = FALSE])
    df_full$group <- y
    fit_full <- glm(group ~ ., data = df_full, family = binomial)
    prob_full <- predict(fit_full, type = "response")
    roc_full <- roc(y, prob_full, quiet = TRUE)
    best_auc_prev <- as.numeric(auc(roc_full))

    backward_df <- rbind(backward_df, data.frame(
        Step = 0,
        RemovedGene = "None",
        GeneSet = paste(current_genes, collapse = ","),
        AUC = round(best_auc_prev, 4),
        N_feature = length(current_genes),
        stringsAsFactors = FALSE
    ))

    message("[5/6] Running backward elimination")
    repeat {
        best_remove <- NA
        best_auc <- -Inf
        
        for (g in current_genes) {
            genes_try <- setdiff(current_genes, g)
            if (length(genes_try) < min_features) next
            df <- as.data.frame(expr_mat[, genes_try, drop = FALSE])
            df$group <- y
            fit <- glm(group ~ ., data = df, family = binomial)
            prob <- predict(fit, type = "response")
            roc_obj <- roc(y, prob, quiet = TRUE)
            auc_val <- as.numeric(auc(roc_obj))
            if (auc_val > best_auc) {
                best_auc <- auc_val
                best_remove <- g
            }
        }
        
        auc_loss <- best_auc_prev - best_auc
        if (auc_loss > max_auc_loss) {
            message("  Stop: AUC loss exceeds threshold")
            break
        }

        current_genes <- setdiff(current_genes, best_remove)
        best_auc_prev <- best_auc

        backward_df <- rbind(backward_df, data.frame(
            Step = step,
            RemovedGene = best_remove,
            GeneSet = paste(current_genes, collapse = ","),
            AUC = round(best_auc, 4),
            N_feature = length(current_genes),
            stringsAsFactors = FALSE
        ))

        step <- step + 1
        if (length(current_genes) <= min_features) break
    }

    message("[6/6] Generating AUC trend plot")
    backward_df$AUC_loss <- c(NA, diff(backward_df$AUC) * (-1))

    plot <- ggplot(backward_df, aes(x = N_feature, y = AUC)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        geom_text(aes(label = sprintf("%.4f", AUC)), vjust = -0.8, size = 3) +
        scale_x_reverse(breaks = backward_df$N_feature) +
        scale_y_continuous(breaks = seq(0, 1, by = 0.1), 
                           labels = sprintf("%.2f", seq(0, 1, by = 0.1))) +
        coord_cartesian(ylim = c(0, 1)) +
        labs(title = "Backward Selection: AUC vs Number of Features",
             x = "Number of Remaining Features",
             y = "AUC",
             subtitle = paste("Stopping threshold (max AUC loss) =", max_auc_loss)) +
        theme_bw() +
        theme(plot.title = element_text(face = "bold"),
              panel.grid.minor = element_blank())

    message("Backward selection finished")
    return(list(resultdf = backward_df, plot = plot))
}
