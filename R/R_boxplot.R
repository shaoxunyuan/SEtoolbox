#' @title R_boxplot: Boxplot with Scatter and a/b/c Letter Significance
#'
#' @description
#' Draws boxplots with overlaid points (no jitter) and Tukey HSD letter
#' markers (a, b, c, ab, ...) for group comparisons. Input is in "universal"
#' format: each row = sample, each column = feature, last column = \code{class}.
#'
#' @param data A data frame: rows = samples, columns = features, last column =
#'   \code{class} (grouping). All columns except \code{class} are plotted.
#' @param class_col Character. Name of the grouping column. Default \code{"class"}.
#' @param feature_cols Character vector. Which columns to plot as y. Default
#'   \code{NULL} = all columns except \code{class_col}.
#' @param point_alpha Numeric. Point transparency. Default \code{0.5}.
#' @param point_size Numeric. Point size. Default \code{1.2}.
#' @param box_width Numeric. Box width. Default \code{0.5}.
#' @param box_alpha Numeric. Box fill transparency. Default \code{0.6}.
#' @param letter_size Numeric. Significance letter text size. Default \code{4}.
#' @param letter_ymin_frac Numeric. Letter y position as fraction of y range
#'   above each group max. Default \code{0.05}.
#' @param ylab Character. Y-axis label. Default \code{"value"}.
#'
#' @return A \code{ggplot} object (boxplot + points, no jitter, with letters).
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(
#'   gene1 = c(rnorm(20), rnorm(20, 1), rnorm(20, -0.5)),
#'   gene2 = c(rnorm(20, 2), rnorm(20), rnorm(20, 0.5)),
#'   class = rep(c("A", "B", "C"), each = 20)
#' )
#' R_boxplot(dat)
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point geom_text
#'   theme_bw theme element_blank position_identity labs facet_wrap
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by summarise left_join n_distinct
#' @importFrom multcompView multcompLetters
#' @importFrom stats aov TukeyHSD
#' @export
R_boxplot <- function(data,
                      class_col = "class",
                      feature_cols = NULL,
                      point_alpha = 0.5,
                      point_size = 1.2,
                      box_width = 0.5,
                      box_alpha = 0.6,
                      letter_size = 4,
                      letter_ymin_frac = 0.05,
                      ylab = NULL) {
    if (!class_col %in% colnames(data)) {
        stop("Last column or 'class_col' must be the grouping column (e.g. 'class').")
    }

    if (is.null(feature_cols)) {
        feature_cols <- setdiff(colnames(data), class_col)
    } else {
        feature_cols <- intersect(feature_cols, colnames(data))
    }
    if (length(feature_cols) == 0) {
        stop("No feature columns to plot.")
    }

    long <- pivot_longer(
        data,
        cols = all_of(feature_cols),
        names_to = "feature",
        values_to = "value"
    )
    long[[class_col]] <- factor(long[[class_col]])

    letter_list <- list()
    for (f in unique(long$feature)) {
        sub <- long[long$feature == f, , drop = FALSE]
        grps <- sub[[class_col]]
        n_grp <- n_distinct(grps)
        if (n_grp < 2) {
            letter_list[[f]] <- data.frame(
                .class = unique(grps),
                label = NA_character_,
                y_position = NA_real_,
                feature = f,
                stringsAsFactors = FALSE
            )
            names(letter_list[[f]])[1] <- class_col
            next
        }
        fm <- as.formula(paste0("value ~ ", class_col))
        fit <- aov(fm, data = sub)
        tukey <- TukeyHSD(fit)
        pmat <- tukey[[class_col]][, "p adj", drop = TRUE]
        names(pmat) <- rownames(tukey[[class_col]])
        letters_vec <- multcompLetters(pmat)$Letters
        letter_df_grp <- data.frame(
            group = names(letters_vec),
            label = as.character(letters_vec),
            stringsAsFactors = FALSE
        )
        names(letter_df_grp)[1] <- class_col
        y_pos <- sub %>%
            group_by(.data[[class_col]]) %>%
            summarise(y_position = max(value, na.rm = TRUE), .groups = "drop")
        y_range <- diff(range(sub$value, na.rm = TRUE))
        if (is.na(y_range) || y_range == 0) y_range <- 1
        y_pos$y_position <- y_pos$y_position + letter_ymin_frac * y_range
        letter_df_grp <- left_join(letter_df_grp, y_pos, by = class_col)
        letter_df_grp$feature <- f
        letter_list[[f]] <- letter_df_grp
    }
    letter_df <- do.call(rbind, letter_list)
    letter_df <- letter_df[!is.na(letter_df$label), , drop = FALSE]

    p <- ggplot(long, aes(x = .data[[class_col]], y = .data$value, fill = .data[[class_col]])) +
        geom_boxplot(
            width = box_width,
            alpha = box_alpha,
            outlier.shape = NA
        ) +
        geom_point(
            position = position_identity(),
            alpha = point_alpha,
            size = point_size
        )

    if (nrow(letter_df) > 0) {
        p <- p + geom_text(
            data = letter_df,
            aes(x = .data[[class_col]], y = y_position, label = label),
            inherit.aes = FALSE,
            size = letter_size
        )
    }

    if (is.null(ylab)) ylab <- "value"
    p <- p + labs(x = "", y = ylab)

    if (length(unique(long$feature)) > 1) {
        p <- p + facet_wrap(~ feature, scales = "free_y")
    }

    p <- p + theme_bw() +
        theme(
            panel.grid = element_blank(),
            legend.position = "none"
        )

    p
}
