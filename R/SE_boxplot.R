#' @title SE_boxplot: Create Boxplots for SummarizedExperiment Object
#' 
#' @description 
#' This function generates violin and box plots for specified genes within a 
#' \code{SummarizedExperiment} object. It enables the visualization of gene 
#' expression levels across different groups. Optionally, it can perform data 
#' normalization and uses ANOVA with Tukey's post-hoc test to add significance 
#' markers on the generated plot.
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
#' A \code{ggplot} object that represents the violin and box plots, including significance markers.
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
#' @export
SE_boxplot <- function(SE,   
                       feature_of_interest = c("AAGAB", "ABCA13", "ABCC4", "ABHD2"),   
                       assayname = "TPM",   
                       group_colname = "group",   
                       normalization = "none") {  
  
    if (!inherits(SE, "SummarizedExperiment")) {  
        stop("Input SE must be a SummarizedExperiment object.")  
    }  

    exp_data_subset <- assay(SE, assayname)[feature_of_interest, , drop = FALSE]  
     
    missing_genes <- setdiff(feature_of_interest, rownames(exp_data_subset))  
    if (length(missing_genes) > 0) {  
        warning(paste("The following genes are not found in the dataset:",   
                      paste(missing_genes, collapse = ", ")))  
        feature_of_interest <- intersect(feature_of_interest, rownames(exp_data_subset))  
        exp_data_subset <- exp_data_subset[feature_of_interest, , drop = FALSE]  
    }  

    if (normalization == "scale") {  
        exp_data_subset <- t(apply(exp_data_subset, 1, scale))   
    } else if (normalization == "log") {  
        exp_data_subset <- log2(exp_data_subset + 1)   
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
	 
    exp_data_long <- as.data.frame(exp_data_subset) %>%   
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
        
        formula_str <- as.formula("express ~ group")  
        anova_result <- aov(formula_str, data = onedata)  
         
        # Get pvalue   
        tukey_result <- TukeyHSD(anova_result)  
        p_values <- tukey_result[[group_colname]][, "p adj"]  
        names(p_values) <- rownames(tukey_result[[group_colname]])  
         
        letters <- multcompLetters(p_values)$Letters  
        letter_df <- data.frame(group = names(letters), label = letters)  

        max_vals <- onedata %>%  
                    group_by(group) %>%  
                    summarise(y_position = max(express, na.rm = TRUE) + 0.1)   
        
        results <- left_join(max_vals, letter_df, by = "group")  
        results$feature <- unique(onedata$feature)  
        
        return(results)  
    }  
    
    dfout_summary <- data.frame()  
    
    for (feature in unique(exp_data_long$feature)) {  
        onedata <- exp_data_long[exp_data_long$feature == feature,]  
        out <- MakeSigMarker(onedata)  
        dfout_summary <- rbind(dfout_summary, out)  
    }  
    dfout_summary <- dfout_summary[!is.na(dfout_summary$label),]  
   
    plot <- ggplot(exp_data_long, aes(x = group, y = express, fill = feature)) +   
            geom_violin(alpha = 0.5, trim = FALSE, position = position_dodge(width = 0.8), fill = NA, color  = "lightgray") +  
            geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.shape = NA, color = "darkgray", alpha = 0.5) +   
            geom_text(data = dfout_summary, aes(x = group, y = y_position + 0.1 , label = label), size = 5, position = position_nudge(y = 0), inherit.aes = FALSE) +  
            facet_row(~ feature, scales = "free", space = "free") +  
            theme_minimal() +  
            theme(  
                legend.position = "none",  
                legend.title = element_text(size = 12, colour = "black"),  
                strip.text = element_text(size = 12, colour = "black"),  
                panel.border = element_rect(color = "gray", fill = NA, size = 1),  
                axis.text.x = element_text(size = 12, colour = "black"),  
                axis.text.y = element_text(size = 12, colour = "black"),  
                axis.title.x = element_blank(),  
                axis.title.y = element_text(size = 12, colour = "black"),  
                panel.grid.major = element_blank(),  
                panel.grid.minor = element_blank()  
            )  

    return(plot)  
}