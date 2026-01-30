#' @title Plot Density and Boxplots of Gene Expression  
#'  
#' @description  
#' This function generates density plots and boxplots of gene expression values   
#' from a SummarizedExperiment object, allowing for optional grouping of samples.   
#'  
#' @param SE A SummarizedExperiment object that contains assay data.  
#' @param assayname A character string specifying the assay to be used.   
#' Default is "TPM".  
#' @param group_colname A character string specifying the column name in   
#' the sample metadata that will be used for grouping.   
#' If NULL (default), a combined density plot for all samples is created.  
#'  
#' @import RColorBrewer  
#' @import ggplot2  
#' @import dplyr  
#' @import reshape2  
#' @import tibble  
#' @import scale  
#'  
#' @details The expression values are transformed to a log2 scale before   
#' plotting. Density plots can be generated for individual groups if a   
#' grouping column is provided; otherwise, a global density plot is shown.   
#' Boxplots display the distribution of log2-transformed expression values for each sample,   
#' highlighting any outliers.  
#'  
#' @return A list containing two ggplot2 objects:   
#' \item{plot_density}{a density plot of gene expression.}  
#' \item{plot_boxplot}{a boxplot of log2-transformed expression values for each sample.}  
#'  
#' @examples  
#' # Load necessary libraries  
#' library(SummarizedExperiment)  
#'  
#' # Create a sample SummarizedExperiment object  
#' # (Sample code for creation)  
#' # SE <- SummarizedExperiment(assays = list(TPM = matrix(rnorm(100), nrow=10)),   
#' #                             colData = data.frame(BioSample = 1:10, Group = c("A", "B")))   
#'  
#' # Generate density and boxplot  
#' # plots <- SE_density(SE, assayname = "TPM", group_colname = "Group")  
#' @export  
SE_density <- function(SE, assayname = "TPM", group_colname = NULL) {  
    expdata = as.data.frame(assay(SE, assayname))  

    low1_count <- colSums(expdata < 1, na.rm = TRUE)
    zero_count <- colSums(expdata == 0, na.rm = TRUE)
    n_gene <- nrow(expdata)
    df_low <- data.frame(Sample = colnames(expdata),LowExprFraction = low1_count / n_gene,ZeroExprFraction = zero_count / n_gene)
    df_low$Sample <- factor(df_low$Sample, levels = df_low$Sample[order(df_low$ZeroExprFraction)])
    plot_fraction <- ggplot(df_low, aes(x = Sample, group = 1)) +
          geom_line(aes(y = LowExprFraction, color = "express < 1"), linewidth = 0.6) +
          geom_point(aes(y = LowExprFraction, color = "express < 1"), size = 1.5) +
          geom_line(aes(y = ZeroExprFraction, color = "express = 0"), linewidth = 0.6, linetype = "dashed") +
          geom_point(aes(y = ZeroExprFraction, color = "express = 0"), size = 1.5) +
          scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
          scale_color_manual(values = c("express < 1" = "#4DBBD5", "express = 0" = "#E64B35")) +
          labs(x = "Samples (ordered by fraction of express = 0)", y = "Fraction of genes", title = "") +
          theme_bw() +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top")
    
    sample_info <- as.data.frame(colData(SE))   
    expdata = tibble::rownames_to_column(expdata, var = "feature")  
    expdata_long <- reshape2::melt(expdata, id.vars = "feature", variable.name = "sample", value.name = "express")  
    message("Delete all zero expression from records!")
    expdata_long = expdata_long[expdata_long$express > 0,]  
    expdata_long <- expdata_long %>% mutate(log2express = log2(express))   
    colors <- brewer.pal(8, "Set2")  

    if (!is.null(group_colname)) {  
        expdata_long$group = mapvalues(expdata_long$sample, sample_info$BioSample, sample_info[[group_colname]], warn_missing = FALSE)  
        # 如果提供了分组，则根据分组绘制密度图  
        plot_density = ggplot(expdata_long, aes(x = log2express, fill = group)) +  
            geom_density(alpha = 0.5, color = NA) +  
            scale_fill_manual(values = colors) +  # 使用颜色调色板  
            theme_minimal() +  
            labs(title = "", x = "Expression", y = "") +  
            theme(legend.position = c(0.75, 0.75)) +
            theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12)) 
    } else {  
        # 如果没有分组，则绘制整体密度图  
        plot_density = ggplot(expdata_long, aes(x = log2express)) +  
            geom_density(aes(fill = colors[1]), alpha = 0.5, color = NA) +  
            theme_minimal() +  
            labs(title = "", x = "Log2Expression", y = "") +  
            theme(legend.position = "none")  +
            theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12)) 
    }  
      
    # 绘制箱线图  
    plot_boxplot = ggplot(expdata_long, aes(x = sample, y = log2express, group = sample)) +  
        geom_boxplot(outlier.color = "gray", outlier.shape = 16, outlier.size = 2, fill = "lightgray") +  
        theme_minimal() +  
        labs(title = "", x = "Samples", y = "Log2Expression") +  
        theme(axis.text.x = element_blank(),axis.text.y = element_text(size=12))  

    return(list(plot_fraction = plot_fraction, plot_boxplot = plot_boxplot, plot_density = plot_density))  

}
