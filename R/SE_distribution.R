#' Plot Distribution of Non-Zero Entries in an Expression Set  
#'  
#' This function generates two plots: a bar plot showing the count of non-zero entries   
#' for each feature in the given expression matrix, and a histogram representing the   
#' distribution of the fraction of non-zero entries across samples.  
#'  
#' @param SE A SummarizedExperiment object containing the expression data.  
#' @param assayname A string indicating the assay name to use from the SummarizedExperiment object. Default is "TPM".  
#' @param ZeroasNA A logical value indicating whether zeros should be treated as NA. Default is TRUE.  
#'  
#' @return A ggplot object containing a combined plot of the bar plot and histogram.  
#'   
#' @examples  
#' # Assuming `se` is a SummarizedExperiment object  
#' distribution_plot <- SE_distribution(SE = se, assayname = "TPM", ZeroasNA = TRUE)  
#' print(distribution_plot)  
#'  
#' @import ggplot2  
#' @importFrom cowplot plot_grid  
#' @export  
SE_distribution = function(SE, assayname = "TPM", ZeroasNA = TRUE) {  
  expdata = assay(SE, assayname)  
  if (ZeroasNA) {  
      expdata[expdata == 0] <- NA  
  }  
  non_zero_count <- rowSums(!is.na(expdata))  
  plot_data <- data.frame(Feature = names(non_zero_count), Count = non_zero_count)  
  plot_data <- plot_data[order(-plot_data$Count), ]  
  plot_data$Feature <- factor(plot_data$Feature, levels = plot_data$Feature)  
  plot_data$Fraction = round(plot_data$Count / ncol(expdata), 2)  
  plot1 = ggplot(plot_data, aes(x = Feature, y = Count)) +   
    geom_col(fill = "skyblue") +   
    labs(x = paste0("Feature: ", nrow(expdata)), y = paste0("SampleCount: ", ncol(expdata)), title = "") +   
    theme_minimal() +   
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())  
  plot2 = ggplot(plot_data, aes(x = Fraction)) +   
    geom_histogram(breaks = seq(0, 1, by = 0.1), color = "gray", fill = "skyblue") +   
    labs(x = "Fraction", y = "Count", title = "") +   
    theme_minimal() + theme(panel.grid = element_blank()) +   
    scale_x_continuous(breaks = seq(0, 1, by = 0.1))  
  plot = plot_grid(plot1, plot2, ncol = 2, align = "h")  
  return(plot)  
}