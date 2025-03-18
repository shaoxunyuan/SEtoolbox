#' Calculate Detection Ratio and Update SummarizedExperiment Object's rowData  
#'  
#' This function computes the detection ratio of expression data and updates  
#' the rowData of the provided SummarizedExperiment object with detection sample counts  
#' and ratios. It also generates a histogram of detection ratios.  
#'  
#' @param SE A SummarizedExperiment object containing expression data.  
#' @param assayname The name of the assay to be used for calculations. Default is "TPM".  
#' @param group_col The name of the column containing group information. Default is "group".  
#' @return An updated SummarizedExperiment object with detection sample counts  
#'   and ratios in rowData, along with a displayed histogram of detection ratios.  
#' @import SummarizedExperiment  
#' @import ggplot2  
#' @export  
SE_detectratio <- function(SE, assayname = "TPM", group_col = "group") {  
    # Load necessary libraries  
    options(warn = -1)  
    library(SummarizedExperiment)  
    library(ggplot2)  
    
    # Check input validity  
    if (!inherits(SE, "SummarizedExperiment")) {  
        stop("Input SE must be a SummarizedExperiment object.")  
    }  
    if (!assayname %in% names(assays(SE))) {  
        stop(paste("Assay", assayname, "not found in SE."))  
    }  
    
    # Get feature information and expression matrix  
    feature_info <- rowData(SE)  
    expdata <- assay(SE, assayname)  
    
    # Calculate detection samples and ratios  
    detect_samples <- rowSums(expdata != 0)  
    total_samples <- ncol(expdata)  
    
    # Update feature_info with detection sample counts and ratios  
    feature_info[["detectsample"]] <- detect_samples  
    feature_info[["detectratio"]] <- round(detect_samples / total_samples, 4)  

    # Ensure the group column exists  
    if (!group_col %in% colnames(colData(SE))) {  
        stop("The colData must contain a column named 'group'.")  
    }    
    group_info <- colData(SE)[[group_col]]  
    
    # Split expression data by group  
    split_matrices <- split(as.data.frame(expdata), group_info)  
    
    # Function to calculate detection samples and ratios per group  
    Detect_sampleandratio <- function(data, selectgroup) {  
        detect_samples <- rowSums(data != 0)  
        total_samples <- ncol(data)  
        
        feature_info[[paste0("detectsample_", selectgroup)]] <- detect_samples  
        feature_info[[paste0("detectratio_", selectgroup)]] <- round(detect_samples / total_samples, 4)  
        
        return(feature_info)  
    }  

    # Calculate detection ratios for each group  
    for(selectgroup in names(split_matrices)) {  
        feature_info <- Detect_sampleandratio(split_matrices[[selectgroup]], selectgroup)  
    }  

    # Update SE's rowData with the updated feature_info  
    rowData(SE) <- feature_info  

    # Create a plot of detection ratios  
    non_zero_count <- rowSums(!is.na(expdata))  
    plot_data <- data.frame(Feature = names(non_zero_count), Count = non_zero_count)  
    plot_data <- plot_data[order(-plot_data$Count), ]  
    plot_data$Feature <- factor(plot_data$Feature, levels = plot_data$Feature)  
    plot_data$Fraction <- round(plot_data$Count / ncol(expdata), 2)  

    # Plot sample counts  
    plot1 <- ggplot(plot_data, aes(x = Feature, y = Count)) +   
        geom_col(fill = "steelblue") +   
        labs(x = paste0("Feature: ", nrow(expdata)), y = paste0("Sample Count: ", ncol(expdata)), title = "") +      
        theme_minimal() +   
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(),   
              plot.background = element_rect(fill = "white", color = NA))  

    # Plot histogram of detection ratios  
    plot2 <- ggplot(as.data.frame(rowData(SE)), aes(x = detectratio)) +  
        geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black", alpha = 0.7) +  
        geom_density(aes(y = ..count.. * 0.05), color = "salmon", size = 1) +  
        labs(title = "",  x = "Detection Ratio",  y = "Count") +  
        scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +  
        theme_minimal() +  
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())  

    # Combine plots  
    plot <- plot_grid(plot1, plot2, ncol = 2, align = "h")  

    return(list(SE = SE, plot = plot))  
}