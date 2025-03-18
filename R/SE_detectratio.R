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
    detect_samples <- rowSums(expdata != 0)                   
    total_samples <- ncol(expdata)  
    feature_info[,paste0("detectsample")] <- detect_samples  
    feature_info[,paste0("detectratio")] <- round(detect_samples / total_samples, 4)   

    # Detection ratios in groups  
    if (!group_col %in% colnames(colData(SE))) {  
        stop("The colData must contain a column named 'group'.")  
    }    
    group_info <- colData(SE)[[group_col]]     
    split_matrices <- split(as.data.frame(expdata), group_info)  

    Detect_sampleandratio = function(data, selectgroup) {  
        detect_samples <- rowSums(data != 0)                   
        total_samples <- ncol(data)  
        feature_info[,paste0("detectsample_", selectgroup)] <- detect_samples  
        feature_info[,paste0("detectratio_", selectgroup)] <- round(detect_samples / total_samples, 4)   
        return(feature_info)  
    }  

    for(selectgroup in names(split_matrices)){  
        feature_info = Detect_sampleandratio(split_matrices[[selectgroup]], selectgroup)  
    }  
    
    # Update SE's rowData with the updated feature_info  
    rowData(SE) <- feature_info  

    # Generate histogram of detection ratios  
    hist_plot <- ggplot(feature_info, aes(x = detectratio)) +  
                 geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black", alpha = 0.7) +  
                 geom_density(aes(y = ..count.. * 0.05), color = "red", size = 1) +  
                 labs(title = "Density Distribution of Detection Ratios",  
                      x = "Detection Ratio",  
                      y = "Count") +  
                 theme_minimal()  

    # Return the updated SE  
    return(list(SE = SE, plot = hist_plot))  
}