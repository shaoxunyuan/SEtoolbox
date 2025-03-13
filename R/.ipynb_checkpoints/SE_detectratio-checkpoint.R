#' Calculate Detection Ratio and Update SummarizedExperiment Object's rowData  
#'  
#' This function computes the detection ratio of expression data and updates   
#' the rowData of the provided SummarizedExperiment object with detection sample counts   
#' and ratios. It also generates a histogram of detection ratios.  
#'  
#' @param SE A SummarizedExperiment object containing expression data.  
#' @param assayname The name of the assay to be used for calculations. Default is "TPM".  
#' @return An updated SummarizedExperiment object with detection sample counts   
#'   and ratios in rowData, and a histogram of detection ratios is displayed.  
#' @import SummarizedExperiment  
#' @import ggplot2  
#' @export  
SE_detectratio <- function(SE, assayname = "TPM") {  
    # Load necessary libraries  
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

    # Calculate number of non-zero samples for each feature  
    detectsample_counts <- rowSums(expdata != 0)  

    # Update the detectsample column in feature_info  
    feature_info$detectsample <- detectsample_counts  

    # Calculate detection ratio  
    total_samples <- ncol(expdata)  
    feature_info$detectratio <- round(feature_info$detectsample / total_samples, 4)   

    # Update SE's rowData with the updated feature_info  
    rowData(SE) <- feature_info  

    # Generate histogram of detection ratios  
    hist_plot <- ggplot(feature_info, aes(x = detectratio)) +  
                 geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +  
                 geom_density(aes(y = ..count.. * 0.05), color = "red", size = 1) +  
                 labs(title = "Density Distribution of Detection Ratios",  
                      x = "Detection Ratio",  
                      y = "Count") +  
                 theme_minimal()  

    # Print the histogram  
    print(hist_plot)  

    # Return the updated SE  
    return(SE)  
}  