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
    
	# Get feature/sample information and expression matrix  
    feature_info <- as.data.frame(rowData(SE)  )
    sample_info <- as.data.frame(colData(SE))
	sample_info[] <- lapply(sample_info, function(x) {  
    if (inherits(x, "integer64")) {  
        return(as.integer(x))  # 或者根据需要使用as.numeric(x)  
    } else {  
        return(x)  
    }  
})  
    expdata <- as.data.frame(assay(SE, assayname))
    
    # Calculate detection samples and ratios  
    detect_samples <- rowSums(expdata != 0)  
    total_samples <- ncol(expdata)  
    
    # Update feature_info with detection sample counts and ratios  
    feature_info[,"detectsample"] <- detect_samples  
    feature_info[,"detectratio"] <- round(detect_samples / total_samples, 4)  

    # Ensure the group column exists  
    if (!group_col %in% colnames(colData(SE))) {  
        stop("The colData must contain a column named 'group'.")  
    }    
    
    group_info <- unique(colData(SE)[[group_col]]);group_info
	
	split_matrices = list()
	for(selectgroup in group_info){
		expdata_sub = expdata[,colnames(expdata) %in% rownames(sample_info[sample_info[,group_col] == selectgroup,])]
		split_matrices[[selectgroup]] = expdata_sub
	}
	
	# Function to calculate detection samples and ratios per group  
    Detect_sampleandratio <- function(data, selectgroup) {  
        detect_samples <- rowSums(data != 0)  
        total_samples <- ncol(data)  
        
        feature_info[,paste0("detectsample_", selectgroup)] <- detect_samples  
        feature_info[,paste0("detectratio_", selectgroup)] <- round(detect_samples / total_samples, 4)  
        
        return(feature_info)  
    }  

    # Calculate detection ratios for each group  
    for(selectgroup in group_info) {  
        feature_info <- Detect_sampleandratio(split_matrices[[selectgroup]], selectgroup)  
    }  

    # Update SE's rowData with the updated feature_info  
    rowData(SE) <- feature_info   
	
    # Plot histogram of detection ratios  
    plot <- ggplot(as.data.frame(rowData(SE)), aes(x = detectratio)) +  
        geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black", alpha = 0.7) +  
        geom_density(aes(y = ..count.. * 0.05), color = "salmon", size = 1) +  
        labs(title = "",  x = "Detection Ratio",  y = "Count") +  
        scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +  
        theme_minimal() +  
        theme(panel.grid = element_blank())  

    return(list(SE = SE, plot = plot))  
}