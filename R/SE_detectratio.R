#' @title SE_detectratio: Calculate Detection Ratio for SummarizedExperiment Object  
#'   
#' @description   
#' This function computes the detection ratio of expression data in a   
#' \code{SummarizedExperiment} object and updates its \code{rowData} with detection   
#' sample counts and ratios. It also generates a histogram of detection ratios to   
#' visualize the distribution. Detection ratios are calculated based on the number   
#' of non-zero samples for each feature across provided groups.  
#'   
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.  
#' @param assayname A string indicating which assay to use for calculations.   
#' The default value is \code{"TPM"}.  
#' @param group_col A string representing the column name in \code{colData} that   
#' contains group information. The default value is \code{"group"}.  
#'   
#' @return   
#' A list containing:  
#' \item{SE}{An updated \code{SummarizedExperiment} object with detection sample counts  
#' and ratios added to \code{rowData}.}  
#' \item{plot}{A \code{ggplot} object representing the histogram of detection ratios.}  
#'   
#' @examples   
#' # Create a dummy SummarizedExperiment object  
#' data_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)  
#' rownames(data_matrix) <- paste0("Gene", 1:100)  
#' colnames(data_matrix) <- paste0("Sample", 1:10)  
#' sample_info <- DataFrame(group = rep(c("A", "B"), each = 5))  
#' SE <- SummarizedExperiment(assays = list(TPM = data_matrix), colData = sample_info)  
#'   
#' # Call the SE_detectratio function  
#' result <- SE_detectratio(SE, assayname = "TPM", group_col = "group")  
#' print(result$plot)  # Display the histogram  
#'   
#' @export   
SE_detectratio <- function(SE, assayname = "TPM", group_col = NULL) {  
    
    # Check input validity  
    if (!inherits(SE, "SummarizedExperiment")) {  
        stop("Input SE must be a SummarizedExperiment object.")  
    }  
    if (!assayname %in% names(assays(SE))) {  
        stop(paste("Assay", assayname, "not found in SE."))  
    }  
    
    # Get feature/sample information and expression matrix  
    feature_info <- as.data.frame(rowData(SE))  
    sample_info <- colData(SE)  
    sample_info[] <- lapply(sample_info, function(x) {  
        if (inherits(x, "integer64")) {  
            return(as.integer(x))  # integer64 to numeric  
        } else {  
            return(x)  
        }  
    })  
    sample_info = as.data.frame(sample_info)  
    expdata <- as.data.frame(assay(SE, assayname))  
    
    # Calculate detection samples and ratios  
    detect_samples <- rowSums(expdata != 0)  
    total_samples <- ncol(expdata)  
    
    # Update feature_info with detection sample counts and ratios  
    feature_info[,"detectsample"] <- detect_samples  
    feature_info[,"detectratio"] <- round(detect_samples / total_samples, 4)  

    # If group_col is not NULL, calculate detection ratios for each group  
    if (!is.null(group_col) && group_col %in% colnames(colData(SE))) {  
        group_info <- unique(colData(SE)[[group_col]])  
        
        split_matrices = list()  
        for(selectgroup in group_info) {  
            expdata_sub = expdata[, colnames(expdata) %in% rownames(sample_info[sample_info[, group_col] == selectgroup,])]  
            split_matrices[[selectgroup]] = expdata_sub  
        }  
        
        # Function to calculate detection samples and ratios per group  
        Detect_sampleandratio <- function(data, selectgroup) {  
            detect_samples <- rowSums(data != 0)  
            total_samples <- ncol(data)  
            
            feature_info[, paste0("detectsample_", selectgroup)] <- detect_samples  
            feature_info[, paste0("detectratio_", selectgroup)] <- round(detect_samples / total_samples, 4)  
            
            return(feature_info)  
        }  

        # Calculate detection ratios for each group  
        for(selectgroup in group_info) {  
            feature_info <- Detect_sampleandratio(split_matrices[[selectgroup]], selectgroup)  
        }  
    }  
    
    # Update SE's rowData with the updated feature_info  
    rowData(SE) <- feature_info   
    
    # Plot histogram of detection ratios  
    plot <- ggplot(as.data.frame(rowData(SE)), aes(x = detectratio)) +  
        geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black", alpha = 0.7) +  
        geom_density(aes(y = ..count.. * 0.1), color = "salmon", size = 1) +  
        labs(title = "", x = "Detection Ratio", y = "Count") +  
        scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +  # 添加简单的刻度  
        theme_minimal() +  
        theme(panel.grid = element_blank())  

    return(list(SE = SE, plot = plot))  
}