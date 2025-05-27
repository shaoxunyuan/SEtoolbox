#' @title SE_detectratio: Calculate Detection Ratio for SummarizedExperiment Object  
#'   
#' @description   
#' This function computes the detection ratio of expression data in a   
#' \code{SummarizedExperiment} object and updates its \code{rowData} with detection   
#' sample counts and ratios. It also generates a histogram of detection ratios to   
#' visualize the distribution. Detection ratios are calculated based on the number   
#' of non-zero samples for each feature across provided groups. The function can also   
#' calculate detection ratios per group if specified.  
#'   
#' @param SE A \code{SummarizedExperiment} object containing gene expression data. This   
#' should include the assay for which detection ratios will be calculated.  
#' @param assayname A string indicating which assay to use for calculations.   
#' The default value is \code{"TPM"}. The assay must exist in the \code{SummarizedExperiment} object.  
#' @param group_colname A string representing the column name in \code{colData} that   
#' contains group information. The default value is \code{"group"}. This is optional; if not provided,   
#' the function will calculate overall detection ratios only.  
#'   
#' @return   
#' A list containing:  
#' \item{SE}{An updated \code{SummarizedExperiment} object with detection sample counts  
#' and ratios added to \code{rowData}. The detection sample counts are added to the   
#' \code{"detectsample"} and corresponding ratios to the \code{"detectratio"} columns.}  
#' \item{plot_feature}{A \code{ggplot} object representing the histogram of detection ratios   
#' for features. This visualizes the distribution of detection ratios across all features   
#' and groups if applicable.}  
#' \item{plot_sample}{A \code{ggplot} object representing the bar plot of detection ratios   
#' for samples. The Y-axis displays the expression fraction formatted as a percentage,   
#' and the samples are displayed on the X-axis.}  
#'
#' @import SummarizedExperiment  
#' @import dplyr  
#' @import ggplot2  
#' @import tidyr  
#' @import cowplot  
#' @import RColorBrewer  
#' @import scales   
#' 
#' @examples   
#' # Load example SummarizedExperiment object 
#' 
#' SE <- loadSE()   
#'
#' # Call the SE_detectratio function
#'  
#' result <- SE_detectratio(SE, assayname = "TPM", group_colname = "group")   
#'
#' SEdectectratio <-  result$SE # Return an SummarizedExperiment object after detectratio calculate 
#'  
#' print(result$plot_sample)  # Display the sample histogram  
#'
#' print(result$plot_feature)  # Display the sample histogram
#'  
#' @export  

SE_detectratio <- function(SE, assayname = "TPM", group_colname = NULL) {  
    # Check input validity  
    if (!inherits(SE, "SummarizedExperiment")) {  
        stop("Input SE must be a SummarizedExperiment object.")  
    }  
    if (!assayname %in% names(assays(SE))) {  
        stop(paste("Assay", assayname, "not found in SE."))  
    }  
    
    feature_info <- as.data.frame(rowData(SE))  
    sample_info <- as.data.frame(colData(SE))  
    expdata <- as.data.frame(assay(SE, assayname))  
    
	# zero express plot
	zero_counts <- rowSums(expdata == 0)  
	zero_counts_df <- data.frame(Column = names(zero_counts), Zero_Counts = zero_counts)
	zero_counts_df <- zero_counts_df[order(zero_counts_df$Zero_Counts), ]
	zero_counts_df$Column <- factor(zero_counts_df$Column,   
									 levels = zero_counts_df$Column[order(zero_counts_df$Zero_Counts,decreasing = T)])  
	plot_feature_distribution = ggplot(zero_counts_df, aes(x = Column, y = Zero_Counts)) +  
	  geom_bar(stat = "identity") +  
	  labs(title = paste0("Zero express counts (",ncol(expdata)," samples,",nrow(expdata)," features)"),  
		   x = "Genes",  
		   y = "Zero express counts") +  
	  theme(axis.text.x = element_blank(),
			axis.text.y = element_text(size=12),
			axis.ticks.x = element_blank(), 
			axis.line.x = element_blank()) 
    # Calculate detection samples and ratios  
    detect_samples <- rowSums(expdata != 0)  
    total_samples <- ncol(expdata)  
    feature_info[,"detectsample"] <- detect_samples  
    feature_info[,"detectratio"] <- round(detect_samples / total_samples, 4)  
    
    # Group calculations  
    if (!is.null(group_colname) && group_colname %in% colnames(sample_info)) {  
        group_info <- unique(sample_info[[group_colname]])  
        for(selectgroup in group_info) {  
            expdata_sub <- expdata[, sample_info[, group_colname] == selectgroup, drop = FALSE]  
            detect_samples <- rowSums(expdata_sub != 0)  
            total_samples <- ncol(expdata_sub)  
            feature_info[, paste0("detectsample_", selectgroup)] <- detect_samples  
            feature_info[, paste0("detectratio_", selectgroup)] <- round(detect_samples / total_samples, 4)  
        }  
    }  
    
    rowData(SE) <- feature_info   
    
    express_counts <- colSums(expdata != 0)  
    express_counts_df <- data.frame(ExpressCount = express_counts, ExpressFraction = round(express_counts / nrow(expdata), 4))  
    
    sample_info$ExpressCount = plyr::mapvalues(rownames(sample_info),rownames(express_counts_df),express_counts_df$ExpressCount,warn_missing = F)
	
	sample_info$ExpressFraction = plyr::mapvalues(rownames(sample_info),rownames(express_counts_df),express_counts_df$ExpressFraction,warn_missing = F)
	
	sample_info$ExpressCount = as.numeric(sample_info$ExpressCount)
	
	sample_info$ExpressFraction = as.numeric(sample_info$ExpressFraction)
	
    colData(SE) <- S4Vectors::DataFrame(sample_info)   
    
    # Histogram of detection ratios  
    detectratio_long <- feature_info %>% select(starts_with("detectratio")) %>% pivot_longer(everything(), names_to = "Group", values_to = "value")  
    
    total <- detectratio_long[detectratio_long$Group == "detectratio", ]  
    plot_feature_fraction <- ggplot(total, aes(x = value)) +  
        geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black", alpha = 0.7) +  
        geom_density(aes(y = after_stat(count) * 0.1), color = "salmon", size = 1) +  
        labs(title = "", x = "Detection Ratio", y = "Count") +  
        theme_minimal()  
  
    # Group comparisons  
    #group <- detectratio_long[!detectratio_long$Group == "detectratio", ]  
    #if (nrow(group) > 0) {  
    #    plot_feature_group = ggplot(group, aes(x = value, fill = Group)) +  
    #        geom_histogram(aes(y = after_stat(count)), position = "identity", binwidth = 0.1, color = "black", alpha = 0.5) +  
    #        geom_density(aes(y = after_stat(count) * 0.1, color = Group), size = 0.5, alpha = 0.5) +  
    #        labs(title = "", x = "Detection Ratio", y = "Count") +  
    #        scale_fill_manual(values = brewer.pal(8, "Set2")) +   
    #        scale_color_manual(values = brewer.pal(8, "Set2")) + 
    #        theme_minimal()  
	#} 
    
    # Sample expression plot  
    plot_sample <- ggplot(sample_info, aes(x = reorder(BioSample, ExpressCount), y = ExpressCount, color = group)) +  
				geom_point(size = 1, shape = 21, fill = "white") +  
				geom_smooth(method = "loess", color = "black", size = 1.2) + 
				labs(title = "", x = "Sample", y = "Expression Count") +  
				scale_y_continuous(breaks = pretty(range(sample_info$ExpressCount), n = 10)) + 
				theme_minimal() +   
				theme(axis.text.x = element_blank(),          
						axis.text.y = element_text(size = 10),                                  
						plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
						panel.grid.major = element_blank(),                                 
						panel.grid.minor = element_blank(),
						panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
						legend.position = c(0.95, 0.05), 
						legend.justification = c(1, 0)) 


		return(list(SE = SE, 
		plot_feature_fraction = plot_feature_fraction, 
		plot_feature_distribution = plot_feature_distribution, data_feature_distribution = zero_counts_df,
		plot_sample = plot_sample))  
}