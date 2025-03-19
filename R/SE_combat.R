#' Batch Effect Correction for SummarizedExperiment Objects  
#'  
#' This function applies the ComBat algorithm to correct batch effects in   
#' a SummarizedExperiment object based on a specified batch variable and optional experimental design comparison.  
#'  
#' @param SE A SummarizedExperiment object containing the assay data.  
#' @param col_for_combat A character string specifying the column name of   
#'   the batch variable in the sample data. This must be provided and valid.  
#' @param col_for_compare A character string specifying the column name of  
#'   the condition variable in the sample data for controlling comparisons.  
#'   This can be NULL, in which case it will not be used.  
#'  
#' @return A SummarizedExperiment object with the assay data corrected for   
#'   batch effects.  
#'   
#' @import SummarizedExperiment  
#' @import sva  
#' @export  

SE_combat = function(SE, col_for_combat, col_for_compare = NULL) {  
    library(SummarizedExperiment)  
    library(sva)  
    
    # Extract feature and sample information  
    feature_info = rowData(SE)  
    sample_info = colData(SE)  
    
    # Ensure col_for_combat is provided and valid  
    if (is.null(col_for_combat) || !col_for_combat %in% colnames(sample_info)) {  
        stop("Error: 'col_for_combat' must be provided and must be a valid column name in sample_info.")  
    }  
    
    # Extract batch information and convert to factor  
    batch <- sample_info[[col_for_combat]]  
    sample_info[[col_for_combat]] <- factor(sample_info[[col_for_combat]], levels = unique(sample_info[[col_for_combat]]))  
    
    # If col_for_compare is provided, ensure it is valid  
    if (!is.null(col_for_compare) && !col_for_compare %in% colnames(sample_info)) {  
        stop("Error: 'col_for_compare' must be a valid column name in sample_info if provided.")  
    }  

    # Function to perform batch effect correction  
    ExpdataBatch = function(assayname) {  
        expdata = assay(SE, assayname)  # Extract expression data for the given assay  
        
        # Perform ComBat correction  
        if (!is.null(col_for_compare)) {  
            condition <- sample_info[[col_for_compare]]  # Extract condition information  
            expdata_combat <- ComBat(dat = expdata, batch = batch, mod = model.matrix(~ condition))  # Remove batch effects using ComBat with condition  
        } else {  
            expdata_combat <- ComBat(dat = expdata, batch = batch)  # Remove batch effects using ComBat without condition  
        }  
        
        return(expdata_combat)  # Return the corrected expression data  
    }  
    
    # Initialize a list to store corrected expression data  
    expdata_combat = list()  
    
    # Loop through each assay and perform batch effect correction  
    for (assayname in assayNames(SE)) {  
        expdata_combat[[assayname]] = ExpdataBatch(assayname)  # Store the result for each assay in the list  
    }  
    
    # Create a new SummarizedExperiment object with the corrected assays  
    SEcombat <- SummarizedExperiment(assays = expdata_combat, colData = sample_info, rowData = feature_info)  
    return(SEcombat)  
}