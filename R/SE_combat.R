#' Perform Batch Effect Correction using ComBat  
#'  
#' This function adjusts for batch effects in a SummarizedExperiment object using the ComBat method.  
#'  
#' @param SE A SummarizedExperiment object containing expression data.  
#' @param col_for_combat A string specifying the column name in the sample information that contains batch information.  
#' @param col_for_compare An optional string specifying the column name in the sample information that contains the condition for comparison.  
#' import sva
#' @return A SummarizedExperiment object with batch effects corrected.  
#'   
#' @examples  
#' # Assuming 'se' is a SummarizedExperiment object  
#' corrected_se <- SE_combat(se, col_for_combat = "batch", col_for_compare = "condition")  
#'  
#' @export
SE_combat = function(SE, col_for_combat, col_for_compare = NULL) {  	 
    
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

    # Check for col_for_compare correlation with batch  
    if (!is.null(col_for_compare)) {  
        compare <- sample_info[[col_for_compare]]  # Extract compare information  
        
        # Create a contingency table and perform a chi-squared test  
        contingency_table <- table(batch, compare)  
        chisq_test <- chisq.test(contingency_table)  

        # If p-value is significant, indicating correlation  
        if (chisq_test$p.value < 0.05) {  
            message("Warning: compare information is strongly correlated with the batch. Please consider removing the compare information before batch effect correction.")  
            return(NULL)  # Stop further execution and return NULL  
        }  
    }  

    # Function to perform batch effect correction  
    ExpdataBatch = function(assayname) {  
        expdata = assay(SE, assayname)  # Extract expression data for the given assay  
        
        # Perform ComBat correction  
        if (!is.null(col_for_compare)) {  
            compare <- sample_info[[col_for_compare]]  # Extract compare information  
            expdata_combat <- ComBat(dat = expdata, batch = batch, mod = model.matrix(~ compare))  # Remove batch effects using ComBat with compare  
        } else {  
            expdata_combat <- ComBat(dat = expdata, batch = batch)  # Remove batch effects using ComBat without compare  
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