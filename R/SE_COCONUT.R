#' Apply COCONUT to SummarizedExperiment List  
#'  
#' This function applies the COCONUT (COmbat CO-Normalization Using conTrols)   
#' algorithm to a list of SummarizedExperiment objects to remove batch effects  
#' while preserving biological signals.  
#'  
#' @param SElist A list of SummarizedExperiment objects to be batch-corrected.  
#' @param assayname Character string specifying which assay in the SummarizedExperiment   
#'   objects to use for batch correction. This must be an existing assay in all   
#'   SummarizedExperiment objects.  
#' @param group_col Character string specifying the column name in colData that   
#'   contains group information (e.g., disease status).  
#' @param label_healthy Character string specifying the value in the group_col   
#'   that represents healthy controls (e.g., "HC").  
#'  
#' @return A list of SummarizedExperiment objects with batch-corrected data.  
#'  
#' @details The function uses COCONUT to perform batch correction on expression   
#'   data from multiple sources. It requires a column in colData to distinguish   
#'   between disease and healthy samples. The algorithm uses healthy controls to   
#'   learn batch effects and applies the correction to all samples.  
#'   
#'   Note that the expression matrices must not contain zeros, as this can cause  
#'   numerical issues in the COCONUT algorithm. Please impute zeros before using  
#'   this function.  
#'  
#' @examples  
#' \dontrun{  
#' # Create a list of SummarizedExperiment objects  
#' data_list <- list(SE_batch1, SE_batch2, SE_batch3)  
#'   
#' # Apply COCONUT correction  
#' corrected_data <- SE_COCONUT(  
#'   SElist = data_list,  
#'   assayname = "counts",  
#'   group_col = "condition",  
#'   label_healthy = "control"  
#' )  
#' }  
#'  
#' @import SummarizedExperiment 
#' @import COCONUT 
#' @export  
SE_COCONUT = function(SElist, assayname = NULL, group_col = "group", label_healthy = "HC") {  
    options(warn = -1)   
    message("Starting SE_COCONUT function...")  
    
    # Check if assayname is provided  
    if (is.null(assayname)) {  
        stop("assayname must exist in the original data and must be specified")  
    }  
    message("Assay name provided: ", assayname)  
    
    # Check if assayname exists in all SummarizedExperiment objects  
    for (i in seq_along(SElist)) {  
        if (!assayname %in% assayNames(SElist[[i]])) {  
            stop(paste0("assayname '", assayname, "' does not exist in dataset ", i,   
                         ". Please specify an assay that exists in all datasets."))  
        }  
        message("Assay name found in dataset ", i)  
        
        # Check if expression matrix contains zeros  
        expr_data <- assay(SElist[[i]], assayname)  
        if (any(expr_data == 0, na.rm = TRUE)) {  
            message("Zero detected in dataset ", names(SElist)[i], ". Adding 0.0001 to expression matrix.")  
            expr_data <- expr_data + 0.0001  
            assay(SElist[[i]], assayname) <- expr_data  # Update the assay with adjusted values   
        }  	  
    }  
  
    create_COCOobj_type <- function(SEinput, assayname, group_col, label_healthy) {  
        expdata = as.data.frame(assay(SEinput, assayname))  
		sample_info = colData(SEinput)  
		sample_info[] <- lapply(sample_info, function(x) {  
		if (inherits(x, "integer64")) {  
			return(as.integer(x))  # integer64 to numeric  
		} else {  
			return(x)  
		}  
		sample_info = as.data.frame(sample_info)  
	})   
        sample_info[sample_info[,group_col] != label_healthy,]$group = 1   
        sample_info[sample_info[,group_col] == label_healthy,]$group = 0  
        return(list(pheno = sample_info, genes = expdata))  
    }  
    
    COCOobjs = list()  
    for (i in seq_along(SElist)) {  
        COCOobjs[[i]] <- create_COCOobj_type(SElist[[i]], assayname = assayname,   
                                             group_col = group_col, label_healthy = label_healthy)  
        message("Created COCO object for dataset ", i)  
    }  
    
    GSEs.COCONUT <- COCONUT(GSEs = COCOobjs, control.0.col = group_col, byPlatform = FALSE)  
    message("COCONUT correction applied.")  
    
    expdata_list = list()  
    feature_list = list()  
    sample_list = list()  
    
    for (j in seq_along(SElist)) {  
        expdata_batch1 = GSEs.COCONUT[["COCONUTList"]][[j]][["genes"]]  
        expdata_batch0 = GSEs.COCONUT$controlList$GSEs[[j]]$genes  
        expdata_batch = cbind(expdata_batch0, expdata_batch1)  
        
        # Ensure correct ordering  
        expdata_list[[j]] = expdata_batch  
        
        # Extract feature info and match with expression data rows  
        feature_info = rowData(SElist[[j]])  
        feature_info = feature_info[rownames(feature_info) %in% rownames(expdata_batch),]  
        feature_info = feature_info[match(rownames(expdata_batch), rownames(feature_info)),]  
        feature_list[[j]] = feature_info  
        
        # Extract sample info and match with expression data columns  
        sample_info = colData(SElist[[j]])  
        sample_info = sample_info[rownames(sample_info) %in% colnames(expdata_batch),]  
        sample_info = sample_info[match(colnames(expdata_batch), rownames(sample_info)),]  
        sample_list[[j]] = sample_info  
    }  
    
    SElist_batch = list()  
    for(k in seq_along(SElist)){  
        assays_list <- list()  
        assays_list[[assayname]] <- as.matrix(expdata_list[[k]])  
        
        SEbatch <- SummarizedExperiment(  
            assays = assays_list,  
            rowData = feature_list[[k]],  
            colData = sample_list[[k]]  
        )  
        
        SElist_batch[[k]] = SEbatch  
        message("Created SummarizedExperiment object for dataset ", k)  
    }   
    names(SElist_batch) = names(SElist)  
    message("Completed SE_COCONUT function.")  
    
    return(SElist_batch)  
}