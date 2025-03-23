#' Apply COCONUT method for batch effect removal on SummarizedExperiment objects
#'
#' This function takes a list of SummarizedExperiment objects and applies the COCONUT method
#' to remove batch effects. It creates COCO objects from the input SummarizedExperiment objects,
#' processes them using COCONUT, and then combines the processed expression data to create a
#' new SummarizedExperiment object with batch - corrected data.
#'
#' @param SElist A list of SummarizedExperiment objects that contain expression data and sample information.
#' @param assayname The name of the assay in the SummarizedExperiment objects to be used for analysis. Default is "TPM".
#' @param group_col The column name in the sample information that indicates the group of each sample. Default is "group".
#' @param label_healthy The label in the `group_col` that represents the healthy samples. Default is "HC".
#'
#' @return A new SummarizedExperiment object with batch - corrected expression data.
#'
#' @import SummarizedExperiment
#' @importFrom some_package COCONUT # Replace some_package with the actual package name
#'
#' @examples
#' # Assume selist is a list of SummarizedExperiment objects
#' # SE_corrected <- SE_COCONUT(selist)
#'
#' @export
SE_COCONUT = function(SElist, assayname = "TPM", group_col = "group", label_healthy = "HC") {  
    create_COCOobj_type <- function(SEinput, assayname, group_col, label_healthy) {  
        expdata = as.data.frame(assay(SEinput, assayname))  
        sample_info = as.data.frame(colData(SEinput))  
        sample_info[sample_info[[group_col]] == label_healthy, "Healthy0"] <- 0  
        sample_info[sample_info[[group_col]] != label_healthy, "Healthy0"] <- 1  
        return(list(pheno = sample_info, genes = expdata))  
    }  
    
    COCOobjs = list()  
    for (i in seq_along(SElist)) {  
        COCOobjs[[i]] <- create_COCOobj_type(SElist[[i]], assayname = assayname, group_col = group_col, label_healthy = label_healthy)  
    }  
    
    GSEs.COCONUT <- COCONUT(GSEs = COCOobjs, control.0.col = "Healthy0", byPlatform = FALSE)  
    
    expdata_batch.list = list()  
    for (j in seq_along(SElist)) {  
        expdata_batch = GSEs.COCONUT[["COCONUTList"]][[j]][["genes"]]  
        expdata_batch.list[[j]] = expdata_batch  
    }  
    expdata_batch = do.call(cbind, expdata_batch.list)  
    
    feature_info = data.frame(feature = rownames(expdata_batch))  
    rownames(feature_info) = feature_info$feature  
    
    sample_info_merge <- data.frame()  
    for (k in seq_along(SElist)) {  
        sample_info = as.data.frame(colData(SElist[[k]]))  
        sample_info = sample_info[, group_col, drop = FALSE]  # Only extract group column  
        sample_info_merge = rbind(sample_info_merge, sample_info)  
    }  
    colnames(sample_info_merge) = group_col  
    
    sample_info_merge = sample_info_merge[rownames(sample_info_merge) %in% colnames(expdata_batch), , drop = FALSE]  
    
    SE <- SummarizedExperiment(assays = list(TPM = expdata_batch), rowData = feature_info, colData = sample_info_merge)  
    return(SE)  
} 