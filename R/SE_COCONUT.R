#' @title Create a SummarizedExperiment object for COCONUT analysis  
#'  
#' @description  
#' This function processes a list of SummarizedExperiment objects, extracts  
#' the assay data and corresponding sample information, and prepares a new  
#' SummarizedExperiment object suitable for COCONUT analysis. It assigns a  
#' binary label for healthy samples based on the specified group column.  
#'  
#' @param SElist A list of SummarizedExperiment objects. Each object should contain  
#'               the assay data and the sample information needed for analysis.  
#' @param assayname A string specifying the name of the assay to extract from the  
#'                  SummarizedExperiment objects (default is "TPM").  
#' @param group_col A string specifying the name of the column in the colData that   
#'                  contains group information (default is "group").  
#' @param label_healthy A string used to identify healthy samples in the group column   
#'                      (default is "HC").  
#'   
#' @return A SummarizedExperiment object containing the combined assay data,  
#'         feature information, and merged sample information.  
#'  
#' @examples  
#' # Load required package  
#' library(SummarizedExperiment)  
#'   
#' # Create a list of example SummarizedExperiment objects  
#' SE1 <- SummarizedExperiment(assays = list(TPM = matrix(rnorm(100), nrow=10)),  
#'                             colData = DataFrame(group = rep(c("HC", "disease"), each = 5)))  
#' SE2 <- SummarizedExperiment(assays = list(TPM = matrix(rnorm(100), nrow=10)),  
#'                             colData = DataFrame(group = rep(c("HC", "disease"), each = 5)))  
#' SElist <- list(SE1, SE2)  
#'   
#' # Create the combined SummarizedExperiment object for COCONUT analysis  
#' SE_COCONUT_result <- SE_COCONUT(SElist, assayname = "TPM", group_col = "group", label_healthy = "HC")  
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