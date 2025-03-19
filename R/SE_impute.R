#' Fill Missing Values in a SummarizedExperiment Object  
#'  
#' This function imputes missing values (NA) in the provided SummarizedExperiment object using specified methods.  
#' Various imputation techniques can be used to handle missing values, ensuring the robustness of subsequent analyses.  
#'  
#' @param object A SummarizedExperiment object containing the data to be imputed.  
#' @param assayname The name of the assay in the SummarizedExperiment, specifying the type of data to be imputed.  
#' @param group A character string specifying the grouping variable in the sample data.  
#' @param ZerosAsNA A logical value indicating whether to treat zeros as NA. Default is FALSE.  
#' @param RemoveNA A logical value indicating whether to remove samples with a high percentage of NA values based on the cutoff. Default is TRUE.  
#' @param cutoff A numeric value representing the percentage cutoff for NA samples. Default is 20.  
#' @param method A character string specifying the imputation method to use. Options include "none", "LOD", "half_min",  
#'   "median", "mean", "min", "knn", "rf", "global_mean", "svd", and "QRILC". Default is "none".  
#' @param LOD A numeric value representing the limit of detection (used for the LOD imputation method). Default is NULL.  
#' @param knum An integer representing the number of neighbors in the KNN imputation method. Default is 10.  
#'  
#' @return A SummarizedExperiment object containing the imputed values.  
#'  
#' @examples  
#' # Example usage  
#' # Assuming `se` is a SummarizedExperiment object  
#' se_imputed <- SE_impute(se, assayname = "TPM", group = "my_group_column", method = "median", ZerosAsNA = TRUE)  
#'  
#' @export  
SE_impute <- function(SE, assayname = "TPM", group = "group", ZerosAsNA = FALSE, RemoveNA = TRUE,  
                      cutoff = 20, method = c("none", "LOD", "half_min", "median", "mean", "min", "knn", "rf", "global_mean", "svd", "QRILC"),  
                      LOD = NULL, knum = 10) {  
    library(plyr)  
    library(dplyr)  
    
    # Match method parameter  
    method <- match.arg(method, c("none", "LOD", "half_min", "median", "mean", "min", "knn", "rf", "global_mean", "svd", "QRILC"))  
    
    # Process the SummarizedExperiment object  
    sam_tab <- SummarizedExperiment::colData(SE)   
    sam_tab[] <- lapply(sam_tab, function(x) {   
        if(inherits(x, "integer64")) {   
            return(as.numeric(x))  
        }  
        return(x)  
    })   
    sam_tab <- as.data.frame(sam_tab) %>% tibble::rownames_to_column("TempRowNames")   
    prf_tab <- SummarizedExperiment::assay(SE, assayname) %>% as.data.frame() %>% t()   

    # Find index for sample grouping  
    group_index <- which(colnames(sam_tab) == group)  
    samples_groups <- sam_tab[, group_index]  
    to_imp_data <- prf_tab %>% as.matrix()  

    # Treat zeros as NA if set to TRUE  
    if (ZerosAsNA) {  
        to_imp_data[to_imp_data == 0] <- NA  
    }  
    
    to_imp_data <- data.frame(Group = samples_groups, to_imp_data)  
    colnames(to_imp_data)[2:ncol(to_imp_data)] <- colnames(prf_tab)  

    # Calculate the proportion of NAs in the data  
    percent_na <- sum(is.na(to_imp_data))  
    if (percent_na == 0) {  
        print("No missing values detected")  
        if (method != "none") {  
            method <- "none"  
        }  
    }  

    # Remove samples with high percentages of missing values if needed  
    if (isTRUE(RemoveNA)) {  
        count_NA <- stats::aggregate(. ~ Group, data = to_imp_data,   
            function(x) {  
                100 * (sum(is.na(x)) / (sum(is.na(x)) + sum(!is.na(x))))  
            }, na.action = NULL)  
        count_NA <- count_NA %>% dplyr::select(-Group)  
        supress <- unlist(as.data.frame(lapply(count_NA, function(x) all(x > cutoff))))  
        correct_names <- names(count_NA)  
        names(supress) <- correct_names  
        correct_names <- names(supress[supress == "FALSE"])  
        depurdata <- to_imp_data[, 2:ncol(to_imp_data)][!supress]  
        rownames(depurdata) <- rownames(to_imp_data)[!supress]  
    } else {  
        depurdata <- to_imp_data[, 2:ncol(to_imp_data)]  
        rownames(depurdata) <- rownames(to_imp_data)  
    }  

    # Impute missing values based on the selected method  
    if (method == "none") {  
        depurdata[is.na(depurdata)] <- 0  
    } else if (method == "LOD") {  
        if (is.null(LOD)) {  
            message("LOD not provided; using one-tenth of the minimum value as LOD")  
            depurdata_withoutNA <- depurdata[!is.na(depurdata)]  
            LOD <- min(depvar(data_no_na)[depurdata_withoutNA != 0]) / 10  
        }  
        depurdata[is.na(depurdata)] <- LOD  
        depurdata[depurdata == 0] <- LOD  
    } else if (method == "half_min") {  
        depurdata <- apply(depurdata, 2, function(x) {   
            if (is.numeric(x))   
                ifelse(is.na(x), min(x, na.rm = TRUE) / 2, x)   
            else x  
        })  
    } else if (method == "median") {  
        depurdata <- apply(depurdata, 2, function(x) {  
            if (is.numeric(x))   
                ifelse(is.na(x), median(x, na.rm = TRUE), x)  
            else x  
        })  
    } else if (method == "mean") {  
        depurdata <- apply(depurdata, 2, function(x) {  
            if (is.numeric(x))   
                ifelse(is.na(x), mean(x, na.rm = TRUE), x)  
            else x  
        })  
    } else if (method == "min") {  
        depurdata <- apply(depurdata, 2, function(x) {  
            if (is.numeric(x))   
                ifelse(is.na(x), min(x, na.rm = TRUE), x)  
            else x  
        })  
    } else if (method == "knn") {  
        depurdata <- t(depurdata)  
        datai <- impute::impute.knn(depurdata, k = knum, rowmax = 0.5,   
            colmax = 0.8, maxp = 1500)  
        depurdata <- t(datai$data)  
    } else if (method == "rf") {  
        fit <- missForest::missForest(t(depurdata))  
        depurdata <- fit$ximp %>% t()  
    } else if (method == "global_mean") {  
        depurdata <- .GlobalMean(SE = t(depurdata)) %>% t()  
    } else if (method == "svd") {  
        depurdata <- .SVD_wrapper(depurdata)  
    } else if (method == "QRILC") {  
        fit <- log(t(depurdata)) %>% imputeLCMD::impute.QRILC()  
        depurdata <- t(fit[[1]])  
    }  

    colnames(depurdata) <- colnames(prf_tab)  
    rownames(depurdata) <- rownames(prf_tab)  

    # Return final result   
    if (ncol(depurdata) != ncol(prf_tab)) {  
        rdata <- SummarizedExperiment::rowData(SE)  
        cdata <- SummarizedExperiment::colData(SE)  
        mdata <- if (length(SE@metadata) == 0) { NULL } else { SE@metadata }  
        rdata = rdata[rownames(rdata) %in% colnames(depurdata),]  
        res <- SummarizedExperiment(assays = t(depurdata), rowData = rdata, colData = cdata, metadata = mdata)  
    } else {  
        res <- SE  
        SummarizedExperiment::assay(res) <- t(depurdata)  
    }    

    return(res)  
}