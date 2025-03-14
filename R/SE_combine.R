SE_combine <- function(se_list, merge_type = "intersection") {  
    all_assay_names <- unique(unlist(lapply(se_list, function(se) assayNames(se))))  
    cat("assayNames of all SE obj:", all_assay_names, "\n")  
                                        
    Merge_by_intersection = function(se_list, assay_name) {  
        merged_data <- NULL   
        expdata1 <- assay(se_list[[1]], assay_name)  
        common_rows <- rownames(expdata1)  
        for (se in se_list) {  
            expdata <- assay(se, assay_name)  
            common_rows <- intersect(common_rows, rownames(expdata))  
        }  
        for (se in se_list) {  
            expdata <- assay(se, assay_name)  
            filtered_expdata <- expdata[common_rows, , drop = FALSE]  
            if (is.null(merged_data)) {  
                merged_data <- filtered_expdata  
            } else {  
                merged_data <- cbind(merged_data, filtered_expdata)  
            }  
        }    
        cat(paste0("Dim of merged ", assay_name, " data by intersection: "), dim(merged_data), "\n")  
        return(merged_data)  
    }  

    Merge_by_union = function(se_list, assay_name) {  
        merged_data <- NULL  
        expdata1 <- assay(se_list[[1]], assay_name)  
        merged_data <- as.data.frame(matrix(0, nrow = length(all_rows <- union(rownames(expdata1), rownames(expdata1))), ncol = length(all_cols <- colnames(expdata1))))  
        rownames(merged_data) <- all_rows  
        colnames(merged_data) <- all_cols  
        merged_data[rownames(expdata1), colnames(expdata1)] <- expdata1  
        for (se in se_list[-1]) {   
            expdata <- assay(se, assay_name)  
            all_rows <- union(rownames(merged_data), rownames(expdata))  
            all_cols <- union(colnames(merged_data), colnames(expdata))  
            new_merged_data <- as.data.frame(matrix(0, nrow = length(all_rows), ncol = length(all_cols)))  
            rownames(new_merged_data) <- all_rows  
            colnames(new_merged_data) <- all_cols   
            new_merged_data[rownames(merged_data), colnames(merged_data)] <- merged_data   
            new_merged_data[rownames(expdata), colnames(expdata)] <- expdata  
            merged_data <- new_merged_data  
        }  
        cat(paste0("Dim of merged ", assay_name, " data by union: "), dim(merged_data), "\n")  
        return(merged_data)  
    }  

    merged_data_list = list()  
    for (assay_name in all_assay_names) {  
        if (merge_type == "intersection") {  
            merged_data_list[[assay_name]] <- Merge_by_intersection(se_list, assay_name)  
        } else if (merge_type == "union") {  
            merged_data_list[[assay_name]] <- Merge_by_union(se_list, assay_name)  
        } else {  
            stop("Error: merge_type must be either 'intersection' or 'union'.")  
        }  
    }   

    cat("Sample info make")  
    merged_sample_info <- NULL  
    for (i in seq_along(se_list)) {  
        sample_info <- colData(se_list[[i]])  
        current_columns <- names(sample_info)  
        if (is.null(merged_sample_info)) {  
            merged_sample_info <- sample_info  
        } else {  
            all_columns <- union(names(merged_sample_info), current_columns)   
            sample_info[, setdiff(all_columns, current_columns)] <- NA  
            sample_info_reordered <- sample_info[, colnames(merged_sample_info), drop = FALSE]  
            merged_sample_info <- rbind(merged_sample_info, sample_info_reordered)  
        }  
    }  
    cat("Dim of sample info:", dim(merged_sample_info), "\n")  
    cat("Feature info make")  
    merged_feature_info <- data.frame(feature = rownames(merged_data_list[[1]]))  
    rownames(merged_feature_info) <- merged_feature_info$feature  
    cat("Dim of feature info:", dim(merged_feature_info), "\n")  
    cat("SE make")  
    se <- SummarizedExperiment(  
        assays = lapply(merged_data_list, as.matrix),  
        rowData = merged_feature_info,  
        colData = merged_sample_info  
    )  
    cat("SE finished\n")  
    return(se)  
}