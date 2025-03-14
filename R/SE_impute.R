#' 填补 SummarizedExperiment 或 phyloseq 对象中的缺失值  
#'  
#' 本函数利用指定的方法对给定的 SummarizedExperiment 对象中的缺失值 (NA) 进行填补。  
#' 可以使用多种插补技术来处理缺失值，以确保后续分析的稳健性。  
#'  
#' @param object 一个 SummarizedExperiment 对象，包含需插补的数据。  
#' @param assayname SummarizedExperiment assay名字，选择待插补数据类型。 
#' @param group 一个字符字符串，指定样本数据中的分组变量。  
#' @param ZerosAsNA 一个逻辑值，指示是否将零视为 NA。默认值为 FALSE。  
#' @param RemoveNA 一个逻辑值，指示是否基于 cutoff 删除高 NA 百分比的样本。默认值为 TRUE。  
#' @param cutoff 一个数值，表示 NA 样本的百分比截断值。默认值为 20。  
#' @param method 一个字符字符串，指定使用的插补方法。可选项包括 "none"、"LOD"、"half_min"、   
#'   "median"、"mean"、"min"、"knn"、"rf"、"global_mean"、"svd" 和 "QRILC"。默认值为 "none"。  
#' @param LOD 一个数值，表示检出限（用于 LOD 插补方法）。默认值为 NULL。  
#' @param knum 一个整数值，表示 KNN 插补方法中邻居的数量。默认值为 10。  
#'  
#' @return SummarizedExperiment 对象，包含填补后的值。  
#'  
#' @examples  
#' # 使用示例  
#' # 假设 `se` 是一个 SummarizedExperiment 对象  
#' se_imputed <- SE_impute(se, assayname = "TPM", group = "my_group_column", method = "median", ZerosAsNA = TRUE)  
#'  
#' @export  
SE_impute <- function(object, assayname = "TPM", group, ZerosAsNA = FALSE, RemoveNA = TRUE,  
                      cutoff = 20, method = c("none", "LOD", "half_min", "median", "mean", "min", "knn", "rf", "global_mean", "svd", "QRILC"),  
                      LOD = NULL, knum = 10) {  
    
    # 匹配方法参数  
    method <- match.arg(method, c("none", "LOD", "half_min", "median", "mean", "min", "knn", "rf", "global_mean", "svd", "QRILC"))  
    
    #if (base::missing(method)) {  
    #    message("method 参数为空！将使用 KNN。")  
    #}  

    # 对于 SummarizedExperiment 对象的处理  
    sam_tab <- SummarizedExperiment::colData(object) 
    sam_tab[] <- lapply(sam_tab, function(x){if(inherits(x,"integer64")) {return(as.numeric(x))}  
        return(x)  
    }) 
    sam_tab = as.data.frame(sam_tab) %>% tibble::rownames_to_column("TempRowNames") 
    prf_tab <- SummarizedExperiment::assay(object,assayname) %>% as.data.frame() %>% t() 

    # 找到样本分组的索引  
    group_index <- which(colnames(sam_tab) == group)  
    samples_groups <- sam_tab[, group_index]  
    to_imp_data <- prf_tab %>% as.matrix()  

    # 将零视为 NA（如果设置为 TRUE）  
    if (ZerosAsNA) {  
        to_imp_data[to_imp_data == 0] <- NA  
        to_imp_data <- data.frame(cbind(Group = samples_groups, to_imp_data))  
        colnames(to_imp_data)[2:ncol(to_imp_data)] <- colnames(prf_tab)  
    } else {  
        to_imp_data <- data.frame(cbind(Group = samples_groups, to_imp_data))  
        colnames(to_imp_data)[2:ncol(to_imp_data)] <- colnames(prf_tab)  
    }  

    # 计算数据中 NA 的比例  
    percent_na <- sum(is.na(to_imp_data))  
    if (percent_na == 0) {  
        print("未检测到缺失值")  
        if (method != "none") {  
            method <- "none"  
        }  
    }  

    # 根据需要删除高缺失值样本  
    if (isTRUE(RemoveNA)) {  
        count_NA <- stats::aggregate(. ~ Group, data = to_imp_data,   
            function(x) {  
                100 * (sum(is.na(x)) / (sum(is.na(x)) + sum(!is.na(x))))  
            }, na.action = NULL)  
        count_NA <- count_NA %>% dplyr::select(-Group)  
        correct_names <- names(count_NA)  
        supress <- unlist(as.data.frame(lapply(count_NA, function(x) all(x > cutoff))))  
        names(supress) <- correct_names  
        correct_names <- names(supress[supress == "FALSE"])  
        depurdata <- to_imp_data[, 2:ncol(to_imp_data)][!supress]  
        samplename <- rownames(depurdata)  
        depurdata <- sapply(depurdata, function(x) as.numeric(as.character(x)))  
        rownames(depurdata) <- samplename  
    } else {  
        depurdata <- to_imp_data[, 2:ncol(to_imp_data)]  
        samplename <- rownames(depurdata)  
        depurdata <- sapply(depurdata, function(x) as.numeric(as.character(x)))  
        correct_names <- colnames(prf_tab)  
        rownames(depurdata) <- samplename  
    }  

    # 根据选择的插补方法填补缺失值  
    if (method == "none") {  
        depurdata[is.na(depurdata)] <- 0  
    } else if (method == "LOD") {  
        if (is.null(LOD)) {  
            message("未提供 LOD，将最小值的十分之一作为 LOD")  
            depurdata_withoutNA <- depurdata[!is.na(depurdata)]  
            LOD <- min(depurdata_withoutNA[depurdata_withoutNA != 0]) / 10  
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
        depurdata <- .GlobalMean(object = t(depurdata)) %>% t()  
    } else if (method == "svd") {  
        depurdata <- .SVD_wrapper(depurdata)  
    } else if (method == "QRILC") {  
        fit <- log(t(depurdata)) %>% imputeLCMD::impute.QRILC()  
        depurdata <- t(fit[[1]])  
    }  

    colnames(depurdata) <- correct_names  
    rownames(depurdata) <- rownames(prf_tab)  

    # 最终返回结果   
        if (ncol(depurdata) != ncol(prf_tab)) {  
            rdata <- SummarizedExperiment::rowData(object)  
            cdata <- SummarizedExperiment::colData(object)  
            if (length(object@metadata) == 0) {  
                mdata <- NULL  
            } else {  
                mdata <- object@metadata  
            }  
            rdata = rdata[rownames(rdata) %in% colnames(depurdata),]
            res <- SummarizedExperiment(assays =  t(depurdata), rowData = rdata, colData = cdata, metadata = mdata)
        } else {  
            res <- object  
            SummarizedExperiment::assay(res) <- t(depurdata)  
        }    

    return(res)  
}