#' @title Calculate AUC based on gene expressions in a SummarizedExperiment object  
#'   
#' @description  
#' This function calculates the Area Under the Curve (AUC) for various combinations   
#' of genes in the provided SummarizedExperiment (SE) object. It allows for weighted   
#' genes and enables users to specify the maximum number of features to combine.  
#'  
#' @param SE A SummarizedExperiment object containing gene expression data   
#'   and associated metadata.  
#' @param assayname A string specifying which assay to use from the   
#'   SummarizedExperiment object (default is "TPM").  
#' @param group_colname A string specifying the name of the column   
#'   in the metadata that defines the grouping of samples (default is "group").  
#' @param feature_of_interest A data.frame containing at least two columns:   
#'   `geneList` (the gene names) and `geneWeight` (the weights for each gene).  
#' @param maxfeaturecount An integer specifying the maximum number of features   
#'   (genes) to consider in the combinations (default is 2).  
#'   
#' @return A data.frame containing the best AUC results for the combinations   
#'   of genes tested. Each row corresponds to a particular combination of genes   
#'   along with its AUC score.  
#'   
 
#' @importFrom pROC roc  
#' @importFrom plyr mapvalues  
#' @export  
#'   
SE_AUCcalc = function(SE, assayname = "TPM", group_colname = "group", feature_of_interest, maxfeaturecount = 2) {  

    # 检查 maxfeaturecount  
    if (!is.numeric(maxfeaturecount) || maxfeaturecount <= 0 || maxfeaturecount %% 1 != 0) {  
        stop("maxfeaturecount must be a positive integer.")  
    }  

    # 数据处理  
    expdata = as.data.frame(t(assay(SE, assayname)))  
    expdata = log2(expdata + 1)  
    sample_info = as.data.frame(colData(SE))  
    
    # 映射组别  
    expdata$group = mapvalues(rownames(expdata), rownames(sample_info), sample_info[, group_colname], warn_missing = FALSE)  
    feature_of_interest = feature_of_interest[feature_of_interest$geneList %in% intersect(feature_of_interest$geneList, names(expdata)), ]  
    expdata = expdata[, c(feature_of_interest$geneList, group_colname)]  

    # 加权基因的计算  
    for(gene in feature_of_interest$geneList) {  
        expdata[, gene] = expdata[, gene] * feature_of_interest[feature_of_interest$geneList == gene, "geneWeight"]  
    }  

    # 获取所有基因，单个基因计算AUC  
    genes = setdiff(names(expdata), group_colname)  
    all_combinations <- lapply(1:length(genes), function(r) combn(genes, r, simplify = FALSE))   
    all_combinations <- unlist(all_combinations, recursive = FALSE)  

    # 计算所有组合的AUC值  
    all_combinations1 = all_combinations[sapply(all_combinations, length) == 1]  
    auc_result.df = data.frame()  
    
    for(index in 1:length(all_combinations1)) {  
        expdatasub = expdata[, c(all_combinations1[[index]], group_colname)]  
        if(ncol(expdatasub) > 2) {  
            expdatasub$score <- rowSums(expdatasub[, -ncol(expdatasub)])   
            expdatasub <- expdatasub[, c("score", group_colname)]   
        }  
        names(expdatasub) = c("score", group_colname)  
        roc_curve <- roc(expdatasub[, group_colname], expdatasub[,"score"], quiet = TRUE)  
        results = data.frame(features = paste(all_combinations1[[index]], collapse = ","), auc = roc_curve$auc)  
        auc_result.df = rbind(auc_result.df, results)  
    }  

    loop_auc = auc_result.df[order(auc_result.df$auc, decreasing = TRUE), ]  
    loop_auc.best = loop_auc[which(loop_auc$auc == max(loop_auc$auc)), ]  
    
	final_results = data.frame() 
    final_results.best = data.frame()  
	
    final_results = rbind(final_results, loop_auc)  
	final_results.best = rbind(final_results.best, loop_auc.best)  

    # 循环计算  
    for(loopnumber in 2:maxfeaturecount) {             
        last_best_genes = final_results$features  
        loop_genes = unique(unlist(strsplit(last_best_genes, split = ",")))  
        loop_combines <- Filter(function(subset) {length(subset) == loopnumber && any(loop_genes %in% subset)}, all_combinations)  
        
        auc_result.df = data.frame()  
        for(index in 1:length(loop_combines)) {  
            expdatasub = expdata[, c(loop_combines[[index]], group_colname)]  
            if(ncol(expdatasub) > 2) {  
                expdatasub$score <- rowSums(expdatasub[, -ncol(expdatasub)])   
                expdatasub <- expdatasub[, c("score", group_colname)]   
            }  
            names(expdatasub) = c("score", group_colname)  
            roc_curve <- roc(expdatasub[, group_colname], expdatasub[,"score"], quiet = TRUE)  
            results = data.frame(features = paste(loop_combines[[index]], collapse = ","), auc = roc_curve$auc)  
            auc_result.df = rbind(auc_result.df, results)  
        }  
        
        if (nrow(auc_result.df) == 0) {  
            break # 跳出循环如果没有 AUC 结果  
        }  
        
        loop_auc = auc_result.df[order(auc_result.df$auc, decreasing = TRUE), ]  
        loop_auc.best = loop_auc[which(loop_auc$auc == max(loop_auc$auc)), ]  
        final_results = rbind(final_results, loop_auc)  
	    final_results.best = rbind(final_results.best, loop_auc.best)  
    }  
  
    return(list(AUC_results = final_results,AUC_results.best = final_results.best))  
}