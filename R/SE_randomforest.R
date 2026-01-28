#' @title SE_randomforest: Random Forest Classification
#' @description This function performs random forest classification on a SummarizedExperiment object. It can be used for classification tasks such as predicting sample groups based on gene expression data.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for classification. The default value is \code{"log2"}.
#' @param group_colname A string representing the column name in \code{colData} that contains group information. Default is "group".
#' @param nfeatures Numeric value indicating the number of top variable features to use. Default is 100.
#' @param ntree Numeric value indicating the number of trees in the random forest. Default is 500.
#' @param mtry Numeric value indicating the number of variables randomly sampled as candidates at each split. Default is NULL (will use sqrt of number of features).
#' @param train_ratio Numeric value indicating the proportion of data to use for training. Default is 0.7.
#' @param scale_data Logical value indicating whether to scale the data before classification. Default is TRUE.
#' @return A list containing:  
#' \item{model}{The trained random forest model.}  
#' \item{predictions}{Predictions on the test set.}  
#' \item{confusion_matrix}{Confusion matrix of predictions.}  
#' \item{accuracy}{Classification accuracy.}  
#' \item{importance}{Feature importance scores.}  
#' \item{plot_importance}{A \code{ggplot} object showing feature importance.}
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform random forest classification
#' rf_result <- SE_randomforest(SE, assayname = "log2", group_colname = "group")
#' 
#' # View accuracy
#' print(rf_result$accuracy)
#' 
#' # View feature importance plot
#' print(rf_result$plot_importance)
#' @export
SE_randomforest <- function(SE, assayname = "log2", group_colname = "group", 
                            nfeatures = 100, ntree = 500, mtry = NULL, 
                            train_ratio = 0.7, scale_data = TRUE) {
    
    exp_data <- assay(SE, assayname)
    groups <- colData(SE)[[group_colname]]
    
    variances <- apply(exp_data, 1, var, na.rm = TRUE)
    top_features <- names(sort(variances, decreasing = TRUE))[1:min(nfeatures, length(variances))]
    exp_data_subset <- exp_data[top_features, ]
    
    if (scale_data) {
        exp_data_subset <- t(scale(t(exp_data_subset)))
    }
    
    exp_data_subset <- t(exp_data_subset)
    exp_data_subset <- as.data.frame(exp_data_subset)
    exp_data_subset$group <- factor(groups)
    
    set.seed(123)
    train_indices <- sample(1:nrow(exp_data_subset), size = floor(train_ratio * nrow(exp_data_subset)))
    train_data <- exp_data_subset[train_indices, ]
    test_data <- exp_data_subset[-train_indices, ]
    
    if (is.null(mtry)) {
        mtry <- floor(sqrt(nfeatures))
    }
    
    rf_model <- randomForest::randomForest(group ~ ., data = train_data, 
                                          ntree = ntree, mtry = mtry, 
                                          importance = TRUE)
    
    predictions <- predict(rf_model, test_data)
    confusion_matrix <- table(Actual = test_data$group, Predicted = predictions)
    accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
    
    importance <- randomForest::importance(rf_model)
    importance_df <- data.frame(
        feature = rownames(importance),
        MeanDecreaseGini = importance[, "MeanDecreaseGini"]
    )
    importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]
    importance_df$feature <- factor(importance_df$feature, 
                                    levels = importance_df$feature)
    
    plot_importance <- ggplot(importance_df[1:min(20, nrow(importance_df)), ], 
                            aes(x = feature, y = MeanDecreaseGini)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_flip() +
        labs(title = "Random Forest Feature Importance",
             x = "Feature",
             y = "Mean Decrease in Gini") +
        theme_minimal()
    
    cat("Random forest classification completed\n")
    cat("Number of features used:", nfeatures, "\n")
    cat("Number of trees:", ntree, "\n")
    cat("Training samples:", nrow(train_data), "\n")
    cat("Test samples:", nrow(test_data), "\n")
    cat("Accuracy:", round(accuracy, 4), "\n")
    
    return(list(
        model = rf_model,
        predictions = predictions,
        confusion_matrix = confusion_matrix,
        accuracy = accuracy,
        importance = importance,
        plot_importance = plot_importance
    ))
}
