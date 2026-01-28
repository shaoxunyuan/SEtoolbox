#' @title SE_SVM: Support Vector Machine Classification
#' @description This function performs Support Vector Machine (SVM) classification on a SummarizedExperiment object. It can be used for classification tasks such as predicting sample groups based on gene expression data.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for classification. The default value is \code{"log2"}.
#' @param group_colname A string representing the column name in \code{colData} that contains group information. Default is "group".
#' @param nfeatures Numeric value indicating the number of top variable features to use. Default is 100.
#' @param kernel A character string specifying the kernel type. Options include "linear", "polynomial", "radial", "sigmoid". Default is "radial".
#' @param cost Numeric value for the cost parameter. Default is 1.
#' @param train_ratio Numeric value indicating the proportion of data to use for training. Default is 0.7.
#' @param scale_data Logical value indicating whether to scale the data before classification. Default is TRUE.
#' @return A list containing:  
#' \item{model}{The trained SVM model.}  
#' \item{predictions}{Predictions on the test set.}  
#' \item{confusion_matrix}{Confusion matrix of predictions.}  
#' \item{accuracy}{Classification accuracy.}
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform SVM classification
#' svm_result <- SE_SVM(SE, assayname = "log2", group_colname = "group")
#' 
#' # View accuracy
#' print(svm_result$accuracy)
#' 
#' # View confusion matrix
#' print(svm_result$confusion_matrix)
#' @export
SE_SVM <- function(SE, assayname = "log2", group_colname = "group", 
                   nfeatures = 100, kernel = "radial", cost = 1, 
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
    
    svm_model <- e1071::svm(group ~ ., data = train_data, kernel = kernel, cost = cost)
    
    predictions <- predict(svm_model, test_data)
    confusion_matrix <- table(Actual = test_data$group, Predicted = predictions)
    accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
    
    cat("SVM classification completed\n")
    cat("Number of features used:", nfeatures, "\n")
    cat("Kernel:", kernel, "\n")
    cat("Cost:", cost, "\n")
    cat("Training samples:", nrow(train_data), "\n")
    cat("Test samples:", nrow(test_data), "\n")
    cat("Accuracy:", round(accuracy, 4), "\n")
    
    return(list(
        model = svm_model,
        predictions = predictions,
        confusion_matrix = confusion_matrix,
        accuracy = accuracy
    ))
}
