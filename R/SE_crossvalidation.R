#' @title SE_crossvalidation: Cross-Validation for Machine Learning Models
#' @description This function performs k-fold cross-validation to evaluate the performance of machine learning models on a SummarizedExperiment object.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to use for cross-validation. The default value is \code{"log2"}.
#' @param group_colname A string representing the column name in \code{colData} that contains group information. Default is "group".
#' @param nfeatures Numeric value indicating the number of top variable features to use. Default is 100.
#' @param model_type A character string specifying the model type. Options include "randomforest", "svm", "logistic". Default is "randomforest".
#' @param k_folds Numeric value indicating the number of folds for cross-validation. Default is 5.
#' @param scale_data Logical value indicating whether to scale the data before cross-validation. Default is TRUE.
#' @return A list containing:  
#' \item{accuracies}{A vector of accuracy scores for each fold.}  
#' \item{mean_accuracy}{Mean accuracy across all folds.}  
#' \item{sd_accuracy}{Standard deviation of accuracy across all folds.}  
#' \item{plot}{A \code{ggplot} object showing accuracy across folds.}
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform 5-fold cross-validation with random forest
#' cv_result <- SE_crossvalidation(SE, assayname = "log2", group_colname = "group", 
#'                                 model_type = "randomforest", k_folds = 5)
#' 
#' # View mean accuracy
#' print(cv_result$mean_accuracy)
#' 
#' # View accuracy plot
#' print(cv_result$plot)
#' @export
SE_crossvalidation <- function(SE, assayname = "log2", group_colname = "group", 
                               nfeatures = 100, model_type = "randomforest", 
                               k_folds = 5, scale_data = TRUE) {
    
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
    folds <- sample(rep(1:k_folds, length.out = nrow(exp_data_subset)))
    
    accuracies <- numeric(k_folds)
    
    for (i in 1:k_folds) {
        train_data <- exp_data_subset[folds != i, ]
        test_data <- exp_data_subset[folds == i, ]
        
        if (model_type == "randomforest") {
            model <- randomForest::randomForest(group ~ ., data = train_data, 
                                                ntree = 500, importance = FALSE)
        } else if (model_type == "svm") {
            model <- e1071::svm(group ~ ., data = train_data, kernel = "radial")
        } else if (model_type == "logistic") {
            model <- glm(group ~ ., data = train_data, family = binomial)
        } else {
            stop("Unknown model type. Available options: randomforest, svm, logistic")
        }
        
        if (model_type == "logistic") {
            predictions <- predict(model, test_data, type = "response")
            predictions <- ifelse(predictions > 0.5, levels(test_data$group)[2], 
                                levels(test_data$group)[1])
        } else {
            predictions <- predict(model, test_data)
        }
        
        confusion_matrix <- table(Actual = test_data$group, Predicted = predictions)
        accuracies[i] <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
    }
    
    mean_accuracy <- mean(accuracies)
    sd_accuracy <- sd(accuracies)
    
    cv_df <- data.frame(
        Fold = 1:k_folds,
        Accuracy = accuracies
    )
    
    plot_cv <- ggplot(cv_df, aes(x = Fold, y = Accuracy)) +
        geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
        geom_hline(yintercept = mean_accuracy, linetype = "dashed", 
                   color = "red", size = 1) +
        geom_text(aes(label = round(Accuracy, 3)), vjust = -0.5) +
        labs(title = paste0(k_folds, "-Fold Cross-Validation Accuracy (", model_type, ")"),
             x = "Fold",
             y = "Accuracy") +
        ylim(0, 1) +
        theme_minimal()
    
    cat("Cross-validation completed\n")
    cat("Model type:", model_type, "\n")
    cat("Number of folds:", k_folds, "\n")
    cat("Mean accuracy:", round(mean_accuracy, 4), "\n")
    cat("Standard deviation:", round(sd_accuracy, 4), "\n")
    
    return(list(
        accuracies = accuracies,
        mean_accuracy = mean_accuracy,
        sd_accuracy = sd_accuracy,
        plot = plot_cv
    ))
}
