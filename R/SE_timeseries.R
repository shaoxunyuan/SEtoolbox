#' @title SE_timeseries: Time Series Analysis
#' @description This function performs time series analysis on a SummarizedExperiment object. It can identify temporal patterns, trends, and periodicity in the data.
#' @param SE A \code{SummarizedExperiment} object containing time series data.
#' @param assayname A string indicating which assay to use for time series analysis. The default value is \code{"TPM"}.
#' @param time_col A string indicating the column name in colData containing time information. Default is "time".
#' @param features A character vector of feature names to analyze. If NULL, all features are analyzed. Default is NULL.
#' @param method A character string specifying the analysis method. Options include "trend", "seasonality", "decomposition", "forecast". Default is "trend".
#' @param forecast_periods Numeric value indicating the number of periods to forecast (only used if method = "forecast"). Default is 5.
#' @return A list containing:  
#' \item{results}{Time series analysis results for each feature.}  
#' \item{plots}{A list of ggplot objects showing time series visualizations.}
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Perform trend analysis
#' ts_result <- SE_timeseries(SE, assayname = "TPM", time_col = "time", 
#'                             method = "trend")
#' 
#' # View results
#' print(ts_result$results)
#' 
#' # View plots
#' print(ts_result$plots[[1]])
#' @export
SE_timeseries <- function(SE, assayname = "TPM", time_col = "time", features = NULL, 
                          method = "trend", forecast_periods = 5) {
    
    exp_data <- assay(SE, assayname)
    
    if (!time_col %in% colnames(colData(SE))) {
        stop(paste0("Column '", time_col, "' not found in colData"))
    }
    
    time_points <- colData(SE)[[time_col]]
    
    if (is.null(features)) {
        features <- rownames(exp_data)
    }
    
    results <- list()
    plots <- list()
    
    for (feature in features) {
        if (!(feature %in% rownames(exp_data))) {
            warning(paste0("Feature '", feature, "' not found. Skipping."))
            next
        }
        
        ts_data <- exp_data[feature, ]
        
        if (method == "trend") {
            trend_model <- lm(ts_data ~ time_points)
            trend_slope <- coef(trend_model)[2]
            trend_pvalue <- summary(trend_model)$coefficients[2, 4]
            
            results[[feature]] <- list(
                slope = trend_slope,
                p_value = trend_pvalue,
                trend = ifelse(trend_pvalue < 0.05, 
                               ifelse(trend_slope > 0, "increasing", "decreasing"),
                               "no trend")
            )
            
            plot_df <- data.frame(
                time = time_points,
                expression = ts_data,
                fitted = fitted(trend_model)
            )
            
            ts_plot <- ggplot(plot_df, aes(x = time, y = expression)) +
                geom_point(alpha = 0.5) +
                geom_line(aes(y = fitted), color = "red", size = 1) +
                labs(title = paste("Time Series:", feature),
                     x = "Time",
                     y = "Expression") +
                theme_minimal()
            
        } else if (method == "seasonality") {
            ts_obj <- ts(ts_data, frequency = length(unique(time_points)))
            decomp <- decompose(ts_obj)
            
            results[[feature]] <- list(
                seasonal = decomp$seasonal,
                trend = decomp$trend,
                random = decomp$random
            )
            
            plot_df <- data.frame(
                time = time_points,
                expression = ts_data
            )
            
            ts_plot <- ggplot(plot_df, aes(x = time, y = expression)) +
                geom_line(alpha = 0.5) +
                geom_smooth(method = "loess", se = FALSE, color = "blue") +
                labs(title = paste("Time Series with Seasonality:", feature),
                     x = "Time",
                     y = "Expression") +
                theme_minimal()
            
        } else if (method == "decomposition") {
            ts_obj <- ts(ts_data, frequency = length(unique(time_points)))
            decomp <- decompose(ts_obj)
            
            results[[feature]] <- list(
                observed = decomp$x,
                trend = decomp$trend,
                seasonal = decomp$seasonal,
                random = decomp$random
            )
            
            decomp_df <- data.frame(
                time = time_points,
                observed = as.numeric(decomp$x),
                trend = as.numeric(decomp$trend),
                seasonal = as.numeric(decomp$seasonal)
            )
            
            decomp_long <- reshape2(decomp_df, idvar = "time", timevar = "type", 
                                 direction = "long")
            
            ts_plot <- ggplot(decomp_long, aes(x = time, y = value, color = type)) +
                geom_line() +
                facet_wrap(~ type, scales = "free_y", ncol = 1) +
                labs(title = paste("Time Series Decomposition:", feature),
                     x = "Time",
                     y = "Value",
                     color = "Component") +
                theme_minimal()
            
        } else if (method == "forecast") {
            ts_obj <- ts(ts_data, frequency = length(unique(time_points)))
            fit <- HoltWinters(ts_obj)
            forecast_result <- predict(fit, n.ahead = forecast_periods)
            
            results[[feature]] <- list(
                fitted = fitted(fit),
                forecast = forecast_result,
                model = fit
            )
            
            plot_df <- data.frame(
                time = c(time_points, max(time_points) + 1:forecast_periods),
                expression = c(as.numeric(ts_obj), as.numeric(forecast_result)),
                type = c(rep("observed", length(ts_obj)), 
                          rep("forecast", forecast_periods))
            )
            
            ts_plot <- ggplot(plot_df, aes(x = time, y = expression, color = type)) +
                geom_line() +
                labs(title = paste("Time Series Forecast:", feature),
                     x = "Time",
                     y = "Expression",
                     color = "Type") +
                theme_minimal() +
                theme(legend.position = "bottom")
            
        } else {
            stop("Unknown method. Available options: trend, seasonality, decomposition, forecast")
        }
        
        plots[[feature]] <- ts_plot
    }
    
    cat("Time series analysis completed\n")
    cat("Method:", method, "\n")
    cat("Number of features analyzed:", length(results), "\n")
    
    return(list(
        results = results,
        plots = plots
    ))
}
