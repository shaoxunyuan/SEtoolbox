#' @title SE_normalize: Normalize Expression Data in SummarizedExperiment Object
#' @description This function normalizes expression data in a SummarizedExperiment object using various normalization methods including TPM, FPKM, RPKM, log2 transformation, quantile normalization, and library size normalization.
#' @param SE A \code{SummarizedExperiment} object containing gene expression data.
#' @param assayname A string indicating which assay to normalize. The default value is \code{"Counts"}.
#' @param method A character string specifying the normalization method. Options include "TPM", "FPKM", "RPKM", "log2", "quantile", "library_size", "median", "upper_quartile". Default is "log2".
#' @param pseudocount Numeric value to add before log transformation to avoid log(0). Default is 1.
#' @param gene_length A numeric vector of gene lengths in base pairs. Required for TPM/FPKM/RPKM methods.
#' @return A normalized \code{SummarizedExperiment} object with the normalized data stored in a new assay named after the method (e.g., "TPM", "log2").
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Log2 transformation
#' SE_log2 <- SE_normalize(SE, assayname = "Counts", method = "log2")
#' 
#' # Library size normalization
#' SE_libsize <- SE_normalize(SE, assayname = "Counts", method = "library_size")
#' 
#' # TPM normalization (requires gene lengths)
#' gene_lengths <- rnorm(nrow(SE), mean = 2000, sd = 500)
#' SE_TPM <- SE_normalize(SE, assayname = "Counts", method = "TPM", gene_length = gene_lengths)
#' @export
SE_normalize <- function(SE, assayname = "Counts", method = "log2", pseudocount = 1, gene_length = NULL) {
    
    exp_data <- assay(SE, assayname)
    
    if (method == "log2") {
        normalized_data <- log2(exp_data + pseudocount)
        cat("Log2 transformation applied with pseudocount =", pseudocount, "\n")
    } else if (method == "library_size") {
        lib_sizes <- colSums(exp_data, na.rm = TRUE)
        scaling_factors <- lib_sizes / mean(lib_sizes)
        normalized_data <- sweep(exp_data, 2, scaling_factors, "/")
        cat("Library size normalization applied\n")
    } else if (method == "median") {
        medians <- apply(exp_data, 2, median, na.rm = TRUE)
        scaling_factors <- medians / median(medians)
        normalized_data <- sweep(exp_data, 2, scaling_factors, "/")
        cat("Median normalization applied\n")
    } else if (method == "upper_quartile") {
        upper_quartiles <- apply(exp_data, 2, function(x) quantile(x, 0.75, na.rm = TRUE))
        scaling_factors <- upper_quartiles / median(upper_quartiles)
        normalized_data <- sweep(exp_data, 2, scaling_factors, "/")
        cat("Upper quartile normalization applied\n")
    } else if (method == "quantile") {
        normalized_data <- normalizeBetweenArrays(as.matrix(exp_data), method = "quantile")
        cat("Quantile normalization applied\n")
    } else if (method == "TPM") {
        if (is.null(gene_length)) {
            stop("gene_length is required for TPM normalization")
        }
        rpk <- exp_data / gene_length
        scaling_factor <- colSums(rpk, na.rm = TRUE) / 1e6
        normalized_data <- sweep(rpk, 2, scaling_factor, "/")
        cat("TPM normalization applied\n")
    } else if (method == "FPKM") {
        if (is.null(gene_length)) {
            stop("gene_length is required for FPKM normalization")
        }
        lib_sizes <- colSums(exp_data, na.rm = TRUE)
        normalized_data <- (exp_data / gene_length) / (lib_sizes / 1e6)
        cat("FPKM normalization applied\n")
    } else if (method == "RPKM") {
        if (is.null(gene_length)) {
            stop("gene_length is required for RPKM normalization")
        }
        lib_sizes <- colSums(exp_data, na.rm = TRUE)
        normalized_data <- (exp_data / gene_length) / (lib_sizes / 1e6)
        cat("RPKM normalization applied\n")
    } else {
        stop("Unknown normalization method. Available methods: log2, library_size, median, upper_quartile, quantile, TPM, FPKM, RPKM")
    }
    
    assay(SE, method) <- normalized_data
    cat("Normalized data stored in assay:", method, "\n")
    
    return(SE)
}
