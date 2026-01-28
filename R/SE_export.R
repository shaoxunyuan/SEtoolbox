#' @title SE_export: Export SummarizedExperiment Object to Various Formats
#' @description This function exports a SummarizedExperiment object to various file formats including CSV, Excel, TSV, and RDS.
#' @param SE A \code{SummarizedExperiment} object to export.
#' @param assayname A string indicating which assay to export. The default value is \code{"TPM"}.
#' @param format A character string specifying the output format. Options include "csv", "excel", "tsv", "rds". Default is "csv".
#' @param output_dir A string indicating the directory to save the exported file. Default is current directory.
#' @param filename A string indicating the base filename (without extension). Default is "SE_export".
#' @return A character string with the path to the exported file.
#' @examples 
#' # Load example SummarizedExperiment object
#' SE <- loadSE()
#' 
#' # Export to CSV
#' csv_path <- SE_export(SE, assayname = "TPM", format = "csv", filename = "my_data")
#' 
#' # Export to Excel
#' excel_path <- SE_export(SE, assayname = "TPM", format = "excel", filename = "my_data")
#' @export
SE_export <- function(SE, assayname = "TPM", format = "csv", 
                      output_dir = ".", filename = "SE_export") {
    
    exp_data <- assay(SE, assayname)
    
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    if (format == "csv") {
        output_path <- file.path(output_dir, paste0(filename, ".csv"))
        write.csv(exp_data, output_path)
        cat("Exported to CSV:", output_path, "\n")
        
    } else if (format == "excel") {
        output_path <- file.path(output_dir, paste0(filename, ".xlsx"))
        write.xlsx(exp_data, output_path)
        cat("Exported to Excel:", output_path, "\n")
        
    } else if (format == "tsv") {
        output_path <- file.path(output_dir, paste0(filename, ".tsv"))
        write.table(exp_data, output_path, sep = "\t", quote = FALSE, row.names = TRUE)
        cat("Exported to TSV:", output_path, "\n")
        
    } else if (format == "rds") {
        output_path <- file.path(output_dir, paste0(filename, ".rds"))
        saveRDS(SE, output_path)
        cat("Exported to RDS:", output_path, "\n")
        
    } else {
        stop("Unknown format. Available options: csv, excel, tsv, rds")
    }
    
    return(output_path)
}
