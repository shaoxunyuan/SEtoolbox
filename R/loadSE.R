#' Load Example SummarizedExperiment Object
#'
#' This function loads the example SummarizedExperiment object from the package.  
#' The example file should be located in the `extdata` directory of the package.  
#'
#' @name loadSE
#' @return A SummarizedExperiment object containing example data.
#' @export
#'
#' @examples
#' library(SEtoolbox)
#' SE <- loadSE()
#' # Now you can analyze the example SE object
loadSE <- function() {
  # Set file path
  file_path <- system.file("extdata", "SE.rds", package = "SEtoolbox")

  # Check if file exists
  if (file.exists(file_path)) {
    # Read RDS file
    example_se <- readRDS(file_path)

    # Check if it's a SummarizedExperiment object
    if (!is(example_se, "SummarizedExperiment")) {
      stop("The loaded file is not a SummarizedExperiment object.")
    }

    # Return SummarizedExperiment object
    example_se
  } else {
    stop("Example file not found.")
  }
}
