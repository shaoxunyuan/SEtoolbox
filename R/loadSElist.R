#' Load Example SummarizedExperiment Object List
#'  
#' This function loads the example SummarizedExperiment object list from the package.  
#' The example file should be located in the `extdata` directory of the package.  
#'  
#' @name loadSElist  
#' @return A SummarizedExperiment object list containing example data.  
#' @export  
#'   
#' @examples  
#' library(SEtoolbox)  
#' SElist <- loadSElist()  
#' # Now you can analyze the example SE object  

loadSElist <- function() {  
  library(SummarizedExperiment)
  # 设文件路径  
  file_path <- system.file("extdata", "SElist.rds", package = "SEtoolbox")  
  
  # 检查文件是否存在  
  if (file.exists(file_path)) {  
    # 读取 RDS 文件  
    example_SElist <- readRDS(file_path)   
    return(example_SElist)  
  } else {  
    stop("Example file not found.")  
  }  
}
