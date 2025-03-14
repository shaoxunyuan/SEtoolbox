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
#' example_SE <- loadSE()  
#' # Now you can analyze the example_SE object  
loadSE <- function() {  
  library(SummarizedExperiment)
  # 设文件路径  
  file_path <- system.file("extdata", "example.SE", package = "SEtoolbox")  
  
  # 检查文件是否存在  
  if (file.exists(file_path)) {  
    # 读取 RDS 文件  
    example_SE <- readRDS(file_path)  
    
    # 检查是否为 SummarizedExperiment 对象  
    if (!is(example_SE, "SummarizedExperiment")) {  
      stop("The loaded file is not a SummarizedExperiment object.")  
    }  
    
    # 返回 SummarizedExperiment 对象  
    return(example_SE)  
  } else {  
    stop("Example file not found.")  
  }  
}
