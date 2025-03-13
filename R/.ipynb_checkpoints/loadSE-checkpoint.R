library(SummarizedExperiment)
loadSE <- function() {  
  # 设文件路径  
  file_path <- system.file("extdata", "example.SE", package = "your_package_name")  
  
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
