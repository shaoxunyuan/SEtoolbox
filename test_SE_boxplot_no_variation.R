# 测试 SE_boxplot 函数处理无变异基因的情况

# 加载必要的包
library(SummarizedExperiment)
library(tidyverse)
library(rstatix)
library(ggforce)
library(multcompView)

# 加载函数
source('R/SE_boxplot.R')

# 创建一个虚拟的 SummarizedExperiment 对象
# 包含一个有变异的基因和一个无变异的基因
data_matrix <- matrix(c(
  # Gene1: 有变异
  rnorm(5), rnorm(5),
  # Gene2: 无变异（全为 1）
  rep(1, 10)
), nrow = 2, ncol = 10, byrow = TRUE)
rownames(data_matrix) <- c("Gene1", "Gene2")
colnames(data_matrix) <- paste0("Sample", 1:10)
sample_info <- DataFrame(group = rep(c("A", "B"), each = 5))
SE <- SummarizedExperiment(assays = list(TPM = data_matrix), colData = sample_info)

# 调用 SE_boxplot 函数
print("Running SE_boxplot function with no variation gene...")
tryCatch({
  result <- SE_boxplot(SE, feature_of_interest = c("Gene1", "Gene2"), 
                     group_colname = "group", normalization = "log")
  
  # 查看结果
  print("\n=== Plot ===")
  print(result$plot)
  
  print("\n=== Differential Expression Results ===")
  print(result$diff_results)
  
  print("\n=== Expression Data ===")
  print(head(result$expdata))
  
  print("\n=== Sample Data ===")
  print(result$sampledata)
  
  print("\nTest completed successfully! The function handled no variation gene correctly.")
}, error = function(e) {
  print(paste("Error occurred:", e$message))
  print("Test failed! The function did not handle no variation gene correctly.")
})
