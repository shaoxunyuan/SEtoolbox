# 测试 SE_boxplot 函数的排序功能

# 加载必要的包
library(SummarizedExperiment)
library(tidyverse)
library(rstatix)
library(ggforce)
library(multcompView)

# 加载函数
source('R/SE_boxplot.R')

# 创建一个虚拟的 SummarizedExperiment 对象
data_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(data_matrix) <- paste0("Gene", 1:100)
colnames(data_matrix) <- paste0("Sample", 1:10)
sample_info <- DataFrame(group = rep(c("A", "B"), each = 5))
SE <- SummarizedExperiment(assays = list(TPM = data_matrix), colData = sample_info)

# 调用 SE_boxplot 函数
print("Running SE_boxplot function...")
result <- SE_boxplot(SE, feature_of_interest = c("Gene1", "Gene2", "Gene3"), 
                   group_colname = "group", normalization = "log")

# 查看结果
print("\n=== BH Results (sorted by p-value) ===")
print(result$diff_results$BH)

print("\n=== FDR Results (sorted by p-value) ===")
print(result$diff_results$fdr)

print("\nTest completed successfully!")
