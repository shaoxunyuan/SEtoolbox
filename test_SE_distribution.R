# 测试 SE_distribution 函数

# 加载必要的包
library(SummarizedExperiment)
library(tidyverse)
library(reshape2)
library(plyr)
library(RColorBrewer)
library(multcompView)

# 加载函数
source('R/SE_distribution.R')

# 创建示例数据
set.seed(123)
data_matrix <- matrix(rnorm(1000, mean = 5), nrow = 100, ncol = 10)
# 为不同组设置不同的表达水平
data_matrix[, 6:10] <- data_matrix[, 6:10] + 2  # B组表达更高

rownames(data_matrix) <- paste0("Gene", 1:100)
colnames(data_matrix) <- paste0("Sample", 1:10)
sample_info <- DataFrame(
  BioSample = paste0("Sample", 1:10),
  group = rep(c("A", "B"), each = 5)
)
SE <- SummarizedExperiment(assays = list(TPM = data_matrix), colData = sample_info)

# 测试函数
print("Testing SE_distribution function...")

# 测试有分组的情况
print("\n1. Testing with grouping and significance letters...")
plots <- SE_distribution(SE, group_colname = "group", transform = "log2", add_significance = TRUE)

print("Gene density plot created:")
print(plots$gene_density)

print("\nGene boxplot with significance letters created:")
print(plots$gene_box_scatter)

print("\nSample density plot created:")
print(plots$sample_density)

print("\nSample boxplot with significance letters created:")
print(plots$sample_box_scatter)

# 测试无分组的情况
print("\n2. Testing without grouping...")
plots_no_group <- SE_distribution(SE, transform = "raw", add_significance = FALSE)
print("Plots without grouping created successfully")

print("\nAll tests completed successfully!")
