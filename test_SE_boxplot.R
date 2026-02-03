# 测试 SE_boxplot 函数

# 加载必要的包
library(SEtoolbox)
library(SummarizedExperiment)

# 创建一个虚拟的 SummarizedExperiment 对象
data_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(data_matrix) <- paste0("Gene", 1:100)
colnames(data_matrix) <- paste0("Sample", 1:10)
sample_info <- DataFrame(group = rep(c("A", "B"), each = 5))
SE <- SummarizedExperiment(assays = list(TPM = data_matrix), colData = sample_info)

# 调用 SE_boxplot 函数
print("Running SE_boxplot function...")
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

print("\nTest completed successfully!")
