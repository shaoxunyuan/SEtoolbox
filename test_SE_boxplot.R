# 运行 SE_boxplot 示例
# 在包根目录运行: Rscript test_SE_boxplot.R

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(S4Vectors)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggforce)
  library(multcompView)
})
source("R/SE_boxplot.R")

message("=== SE_boxplot 示例 ===")

# 构造示例 SE：3 组（A/B/C），每个基因设计明显组间差异，便于显示 a/b/c
set.seed(42)
n_A <- 8
n_B <- 8
n_C <- 8
n <- n_A + n_B + n_C
n_genes <- 10

# 按列填：前 8 列 A，中 8 列 B，后 8 列 C
data_matrix <- matrix(NA, nrow = n_genes, ncol = n)
sd_val <- 0.25
# Gene1: A 低、B 高、C 中 → 预期 a / c / b 或类似
data_matrix[1, ] <- c(rnorm(n_A, 1, sd_val), rnorm(n_B, 4, sd_val), rnorm(n_C, 2.5, sd_val))
# Gene5: A 高、B 低、C 中 → 预期 c / a / b 或类似
data_matrix[5, ] <- c(rnorm(n_A, 5, sd_val), rnorm(n_B, 1.5, sd_val), rnorm(n_C, 3, sd_val))
# Gene10: A≈B 低、C 高 → 预期 a / a / b 或 ab / ab / c
data_matrix[10, ] <- c(rnorm(n_A, 1, sd_val), rnorm(n_B, 1.2, sd_val), rnorm(n_C, 4, sd_val))
# 其余基因也设一点差异，避免全一样
for (i in setdiff(seq_len(n_genes), c(1, 5, 10))) {
  m <- runif(3, 1, 4)
  data_matrix[i, ] <- c(rnorm(n_A, m[1], sd_val), rnorm(n_B, m[2], sd_val), rnorm(n_C, m[3], sd_val))
}
data_matrix <- pmax(data_matrix, 0.1)

rownames(data_matrix) <- paste0("Gene", 1:n_genes)
colnames(data_matrix) <- paste0("Sample", 1:n)

sample_info <- DataFrame(
  group = rep(c("A", "B", "C"), times = c(n_A, n_B, n_C))
)
SE <- SummarizedExperiment(
  assays = list(TPM = data_matrix),
  colData = sample_info
)

# 选 3 个有差异的基因，log 标准化后组间差异更清晰
plot <- SE_boxplot(
  SE,
  feature_of_interest = c("Gene1", "Gene5", "Gene10"),
  assayname = "TPM",
  group_colname = "group",
  normalization = "log"
)

message("绘图完成。")
print(plot)

# 保存图
ggplot2::ggsave("test_SE_boxplot_output.png", plot, width = 8, height = 4, dpi = 120)
message("图已保存为 test_SE_boxplot_output.png")
