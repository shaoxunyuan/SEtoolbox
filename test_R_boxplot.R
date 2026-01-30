# 测试 R_boxplot.R
# 在包根目录运行: Rscript test_R_boxplot.R

# 加载依赖（包内已 Import）
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("需要安装 ggplot2")
if (!requireNamespace("tidyr", quietly = TRUE)) stop("需要安装 tidyr")
if (!requireNamespace("dplyr", quietly = TRUE)) stop("需要安装 dplyr")
if (!requireNamespace("multcompView", quietly = TRUE)) stop("需要安装 multcompView")

# 加载包内函数（相当于 load_all）
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(multcompView)
})
source("R/R_boxplot.R")

message("=== 测试 R_boxplot ===")
set.seed(1)
dat <- data.frame(
  gene1 = c(rnorm(20), rnorm(20, 1), rnorm(20, -0.5)),
  gene2 = c(rnorm(20, 2), rnorm(20), rnorm(20, 0.5)),
  class = rep(c("A", "B", "C"), each = 20)
)

# 1) 基本调用（仅 data，无标题、无 X 轴标签）
p <- R_boxplot(dat)
message("1) R_boxplot(dat) 成功")
stopifnot(inherits(p, "ggplot"))
# 检查无标题、X 轴标签为空
stopifnot(is.null(p$labels$title))
stopifnot(p$labels$x == "")

# 2) 指定 feature_cols
p2 <- R_boxplot(dat, feature_cols = "gene1")
message("2) 指定 feature_cols = 'gene1' 成功")
stopifnot(inherits(p2, "ggplot"))

# 3) 自定义 class_col 名称
dat2 <- dat
names(dat2)[names(dat2) == "class"] <- "group"
p3 <- R_boxplot(dat2, class_col = "group")
message("3) class_col = 'group' 成功")
stopifnot(inherits(p3, "ggplot"))

# 4) 两组的边界情况（Tukey 仍可算）
dat3 <- dat[dat$class %in% c("A", "B"), ]
p4 <- R_boxplot(dat3)
message("4) 两组数据 成功")
stopifnot(inherits(p4, "ggplot"))

message("")
message("全部测试通过。")

# 可选：保存图到文件便于目视检查
if (requireNamespace("ggplot2", quietly = TRUE)) {
  tryCatch({
    ggplot2::ggsave("test_R_boxplot_output.png", p, width = 6, height = 4, dpi = 100)
    message("图已保存为 test_R_boxplot_output.png")
  }, error = function(e) message("保存图片跳过: ", conditionMessage(e)))
}
