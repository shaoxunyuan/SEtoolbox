# SEtoolbox

SEtoolbox 是一个用于 SummarizedExperiment 分析的综合工具包

# 安装

您可以直接从 GitHub 安装该包，

```r
# install.packages("devtools")
devtools::install_github("shaoxunyuan/SEtoolbox")
```

要运行我们的 [vignette](
https://shaoxunyuan.github.io/SEtoolbox/
) 中的示例代码，请将 `dependencies` 参数设置为 `TRUE`，

```r
# install.packages("devtools")

devtools::install_github("shaoxunyuan/SEtoolbox", dependencies = TRUE)
```

# 快速开始

## 加载示例 SE 对象

要加载示例 `SummarizedExperiment` 对象，请使用以下命令：  

```r
library(SEtoolbox)  
library(SummarizedExperiment)

# 单个 SE 对象
SE  = loadSE()

# SE 对象列表。包含 3 个 SE 对象的列表。
SElist  = loadSElist()
```

## 使用示例 SE 对象运行函数

现在您可以将 SE 与所有函数一起使用 

### SE_combine

将多个 SummarizedExperiment 对象合并为一个 SummarizedExperiment 对象

```r
# 包含多个 SE 对象的列表。所有 SE 对象应具有相同的组列名（参数：group_col）和健康对照标签（参数：label_healthy）。

SElist = loadSElist()

SE_combine(SElist,merge_type = "intersection")  # 保留 SE 对象之间的交集特征

SE_combine(SElist,merge_type = "union")  # 保留 SE 对象之间的并集特征
```

### SE_filter

基于表达水平、检测率、方差和变异系数过滤特征

```r
library(SEtoolbox)

SE = loadSE()

# 过滤最小表达量为 1 且最小检测率为 0.5 的特征
SE_filtered = SE_filter(SE, assayname = "TPM", min_expr = 1, min_detectratio = 0.5)

# 过滤最小 CV 为 0.1 的特征
SE_filtered = SE_filter(SE, assayname = "TPM", min_cv = 0.1)
```

### SE_subset

按特征或样本子集化 SummarizedExperiment 对象

```r
library(SEtoolbox)

SE = loadSE()

# 按特定特征子集化
SE_subset = SE_subset(SE, features = c("Gene1", "Gene2", "Gene3"))

# 按 colData 中的条件子集化
SE_subset = SE_subset(SE, condition = list(group = "Treatment"))
```

### SE_normalize

使用各种方法（TPM、FPKM、RPKM、log2、quantile 等）标准化表达数据

```r
library(SEtoolbox)

SE = loadSE()

# Log2 转换
SE_log2 = SE_normalize(SE, assayname = "Counts", method = "log2")

# 库大小标准化
SE_libsize = SE_normalize(SE, assayname = "Counts", method = "library_size")
```

### SE_transform

转换表达数据（log2、log10、sqrt 等）

```r
library(SEtoolbox)

SE = loadSE()

# Log10 转换
SE_log10 = SE_transform(SE, assayname = "TPM", method = "log10")

# 平方根转换
SE_sqrt = SE_transform(SE, assayname = "TPM", method = "sqrt")
```

### SE_batchdetect

使用 PCA 可视化和统计检验检测批次效应

```r
library(SEtoolbox)

SE = loadSE()

# 检测批次效应
result = SE_batchdetect(SE, assayname = "TPM", batch_colname = "batch", group_colname = "group")

# 查看 PCA 图
print(result$plot_pca_batch)
print(result$plot_pca_group)
```

### SE_limma

使用 limma 进行差异表达分析

```r
library(SEtoolbox)

SE = loadSE()

# 使用 limma 进行差异表达分析
SE_limma_result = SE_limma(SE, assayname = "log2", group_colname = "group")
```

### SE_edgeR

使用 edgeR 进行差异表达分析

```r
library(SEtoolbox)

SE = loadSE()

# 使用 edgeR 进行差异表达分析
SE_edgeR_result = SE_edgeR(SE, assayname = "Counts", group_colname = "group")
```

### SE_volcano

为差异表达结果创建火山图

```r
library(SEtoolbox)

SE = loadSE()

# 进行差异表达分析
SE_limma_result = SE_limma(SE, assayname = "log2", group_colname = "group")

# 创建火山图
volcano_plot = SE_volcano(SE_limma_result, logFC_col = "logFC", pvalue_col = "adj.P.Val")
print(volcano_plot)
```

### SE_MAplot

为差异表达结果创建 MA 图

```r
library(SEtoolbox)

SE = loadSE()

# 进行差异表达分析
SE_limma_result = SE_limma(SE, assayname = "log2", group_colname = "group")

# 创建 MA 图
ma_plot = SE_MAplot(SE_limma_result, assayname = "log2", logFC_col = "logFC", pvalue_col = "adj.P.Val")
print(ma_plot)
```

### SE_GSEA

基因集富集分析

```r
library(SEtoolbox)

SE = loadSE()

# 定义基因集
gene_sets = list(
    Pathway1 = c("Gene1", "Gene2", "Gene3"),
    Pathway2 = c("Gene4", "Gene5", "Gene6")
)

# 执行 GSEA
gsea_results = SE_GSEA(SE, assayname = "log2", group_colname = "group", gene_sets = gene_sets)
print(gsea_results)
```

### SE_GO

GO 富集分析

```r
library(SEtoolbox)

SE = loadSE()

# 进行差异表达分析
SE_limma_result = SE_limma(SE, assayname = "log2", group_colname = "group")

# 执行 GO 富集分析
go_results = SE_GO(SE_limma_result, ontology = "BP")
print(go_results)
```

### SE_KEGG

KEGG 通路富集分析

```r
library(SEtoolbox)

SE = loadSE()

# 进行差异表达分析
SE_limma_result = SE_limma(SE, assayname = "log2", group_colname = "group")

# 执行 KEGG 富集分析
kegg_results = SE_KEGG(SE_limma_result, organism = "hsa")
print(kegg_results)
```

### SE_enrichplot

可视化富集分析结果

```r
library(SEtoolbox)

SE = loadSE()

# 执行 GO 富集分析
go_results = SE_GO(SE_limma_result, ontology = "BP")

# 创建条形图
bar_plot = SE_enrichplot(go_results, plot_type = "bar", top_n = 15)
print(bar_plot)

# 创建点图
dot_plot = SE_enrichplot(go_results, plot_type = "dot", top_n = 15)
print(dot_plot)
```

### SE_hierarchical

层次聚类分析

```r
library(SEtoolbox)

SE = loadSE()

# 对样本进行层次聚类
SE_clustered = SE_hierarchical(SE, assayname = "log2", cluster_by = "samples", n_clusters = 3)

# 查看聚类分配
print(colData(SE_clustered)$cluster)
```

### SE_kmeans

K-means 聚类分析

```r
library(SEtoolbox)

SE = loadSE()

# 对样本进行 K-means 聚类
SE_clustered = SE_kmeans(SE, assayname = "log2", cluster_by = "samples", n_clusters = 3)

# 查看聚类分配
print(colData(SE_clustered)$cluster)
```

### SE_clusterplot

可视化聚类结果

```r
library(SEtoolbox)

SE = loadSE()

# 对样本进行 K-means 聚类
SE_clustered = SE_kmeans(SE, assayname = "log2", cluster_by = "samples", n_clusters = 3)

# 创建按聚类着色的 PCA 图
pca_plot = SE_clusterplot(SE_clustered, plot_type = "pca", cluster_by = "samples")
print(pca_plot)
```

### SE_silhouette

计算聚类的轮廓系数

```r
library(SEtoolbox)

SE = loadSE()

# 对样本进行 K-means 聚类
SE_clustered = SE_kmeans(SE, assayname = "log2", cluster_by = "samples", n_clusters = 3)

# 计算轮廓宽度
silhouette_result = SE_silhouette(SE_clustered, cluster_by = "samples")

# 查看平均轮廓宽度
print(silhouette_result$mean_silhouette)

# 绘制轮廓宽度
print(silhouette_result$plot)
```

### SE_tSNE

t-SNE 降维和可视化

```r
library(SEtoolbox)

SE = loadSE()

# 执行 t-SNE
tsne_result = SE_tSNE(SE, assayname = "log2", color_by = "group")

# 查看 t-SNE 图
print(tsne_result$plot)
```

### SE_UMAP

UMAP 降维和可视化

```r
library(SEtoolbox)

SE = loadSE()

# 执行 UMAP
umap_result = SE_UMAP(SE, assayname = "log2", color_by = "group")

# 查看 UMAP 图
print(umap_result$plot)
```

### SE_MDS

多维标度分析

```r
library(SEtoolbox)

SE = loadSE()

# 执行 MDS
mds_result = SE_MDS(SE, assayname = "log2", color_by = "group")

# 查看 MDS 图
print(mds_result$plot)
```

### SE_randomforest

随机森林分类

```r
library(SEtoolbox)

SE = loadSE()

# 执行随机森林分类
rf_result = SE_randomforest(SE, assayname = "log2", group_colname = "group")

# 查看准确率
print(rf_result$accuracy)

# 查看特征重要性图
print(rf_result$plot_importance)
```

### SE_SVM

支持向量机分类

```r
library(SEtoolbox)

SE = loadSE()

# 执行 SVM 分类
svm_result = SE_SVM(SE, assayname = "log2", group_colname = "group")

# 查看准确率
print(svm_result$accuracy)
```

### SE_crossvalidation

机器学习模型的交叉验证

```r
library(SEtoolbox)

SE = loadSE()

# 使用随机森林进行 5 折交叉验证
cv_result = SE_crossvalidation(SE, assayname = "log2", group_colname = "group", model_type = "randomforest", k_folds = 5)

# 查看平均准确率
print(cv_result$mean_accuracy)

# 查看准确率图
print(cv_result$plot)
```

### SE_featureselection

机器学习的特征选择

```r
library(SEtoolbox)

SE = loadSE()

# 按方差选择前 50 个特征
SE_selected = SE_featureselection(SE, assayname = "log2", group_colname = "group", method = "variance", nfeatures = 50)

# 按 limma 差异表达选择特征
SE_selected = SE_featureselection(SE, assayname = "log2", group_colname = "group", method = "limma", nfeatures = 50)
```

### SE_WGCNA

加权基因共表达网络分析

```r
library(SEtoolbox)

SE = loadSE()

# 执行 WGCNA
wgcna_result = SE_WGCNA(SE, assayname = "log2", nfeatures = 5000)

# 查看模块分配
print(head(wgcna_result$moduleColors))
```

### SE_correlation

相关性分析

```r
library(SEtoolbox)

SE = loadSE()

# 计算特征-特征相关性
cor_result = SE_correlation(SE, assayname = "log2", method = "features")

# 计算样本-样本相关性
cor_result = SE_correlation(SE, assayname = "log2", method = "samples")

# 计算特征-性状相关性
cor_result = SE_correlation(SE, assayname = "log2", method = "feature_trait", trait_col = "group")
```

### SE_networkplot

从相关矩阵创建网络图

```r
library(SEtoolbox)

SE = loadSE()

# 计算特征-特征相关性
cor_matrix = SE_correlation(SE, assayname = "log2", method = "features")

# 创建网络图
network_plot = SE_networkplot(cor_matrix, threshold = 0.7)
print(network_plot)
```

### SE_QCreport

生成质量控制报告

```r
library(SEtoolbox)

SE = loadSE()

# 生成 QC 报告
qc_report = SE_QCreport(SE, assayname = "TPM")

# 查看 QC 摘要
print(qc_report$summary)

# 查看 QC 图
print(qc_report$plots$boxplot)
print(qc_report$plots$density)
print(qc_report$plots$pca)
```

### SE_sampleQC

样本质量评估

```r
library(SEtoolbox)

SE = loadSE()

# 评估样本质量
SE_qc = SE_sampleQC(SE, assayname = "TPM")

# 查看 QC 指标
print(colData(SE_qc))

# 识别异常样本
outliers = which(colData(SE_qc)$outlier)
print("异常样本：")
print(colnames(SE_qc)[outliers])
```

### SE_featureQC

特征质量评估

```r
library(SEtoolbox)

SE = loadSE()

# 评估特征质量
SE_qc = SE_featureQC(SE, assayname = "TPM")

# 查看 QC 指标
print(rowData(SE_qc))

# 识别低质量特征
low_quality = which(rowData(SE_qc)$low_quality)
print("低质量特征：")
print(rownames(SE_qc)[low_quality])
```

### SE_export

将 SummarizedExperiment 对象导出为各种格式

```r
library(SEtoolbox)

SE = loadSE()

# 导出为 CSV
csv_path = SE_export(SE, assayname = "TPM", format = "csv", filename = "my_data")

# 导出为 Excel
excel_path = SE_export(SE, assayname = "TPM", format = "excel", filename = "my_data")
```

### SE_summary

为 SummarizedExperiment 对象生成摘要统计

```r
library(SEtoolbox)

SE = loadSE()

# 生成摘要统计
summary_result = SE_summary(SE, assayname = "TPM", group_colname = "group")

# 查看总体摘要
print(summary_result$overall)

# 查看样本统计
print(summary_result$sample_stats)

# 查看组统计
print(summary_result$group_stats)
```

### SE_metadata

从 SummarizedExperiment 对象提取和操作元数据

```r
library(SEtoolbox)

SE = loadSE()

# 提取所有元数据
metadata = SE_metadata(SE, metadata_type = "both")

# 查看 colData
print(metadata$coldata)

# 查看 rowData
print(metadata$rowdata)

# 向 colData 添加新列
SE_updated = SE_metadata(SE, add_coldata = list(new_col = c(1, 2, 3)))
```

### SE_timeseries

时间序列分析

```r
library(SEtoolbox)

SE = loadSE()

# 执行趋势分析
ts_result = SE_timeseries(SE, assayname = "TPM", time_col = "time", method = "trend")

# 查看结果
print(ts_result$results)

# 查看图
print(ts_result$plots[[1]])
```

### SE_trend

趋势分析

```r
library(SEtoolbox)

SE = loadSE()

# 分析特征趋势
SE_trend = SE_trend(SE, assayname = "TPM", trend_col = "time", by = "features")

# 查看趋势结果
print(rowData(SE_trend)$trend_direction)
print(rowData(SE_trend)$trend_significant)
```

### SE_merge

从多个 SummarizedExperiment 对象合并 rowData 或 colData

```r
library(SEtoolbox)

SE1 = loadSE()
SE2 = loadSE()

# 从多个 SE 对象合并 colData
SE_merged = SE_merge(list(SE1, SE2), merge_type = "coldata", merge_method = "outer")

# 合并 rowData 和 colData
SE_merged = SE_merge(list(SE1, SE2), merge_type = "both", merge_method = "outer")
```

### SE_rename

重命名 SummarizedExperiment 对象中的特征或样本

```r
library(SEtoolbox)

SE = loadSE()

# 使用映射重命名特征
mapping = c("Gene1" = "NewGene1", "Gene2" = "NewGene2")
SE_renamed = SE_rename(SE, rename_type = "features", mapping = mapping)

# 向样本名称添加前缀
SE_renamed = SE_rename(SE, rename_type = "samples", prefix = "Sample_")
```

### SE_convert

将 SummarizedExperiment 转换为其他数据结构

```r
library(SEtoolbox)

SE = loadSE()

# 转换为数据框
df = SE_convert(SE, assayname = "TPM", target_format = "data.frame")

# 转换为矩阵
mat = SE_convert(SE, assayname = "TPM", target_format = "matrix")

# 转换为列表
lst = SE_convert(SE, assayname = "TPM", target_format = "list")
```

# Vignette

要充分利用我们的包，请参考我们的 [文档](
https://shaoxunyuan.github.io/SEtoolbox/)
（强烈推荐）。 

# Bug 报告

如果您遇到任何问题、有疑问或想提出建议， 
请随时在我们的 

[bug 追踪器](https://github.com/shaoxunyuan/SEtoolbox/issues) 上报告。

# 联系方式

如需更多咨询，请联系我们： 
邮箱：shaoxunyuan@njucm.edu.cn
