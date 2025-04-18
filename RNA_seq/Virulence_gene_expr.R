# 加载所需包
library(readxl)
library(dplyr)
library(readr)
library(openxlsx)

# 读取表格
virulence_genes <- read_excel("D:/NetDrive/OneDrive - wlulab/Work/复旦/课题/MRSA biofilm vs HNP1/20241114 - RNAseq - 29213 VS HNP1/gene_list_virulence.xlsx")
cds_info <- read_csv("D:/NetDrive/OneDrive - wlulab/Work/复旦/课题/MRSA biofilm vs HNP1/20241114 - RNAseq - 29213 VS HNP1/GCF_022870465.1/ncbi_dataset/GCF_022870465.1/CDS_extracted.csv")
rpkm_exp_info <- read_excel("D:/NetDrive/OneDrive - wlulab/Work/复旦/课题/MRSA biofilm vs HNP1/20241114 - RNAseq - 29213 VS HNP1/H_vs_S.rpkm.all.exp.xlsx")

# 读取两个表格
rpkm_data <- read_excel("H_vs_S.rpkm.all.exp.xlsx")
cds_data <- read.csv("./GCF_022870465.1/ncbi_dataset/GCF_022870465.1/CDS_extracted.csv", stringsAsFactors = FALSE)

# 根据 old_locus_tag 和 GeneID 合并表格
merged_data <- rpkm_data %>%
  left_join(select(cds_data, old_locus_tag, gene), by = c("GeneID" = "old_locus_tag"))

# 将结果写入新的 Excel 文件
install.packages("writexl")
library(writexl)
write_xlsx(merged_data, "H_vs_S_with_gene_virulence.xlsx")

# 根据 Gene list 提取 merged_data 中的匹配信息
virulence_gene_data <- merged_data %>%
  filter(gene %in% virulence_genes$Gene)

write_xlsx(virulence_gene_data, "virulence_gene_data.xlsx")


# 加载必要的 R 包
library(readxl)
library(pheatmap)
library(RColorBrewer)

# 读取提取的数据
extracted_data <- virulence_gene_data
extracted_data <- read_excel("./virulence_gene_data.xlsx")

# 选择热图所需的数据列
heatmap_data <- extracted_data[, c("S29-1", "S29-2", "S29-3", "H29-1", "H29-2", "H29-3")]

# 将 gene 列作为行名
rownames(heatmap_data) <- extracted_data$gene

# 对数据进行 log2 转换（避免 log2(0)，加 1 处理）
heatmap_data <- log2(heatmap_data + 1)

# 创建列注释（分组信息）
annotation_col <- data.frame(
  Group = factor(c(rep("Control", 3), rep("HNP1", 3)), levels = c("Control", "HNP1"))
)
rownames(annotation_col) <- colnames(heatmap_data)

# 创建行注释（调整 |logFC| 为三个区间，不考虑 NA）
annotation_row <- data.frame(
  logFC.Abs = factor(
    ifelse(abs(extracted_data$logFC) <= 1, "<=1",
           ifelse(abs(extracted_data$logFC) <= 2, ">1 & <=2", ">2")),
    levels = c("<=1", ">1 & <=2", ">2")
  ),
  FDR = factor(
    ifelse(extracted_data$FDR < 0.05, "<0.05", ">=0.05"),
    levels = c("<0.05", ">=0.05")
  ),
  row.names = extracted_data$gene
)

# 设置颜色方案
color_palette <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
annotation_colors <- list(
  Group = c(Control = "#ea736a", HNP1 = "#1cb4b7"),
  logFC.Abs = c("<=1" = "lightgrey", ">1 & <=2" = "tan1", ">2" = "tomato"),
  FDR = c("<0.05" = "tomato", ">=0.05" = "lightgrey")
)

# 绘制热图（不标准化）
pheatmap(
  heatmap_data,
  scale = "none",  # 不进行标准化，保留原始 log2 值
  cluster_rows = TRUE,  # 行聚类
  cluster_cols = TRUE,  # 列聚类
  annotation_col = annotation_col,  # 列分组注释
  annotation_row = annotation_row,  # 行注释
  annotation_colors = annotation_colors,  # 注释颜色
  color = color_palette,  # 热图颜色
  show_rownames = TRUE,  # 显示行名
  #fontsize_row = 10,  # 行文字大小
  #cellwidth = 20,  # 单元格宽度
  #cellheight = 12,  # 单元格高度
  #main = "Heatmap of Gene Expression (Control vs HNP1, log2 scale)",  # 标题
  labels_row = extracted_data$gene,  # 使用 gene 作为行标签
  annotation_names_row = TRUE  # 显示行注释的名称
)

#------横向作图-----------------------
# 转置 heatmap_data，使列变为行，行变为列
heatmap_data <- t(heatmap_data)  # 使用 t() 函数转置

# 创建新的列注释（原来是行注释，现在用于列）
annotation_col <- data.frame(
  logFC.Abs = factor(
    ifelse(abs(extracted_data$logFC) <= 1, "<=1",
           ifelse(abs(extracted_data$logFC) <= 2, ">1 & <=2", ">2")),
    levels = c("<=1", ">1 & <=2", ">2")
  ),
  FDR = factor(
    ifelse(extracted_data$FDR < 0.05, "<0.05", ">=0.05"),
    levels = c("<0.05", ">=0.05")
  ),
  row.names = extracted_data$gene  # 列名现在是基因名
)

# 创建新的行注释（原来是列注释，现在用于行）
annotation_row <- data.frame(
  Group = factor(c(rep("Control", 3), rep("HNP1", 3)), levels = c("Control", "HNP1")),
  row.names = c("S29-1", "S29-2", "S29-3", "H29-1", "H29-2", "H29-3")  # 行名现在是样本名
)

# 设置颜色方案
color_palette <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
annotation_colors <- list(
  Group = c(Control = "#ea736a", HNP1 = "#1cb4b7"),
  logFC.Abs = c("<=1" = "lightgrey", ">1 & <=2" = "tan1", ">2" = "tomato"),
  FDR = c("<0.05" = "tomato", ">=0.05" = "lightgrey")
)

# 绘制热图（不标准化）
pheatmap(
  heatmap_data,
  scale = "none",  # 不进行标准化，保留原始 log2 值
  cluster_rows = TRUE,  # 行聚类（现在是样本）
  cluster_cols = TRUE,  # 列聚类（现在是基因）
  annotation_col = annotation_col,  # 新的列注释（基因的 logFC 和 FDR）
  annotation_row = annotation_row,  # 新的行注释（样本的 Group）
  annotation_colors = annotation_colors,  # 注释颜色
  color = color_palette,  # 热图颜色
  show_rownames = TRUE,  # 显示行名（样本名）
  show_colnames = TRUE,  # 显示列名（基因名）
  #fontsize_row = 10,  # 行文字大小
  #cellwidth = 20,  # 单元格宽度
  #cellheight = 12,  # 单元格高度
  #main = "Heatmap of Gene Expression (Control vs HNP1, log2 scale)",  # 标题
  labels_row = rownames(heatmap_data),  # 使用样本名作为行标签
  labels_col = colnames(heatmap_data),  # 使用基因名作为列标签
  annotation_names_row = TRUE,  # 显示行注释的名称
  annotation_names_col = TRUE   # 显示列注释的名称
)



# 保存热图为 PDF 文件
dev.copy(pdf, "heatmap_biofilm_genes_fixed.pdf", width = 12, height = 8)
dev.off()


