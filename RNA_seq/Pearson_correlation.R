# 加载必要的 R 包
library(readxl)
library(pheatmap)
library(RColorBrewer)

# 读取数据
rpkm_exp_info <- read_excel("D:/NetDrive/OneDrive - wlulab/Work/复旦/课题/MRSA biofilm vs HNP1/20241114 - RNAseq - 29213 VS HNP1/H_vs_S.rpkm.all.exp.xlsx")

# 选择样品数据列（假设包含 S29-1, S29-2, S29-3, H29-1, H29-2, H29-3）
sample_data <- rpkm_exp_info[, c("S29-1", "S29-2", "S29-3", "H29-1", "H29-2", "H29-3")]

# 检查数据是否为数值型，如果不是则转换
sample_data <- as.data.frame(lapply(sample_data, as.numeric))

# 计算 Pearson 相关性矩阵
cor_matrix <- cor(sample_data, method = "pearson", use = "complete.obs")

# 创建列/行注释（分组信息）
annotation <- data.frame(
  Group = factor(c(rep("Control", 3), rep("HNP1", 3)), levels = c("Control", "HNP1")),
  row.names = colnames(sample_data)
)

# 设置颜色方案
color_palette <- colorRampPalette(c("white","#2165ac"))(100)  # 从红到蓝的渐变
annotation_colors <- list(
  Group = c(Control = "#ea736a", HNP1 = "#1cb4b7")
)

# 绘制相关性热图
pheatmap(
  cor_matrix,
  scale = "none",  # 相关性值无需标准化
  cluster_rows = TRUE,  # 行聚类
  cluster_cols = TRUE,  # 列聚类
  annotation_row = annotation,  # 行注释（分组）
  annotation_col = annotation,  # 列注释（分组）
  annotation_colors = annotation_colors,  # 注释颜色
  color = color_palette,  # 热图颜色
  show_rownames = TRUE,  # 显示行名
  show_colnames = TRUE,  # 显示列名
  fontsize = 10,  # 文字大小
  cellwidth = 30,  # 单元格宽度
  cellheight = 30,  # 单元格高度
  main = "Pearson Correlation Heatmap (Control vs HNP1)",  # 标题
  border_color = "lavenderblush4"  # 单元格边框颜色
)

# 保存热图为 PDF 文件
dev.copy(pdf, "correlation_heatmap.pdf", width = 8, height = 8)
dev.off()