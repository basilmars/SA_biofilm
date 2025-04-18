# 加载必要的库
library(dplyr)
library(plotly)
library(ggplot2)

# 读取基因表达矩阵和分组信息
#设置文件路径
dir <- "E://R data/20241122/"

#读取文件
csvfile <- file.path(dir,"genes.rpkm.anno.xls")
cpm <- read.csv(csvfile,row.names = 1,stringsAsFactors = F, header = T, sep = "\t")
head(cpm)

expression_data <- cpm[,1:6]

coldata <- data.frame(
  SampleID = colnames(cpm)[1:6],
  Group = c(rep("ctrl", 3),rep("HNP1",3))
)
head(coldata)

# 移除全零行
non_zero_variance_rows <- apply(expression_data, 1, function(row) var(row) > 0)
expression_data <- expression_data[non_zero_variance_rows, ]

# 确保没有行包含NA值
expression_data <- expression_data[complete.cases(expression_data), ]


# 进行PCA分析
pca_result <- prcomp(t(expression_data), scale. = TRUE)

pca_data <- data.frame(pca_result$x)
head(pca_data)

#计算Principal Components百分比
pcvar <- apply(pca_result$x,2,var)
percentage <- round(pcvar/sum(pcvar)*100,1)
percentage <- paste(colnames(pca_data), " (",paste(as.character(percentage),"%",")", sep = ""))


# 绘制二维PCA图并添加分组椭圆
ggplot(pca_data, aes(x = PC1, y = PC2, color = coldata$Group),
       #xlim = c(-150,150),
       #ylim = c(-150,150)
) +
  geom_point(size = 4) +
  #stat_ellipse(aes(group = coldata$Group), type = "norm", level = 0.95, geom = "polygon", alpha = 0.2) +
  labs(title = "S. aureus",
       x = percentage[1], y = percentage[2]) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "right") + 
  geom_hline(aes(yintercept = 0),colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = 0),colour = "black", linetype = "dashed") +
  ylim(-75,75) +
  xlim(-75,75) +
  ggforce::geom_mark_ellipse(aes(fill = coldata$Group, colour = coldata$Group))





