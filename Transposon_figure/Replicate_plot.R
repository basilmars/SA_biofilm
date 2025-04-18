rm(list = ls(all.names = TRUE))

library(readxl)
library(ggplot2)
library(ggrepel)    # 智能标注避免重叠
library(scales)     # 对数坐标轴格式化

data1 = read_excel("F:\\数据/20230603 - 文库筛选/20230603 - 文库筛选 - 生物膜染色/20250402-文库筛选合集 - 副本.xlsx")  # 读取第一个工作表

# 标记目标区域内的点
data1$in_region <- ifelse(data1$x < 25 & data1$y < 25, "Yes", "No")

# 创建散点图
ggplot(data1, aes(x = x, y = y)) +
  # 1. 散点图层：大小2，黑色描边，按区域填充颜色
  geom_point(aes(fill = in_region), 
             size = 2, 
             shape = 21,         # 带描边的圆形
             color = "black",    # 描边颜色
             stroke = 0.5) +     # 描边粗细
  
  # 2. 为区域内的点添加名称标注（红色点）
  geom_text_repel(
    data = subset(data1, in_region == "Yes"),
    aes(label = 名称),
    color = "black",
    size = 3.5,
    box.padding = 0.5,
    segment.color = "grey",
    max.overlaps = Inf
  ) +
  
  # 3. 填充颜色设置：区域内红色，区域外白色
  scale_fill_manual(values = c("Yes" = "red", "No" = "gray"),
                    guide = "none") + # 不显示图例
  
  # 4. 添加1pt虚线框标出0<x<25且0<y<25的正方形区域
  annotate("rect", 
           xmin = 0, xmax = 25, ymin = 0, ymax = 25,
           fill = NA, 
           color = "black", 
           linetype = "dashed", 
           linewidth = 1) +
  
  # 5. 对数坐标轴设置
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  
  # 6. 坐标轴样式设置
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14, color = "black", face = "bold"),  # 增大字号并加粗
    # 刻度标签设置
    axis.text = element_text(size = 14, color = "black"),    
    # 1pt黑色坐标轴线
    axis.line = element_line(color = "black", linewidth = 1),
    # 1pt黑色刻度线
    axis.ticks = element_line(color = "black", linewidth = 1),
    # 刻度线长度
    axis.ticks.length = unit(2, "mm"),
    # 移除所有网格线
    panel.grid = element_blank(),
    # 保持图形为正方形
    aspect.ratio = 1,
    # 移除图例
    legend.position = "none"
  ) +
  labs(x = "Biofilm enhancement, replicate 1 (% of control)", y = "Biofilm enhancement, replicate 2 (% of control)") +
  # 确保坐标轴从0开始且保持正方形比例
  coord_cartesian(xlim = c(1, 100), ylim = c(1, 100)) +
  # 强制图形为正方形比例
  coord_fixed(ratio = 1)

# 保存图形
ggsave("scatter_plot_with_labels.png", width = 5.2, height = 5.2, dpi = 300)
