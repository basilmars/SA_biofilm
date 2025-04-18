##准备用于功能注释的表格
##将参考基因组的fasta格式的CDS序列输入eggNOG-Mapper在线工具，进行转换，下载excel文件

#读取eggNOG-Mapper转换的注释（http://eggnog-mapper.embl.de/job_status?jobname=MM_7tq7dp7q）
library(readxl)
library(dplyr)

eggNOG <- read_excel("E://R data/20241122/MM_acw2c_o8.emapper.annotations.xlsx")

#重命名并移除空白和标题行
colnames(eggNOG) <- eggNOG[2,]
eggNOG <- eggNOG[-c(1,2),]

# 移除倒数3行说明，保留前 n-3 行
eggNOG <- head(eggNOG, n = nrow(eggNOG) - 3)

# 创建 COG 字母到功能层级的映射表
cog_function_hierarchy <- data.frame(
  COG_category = c("J", "A", "K", "L", "B", "D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O", "X",
                   "C", "G", "E", "F", "H", "I", "P", "Q", "R", "S"),
  Function = c("Translation, ribosomal structure and biogenesis",
               "RNA processing and modification",
               "Transcription",
               "Replication, recombination and repair",
               "Chromatin structure and dynamics",
               "Cell cycle control, cell division, chromosome partitioning",
               "Nuclear structure",
               "Defense mechanisms",
               "Signal transduction mechanisms",
               "Cell wall/membrane/envelope biogenesis",
               "Cell motility",
               "Cytoskeleton",
               "Extracellular structures",
               "Intracellular trafficking, secretion, and vesicular transport",
               "Posttranslational modification, protein turnover, chaperones",
               "Mobilome: prophages, transposons",
               "Energy production and conversion",
               "Carbohydrate transport and metabolism",
               "Amino acid transport and metabolism",
               "Nucleotide transport and metabolism",
               "Coenzyme transport and metabolism",
               "Lipid transport and metabolism",
               "Inorganic ion transport and metabolism",
               "Secondary metabolites biosynthesis, transport and catabolism",
               "General function prediction only",
               "Function unknown"),
  Function_Group = c(rep("Information Storage and Processing", 5),
                     rep("Cellular Processes and Signaling", 11),
                     rep("Metabolism", 8),
                     rep("Poorly Characterized", 2))
)

# 查看映射表
head(cog_function_hierarchy)

# 将 COG 分类与功能层级数据合并
cog_data_with_hierarchy <- left_join(eggNOG, cog_function_hierarchy, by = "COG_category")

#加载原cds文件，以匹配Locus_Tag
library(Biostrings)

# 读取 .fna 文件
cds_file <- "E://R data/20241122/ncbi_dataset/GCA_022870465.1/cds_from_genomic.fna"
cds_sequences <- readDNAStringSet(cds_file)

# 转换为只包含ID的数据框
cds_df <- data.frame(
  ID = names(cds_sequences),
  stringsAsFactors = FALSE
)

# 安装并加载必要的包
library(stringr)
library(tidyverse)

# 将结果转换为数据框
cds_df_id <- data.frame(
  query = str_match(cds_df$ID, "(.*?)\\ ")[, 2],
  Gene = str_match(cds_df$ID, "\\[gene=(.*?)\\]")[, 2],
  Locus.Tag = str_match(cds_df$ID, "\\[locus_tag=(.*?)\\]")[, 2],
  Protein = str_match(cds_df$ID, "\\[protein=(.*?)\\]")[, 2],
  Protein_ID = str_match(cds_df$ID, "\\[protein_id=(.*?)\\]")[, 2],
  Location = str_match(cds_df$ID, "\\[location=(.*?)\\]")[, 2],
  GB_Key = str_match(cds_df$ID, "\\[gbkey=(.*?)\\]")[, 2],
  stringsAsFactors = FALSE)

##准备用于分析的差异基因数据

# 读取基因表达矩阵和分组信息
#设置文件路径
dir <- "E://R data/20241122/"

#读取文件
HNP1_vs_ctrl <- read.csv(file.path(dir,"H_vs_S.rpkm.DEG.xls"),
                         row.names = 1,stringsAsFactors = F, header = T, sep = "\t")

# 筛选每个比较组中显著上调的基因
degs_HNP1_vs_ctrl <- HNP1_vs_ctrl %>%
  filter(logFC > 1, Pvalue < 0.05) %>%
  rownames()

# 筛选显著下调的基因
degs_down_HNP1_vs_ctrl <- HNP1_vs_ctrl %>%
  filter(logFC < -1, Pvalue < 0.05) %>%
  rownames()

#将基因的Locus.tag转换为对应的query
up_HNP1_vs_ctrl <- cds_df_id$query[match(degs_HNP1_vs_ctrl, cds_df_id$Locus.Tag)]
down_HNP1_vs_ctrl <- cds_df_id$query[match(degs_down_HNP1_vs_ctrl, cds_df_id$Locus.Tag)]

## COG进行分析

gene_list_up <- up_HNP1_vs_ctrl
gene_list_down <- down_HNP1_vs_ctrl

# 提取和匹配基因 ID
cog_up <- eggNOG %>% filter(query %in% gene_list_up) 
cog_up <- cog_up[,c("query", "COG_category", "Description")]

cog_down <- eggNOG %>% filter(query %in% gene_list_down) 
cog_down <- cog_down[,c("query", "COG_category", "Description")]


library(tidyr)
library(dplyr)

# 使用 separate_rows 将 COG_category 按照单个字母拆分
cog_up_separated <- cog_up %>%
  separate_rows(COG_category, sep = "") %>%
  filter(COG_category != "")

cog_down_separated <- cog_down %>%
  separate_rows(COG_category, sep = "") %>%
  filter(COG_category != "")

# 统计上调基因的 COG 分类
cog_up_count <- cog_up_separated %>%
  group_by(COG_category) %>%
  summarise(up_count = n())

# 统计下调基因的 COG 分类
cog_down_count <- cog_down_separated %>%
  group_by(COG_category) %>%
  summarise(down_count = n())

# 合并上调和下调基因的 COG 分类计数
cog_counts <- full_join(cog_up_count, cog_down_count, by = "COG_category") %>%
  replace_na(list(up_count = 0, down_count = 0))  # 将 NA 替换为 0

library(ggplot2)
library(dplyr)

# 确保 cog_counts 的 COG_category 列是 character 类型 
cog_counts$COG_category <- as.character(cog_counts$COG_category)

# 首先，确保 cog_counts 数据框中合并了 cog_function_hierarchy 中的 Function 列
cog_counts <- left_join(cog_counts, cog_function_hierarchy, by = "COG_category")

# 将 COG_category 中的 - 替换为 "Unknown"
cog_counts <- cog_counts %>%
  mutate(COG_category = ifelse(COG_category == "-", "Unknown", COG_category))

# 然后，将 Y 轴按照 Function 列进行排序
cog_counts$COG_category <- factor(cog_counts$COG_category, levels = rev(c(cog_function_hierarchy$COG_category,"Unknown")))

# 绘制对称的 COG 分类柱状图
ggplot(cog_counts, aes(y = COG_category)) +
  geom_bar(aes(x = -down_count), stat = "identity", fill = "#75b1d3") +  # 左侧为下调基因
  geom_bar(aes(x = up_count), stat = "identity", fill = "#e59415") +      # 右侧为上调基因
  
  # 在 x=0 处添加 Y 轴的黑色实线
  geom_vline(xintercept = 0, color = "black") +
  
  labs(title = "S. aureus + HNP1 vs S. aureus", x = "Number of genes", y = "COG Category") +

  # 设置轴的样式
  theme_minimal() +
  theme(
    axis.line.y = element_blank(),                      # 移除默认的 Y 轴线
    axis.line.x = element_blank(),                      # 移除默认的 X 轴线
    #panel.grid = element_blank(),                       # 去掉网格线
    axis.ticks.x = element_line(color = "black"),       # X轴刻度线
    axis.text.x.top = element_text(),                   # 启用顶部的 X 轴标签
    axis.line.x.top = element_line(color = "black")     # 顶部 X 轴黑色实线
  ) +
  scale_x_continuous(position = "top",                  # 将 X 轴放在上方，设置 X 轴的间隔为 5
                     breaks = seq(-5 * ceiling(max(cog_counts$down_count) / 5),  # 负最大值
                                  5 * ceiling(max(cog_counts$up_count) / 5),    # 正最大值
                                  by = 10))



