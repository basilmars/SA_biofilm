# 加载必要的包
library(dplyr)
library(readxl)

# 读取 genomic.xlsx 文件(genomic.gbff文件另存得到)
xlsx_data <- read_excel("./GCF_022870465.1/ncbi_dataset/GCF_022870465.1/genomic.xlsx",
                        col_names = FALSE, trim_ws = FALSE)
lines <- as.character(xlsx_data[[1]])  # 将数据转换为字符向量，每行一个元素

# 目标属性列表
target_attrs <- c("gene", "locus_tag", "old_locus_tag", "GO_function", "GO_component", "GO_process", "product", "protein_id", "translation")

# 找到所有CDS记录的起始行索引（支持所有形式）
cds_start_indices <- which(grepl("^\\s+CDS\\s+([<>]?[0-9]+\\.\\.[<>]?[0-9]+|complement\\([<>]?[0-9]+\\.\\.[<>]?[0-9]+\\))$", lines))

# 初始化存储CDS数据的列表
cds_data <- list()

# 遍历每个CDS记录
for (idx in seq_along(cds_start_indices)) {
  # 提取当前CDS的索引
  current_idx <- cds_start_indices[idx]
  
  # 提取位置（num）
  location <- sub("^\\s+CDS\\s+([<>]?[0-9]+\\.\\.[<>]?[0-9]+|complement\\([<>]?[0-9]+\\.\\.[<>]?[0-9]+\\))", "\\1", lines[current_idx])
  
  # 确定属性块的结束位置（支持所有形式）
  end <- current_idx + 1
  while (end <= length(lines) && !grepl("^\\s*(gene|tRNA|rRNA|misc_feature|CDS)\\s+(complement\\([<>]?[0-9]+\\.\\.[<>]?[0-9]+\\)|[<>]?[0-9]+\\.\\.[<>]?[0-9]+)", lines[end])) {
    end <- end + 1
  }
  end <- end - 1
  
  # 收集所有属性行
  attribute_lines <- lines[(current_idx + 1):end]

  # 输出当前CDS和下一个CDS之间的行（如果下一个CDS也在目标记录中）
  if (idx < length(cds_start_indices)) {
    next_idx <- cds_start_indices[idx + 1]
    next_location <- sub("^\\s+CDS\\s+([<>]?[0-9]+\\.\\.[<>]?[0-9]+|complement\\([<>]?[0-9]+\\.\\.[<>]?[0-9]+\\))", "\\1", lines[next_idx])
    if (location %in% target_locations && next_location %in% target_locations) {
      cat("Lines between", location, "and", next_location, "(index", end + 1, "to", next_idx, "):\n", file = stdout())
      print(lines[(end + 1):next_idx])
      cat("Lines between", location, "and", next_location, "(index", end + 1, "to", next_idx, "):\n", file = debug_file)
      writeLines(lines[(end + 1):next_idx], debug_file)
      cat("\n", file = debug_file)
    }
  }
  
  # 初始化属性列表
  attr_list <- list()
  current_attr <- NULL
  
  # 逐行解析属性
  for (line in attribute_lines) {
    if (grepl("^\\s*/", line)) {
      # 新属性开始
      attr_name <- sub("=.*", "", sub("^\\s*/", "", line))
      value_part <- sub("^\\s*/\\w+=\"?", "", line)
      value_part <- sub("\"$", "", value_part)
      current_attr <- attr_name
      if (current_attr %in% names(attr_list)) {
        attr_list[[current_attr]] <- c(attr_list[[current_attr]], value_part)
      } else {
        attr_list[[current_attr]] <- value_part
      }
    } else if (!is.null(current_attr)) {
      # 多行属性的续行
      value_part <- sub("^\\s+", "", line)
      value_part <- sub("\"$", "", value_part)
      attr_list[[current_attr]] <- paste(attr_list[[current_attr]], value_part, collapse = " ")
    }
  }
  
  # 创建CDS记录行
  cds_row <- list(num = location)
  for (attr_name in target_attrs) {
    if (attr_name %in% names(attr_list)) {
      cds_row[[attr_name]] <- paste(attr_list[[attr_name]], collapse = ", ")
    } else {
      cds_row[[attr_name]] <- NA
    }
  }
  
  # 将记录添加到cds_data
  cds_data <- append(cds_data, list(cds_row))
}

# 将所有CDS记录合并成数据框
df <- bind_rows(cds_data)

#输出
write.csv(df, "CDS_extracted.csv", row.names = FALSE)
