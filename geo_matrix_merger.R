#!/usr/bin/env Rscript
# 文件名：geo_matrix_merger.R
# 增强版：基于ID交集的安全合并

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(cli)
  library(dplyr)
})

option_list <- list(
  make_option(c("-g", "--gse_ids"), type = "character",
              help = "逗号分隔的GSE编号列表 (例: GSE12345,GSE67890)"),
  make_option(c("-d", "--data_dir"), type = "character", default = ".",
              help = "矩阵文件存储目录 [默认当前目录]"),
  make_option(c("-o", "--output"), type = "character", default = "merged_matrix.csv",
              help = "输出文件名 [默认: merged_matrix.csv]"),
  make_option(c("-s", "--suffix"), type = "character", default = "_matrix.csv",
              help = "矩阵文件后缀模式 [默认: '_matrix.csv']")
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser)

# 参数验证
if (is.null(args$gse_ids)) stop("必须提供GSE编号列表", call. = FALSE)

tryCatch({
  # 解析GSE列表
  gse_list <- unlist(strsplit(args$gse_ids, ","))
  cli_alert_info("待合并GSE数量: {length(gse_list)}")
  
  # 构建文件路径
  file_paths <- file.path(args$data_dir, paste0(gse_list, args$suffix))
  
  # 第一阶段：获取所有ID的交集
  cli_h2("正在计算ID交集")
  id_list <- lapply(file_paths, function(path) {
    dt <- fread(path, select = 1, header = TRUE)  # 只读取第一列
    if (names(dt)[1] == "V1") setnames(dt, "V1", "ID")  # 处理无名首列
    dt$ID
  })
  
  common_ids <- Reduce(intersect, id_list)
  if (length(common_ids) == 0) stop("没有共同的ID，无法合并", call. = FALSE)
  cli_alert_success("ID交集数量: {length(common_ids)}")
  
  # 第二阶段：基于交集合并数据
  cli_h2("开始合并表达矩阵")
  merged_dt <- lapply(file_paths, function(path) {
    dt <- fread(path, header = TRUE)
    if (names(dt)[1] == "V1") setnames(dt, "V1", "ID")
    
    # 筛选交集ID并设置键
    dt_subset <- dt[ID %chin% common_ids]
    setkey(dt_subset, ID)
    
    # 添加GSE前缀
    gse_id <- gsub(args$suffix, "", basename(path))
    setnames(dt_subset, names(dt_subset)[-1], paste0(gse_id, "_", names(dt_subset)[-1]))
    dt_subset
  }) %>% Reduce(function(x,y) merge(x,y, by = "ID", all = FALSE), .)
  
  # 输出统计信息
  cli_h2("合并结果统计")
  cli_li("原始数据集数量: {length(gse_list)}")
  cli_li("合并前总ID数量: {sum(sapply(id_list, length))}")
  cli_li("最终保留ID数量: {nrow(merged_dt)}")
  cli_li("样本数量: {ncol(merged_dt)-1}")
  
  # 保存结果
  fwrite(merged_dt, args$output)
  cli_alert_success("合并文件已保存至: {args$output}")
  
}, error = function(e) {
  cli_alert_danger("合并失败: {e$message}")
  quit(status = 1)
})
