#!/usr/bin/env Rscript
# 文件名：geo_auto_group.R
# 功能：自动筛选二分类列并保存分组信息
# 用法示例：
#cd GSE12345保存路径
#Rscript ./GEO/geo_auto_group.R -n GSE12345

suppressPackageStartupMessages({
  library(optparse)
  library(GEOquery)
  library(dplyr)
  library(cli)
})

# 定义命令行参数
option_list <- list(make_option(c("-n", "--GSEnumber"), type = "character", help = "GSE编号(例如：GSE118370)"))

# 解析参数
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# 参数验证
if (is.null(args$GSEnumber)) {
  print_help(parser)
  stop("必须提供GSE编号", call. = FALSE)
}
args$storage_dir=getwd()
tryCatch({
  # 设置工作目录
  setwd(args$storage_dir)
  
  # 构建文件路径
  series_matrix <- paste0(args$GSEnumber, "_series_matrix.txt.gz")
  output_file <- paste0(args$GSEnumber, "_group.xls")
  
  # 检查文件是否存在
  if (!file.exists(series_matrix)) {
    stop(paste("未找到系列矩阵文件:", series_matrix, 
               "请确认：", 
               "1. 文件命名符合 GSEXXXX_series_matrix.txt.gz 格式",
               "2. 文件储存位置", getwd(),
               sep = "\n"))
  }
  
  # 加载GEO数据
  eSet <- getGEO(filename = series_matrix, 
                 destdir = args$storage_dir,
                 AnnotGPL = TRUE,
                 getGPL = FALSE)
  
  pd <- pData(eSet)
  pd$sample_id=rownames(pd)
  
  # 保存结果
  write.table(pd, file = output_file,
              sep = "\t", na = "",
              row.names = FALSE, quote = FALSE)
  
  cli_alert_success(
    "成功保存")
}, error = function(e) {
  cli_alert_danger(paste("运行失败:", e$message))
  quit(status = 1)
})
