#!/usr/bin/env Rscript
# 文件名：merge_pheno_data.R
# 功能：按列名合并多个分组文件

suppressPackageStartupMessages({
  library(data.table)
  library(cli)
  library(optparse)
})

# 命令行参数配置
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character",
              help = "输入文件目录"),
  make_option(c("-o", "--output_file"), type = "character",
              help = "输出文件路径"),
  make_option(c("-g", "--gse_ids"), type = "character",
              help = "逗号分隔的GSE编号列表"),
  make_option("--id_col", type = "character", default = "sample_id",
              help = "样本ID列名 [默认: sample_id]"),
  make_option("--group_col", type = "character", default = "group",
              help = "分组列名 [默认: group]")
)

parser <- OptionParser(usage = "%prog -i /input/dir -o merged.csv -g GSE1,GSE2", 
                       option_list = option_list)
args <- parse_args(parser)

# 主函数
merge_pheno_files <- function() {
  tryCatch({
    # 参数校验
    mandatory <- c("input_dir", "output_file", "gse_ids")
    missing <- mandatory[!mandatory %in% names(args)]
    if (length(missing) > 0) stop("缺少必要参数: ", paste(missing, collapse = ", "))
    
    gse_list <- unlist(strsplit(args$gse_ids, ","))
    if (length(gse_list) == 0) stop("未指定GSE编号")
    
    # 合并数据
    merged_data <- rbindlist(lapply(gse_list, function(gse_id) {
      file_path <- file.path(args$input_dir, paste0(gse_id, "_group.xls"))
      if (!file.exists(file_path)) stop("文件不存在: ", file_path)
      
      dt <- fread(file_path)
      
      # 列名校验
      required_cols <- c(args$id_col, args$group_col)
      missing_cols <- setdiff(required_cols, colnames(dt))
      if (length(missing_cols) > 0) {
        stop(paste0("文件 ", basename(file_path), " 缺少列: ", paste(missing_cols, collapse = ", ")))
      }
      
      # 添加前缀并选择列
      dt[, (args$id_col) := paste(gse_id, get(args$id_col), sep = "_")]
      dt[, .(sample_id = get(args$id_col), group = get(args$group_col))]
    }))
    
    # 保存结果
    fwrite(merged_data, args$output_file)
    cli_alert_success("成功合并 {length(gse_list)} 个数据集到: {args$output_file}")
    
  }, error = function(e) {
    cli_alert_danger("合并失败: {e$message}")
    quit(status = 1)
  })
}

merge_pheno_files()
