#!/usr/bin/env Rscript
# 文件名：geo_group_analysis2.R
# 增强版：完善的模糊匹配参数处理
# 用法示例：Rscript geo_group_analysis2.R -i input.xls -c diagnosis --fuzzy "tumor*,*metastasis" --show-all

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(cli)
  library(stringr)
})

# 增强版模糊匹配函数
fuzzy_match <- function(values, patterns, ignore_case = TRUE) {
  # 参数校验
  if (length(patterns) == 0) stop("模糊匹配模式不能为空", call. = FALSE)
  
  # 预处理
  values <- as.character(values)
  orig_case <- unique(values)
  if (ignore_case) {
    values <- tolower(values)
    patterns <- tolower(patterns)
  }
  
  # 生成正则表达式
  regex_patterns <- sapply(patterns, function(p) {
    p <- str_trim(p)
    if (p == "") return(NULL)
    
    # 转换逻辑
    if (str_detect(p, "^\\*.*\\*$")) {    # *contain*
      paste0(".*", str_replace_all(p, "\\*", ""), ".*")
    } else if (str_starts(p, "\\*")) {    # *end
      paste0(str_replace(p, "\\*", ""), "$")
    } else if (str_ends(p, "\\*")) {      # start*
      paste0("^", str_replace(p, "\\*", ""))
    } else {                               # exact match
      paste0("^", p, "$")
    }
  }) %>% unlist() %>% unique()
  
  # 执行匹配
  matched <- rep(FALSE, length(values))
  for (rp in regex_patterns) {
    matched <- matched | str_detect(values, regex(rp))
  }
  
  # 返回原始大小写的结果
  orig_case[matched]
}

# 命令行参数配置
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "输入文件路径"),
  make_option(c("-c", "--column"), type = "character", help = "分组列名"),
  make_option(c("-f", "--fuzzy"), type = "character", 
              help = "模糊匹配模式（支持*通配符，多个模式用逗号分隔）"),
  make_option("--strict", action = "store_false", dest = "ignore_case",
              default = TRUE, help = "禁用大小写忽略"),
  make_option("--show-all", action = "store_true", 
              help = "显示所有可能匹配的值"),
  make_option("--overwrite", action = "store_true", help = "覆盖原文件")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# 参数校验增强
mandatory <- c("input", "column", "fuzzy")
missing_args <- mandatory[!mandatory %in% names(args)]
if (length(missing_args) > 0) {
  stop("缺少必要参数: ", paste(missing_args, collapse = ", "), 
       "\n使用 --help 查看帮助", call. = FALSE)
}

tryCatch({
  # 数据读取
  df <- read_tsv(args$input, col_types = cols(), show_col_types = FALSE)
  
  # 列名校验
  if (!args$column %in% colnames(df)) {
    stop(sprintf("列名 '%s' 不存在，可用列：%s", 
                 args$column, paste(colnames(df), collapse = ", ")), 
         call. = FALSE)
  }
  
  # 解析模糊模式
  patterns <- unlist(str_split(args$fuzzy, ",\\s*"))
  
  # 执行模糊匹配
  unique_vals <- unique(df[[args$column]])
  matched_values <- fuzzy_match(unique_vals, patterns, args$ignore_case)
  
  # 匹配结果校验
  if (length(matched_values) == 0) {
    msg <- sprintf("未找到匹配模式 '%s' 的值", args$fuzzy)
    if (args$show_all) {
      msg <- paste(msg, "\n所有可能的值：\n", paste(unique_vals, collapse = "\n"))
    }
    stop(msg, call. = FALSE)
  }
  
  # 显示匹配结果
  cli_h1("匹配结果")
  cli_alert_info("匹配模式: {args$fuzzy}")
  cli_alert_success("匹配到 {length(matched_values)} 个值：")
  cli_ul(matched_values)
  
  # 创建分组列
  df <- df %>%
    mutate(
      group = if_else(
        .data[[args$column]] %in% matched_values,
        "case", "control",
        missing = "unknown"
      ),
      .after = all_of(args$column)
    )
  
  # 文件保存
  output_file <- if (args$overwrite) args$input else sub(".xls", "_grouped.xls", args$input)
  write_tsv(df, output_file)
  
  # 统计报告
  stats <- df %>%
    count(group) %>%
    mutate(percentage = scales::percent(n / sum(n), accuracy = 0.1))
  
  cli_h1("最终分组分布")
  print(stats)
  cli_alert_success("结果已保存至：{output_file}")
  
}, error = function(e) {
  cli_alert_danger(e$message)
  quit(status = 1)
})
