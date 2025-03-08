#!/usr/bin/env Rscript
# 文件名：geo_batch_correction.R
# 修复版：增强样本匹配与错误处理

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(sva)
  library(ggplot2)
  library(cli)
})

# 命令行参数配置
option_list <- list(
  make_option(c("-i", "--input_matrix"), type = "character",
              help = "合并后的表达矩阵文件路径"),
  make_option(c("-g", "--group_files"), type = "character",
              help = "逗号分隔的分组文件路径列表"),
  make_option(c("-o", "--output"), type = "character",
              default = "batch_corrected_matrix.csv",
              help = "输出文件名 [默认: batch_corrected_matrix.csv]"),
  make_option(c("-p", "--plot_prefix"), type = "character",
              default = "qc_plot",
              help = "质控图前缀 [默认: qc_plot]")
)

parser <- OptionParser(usage = "%prog -i merged_matrix.csv -g 'group1.xls,group2.xls'", 
                       option_list = option_list)
args <- parse_args(parser)
# 参数验证
mandatory <- c("input_matrix", "group_files")
missing <- mandatory[!mandatory %in% names(args)]
if (length(missing) > 0) {
  stop("缺少必要参数: ", paste(missing, collapse = ", "), call. = FALSE)
}

tryCatch({
  # 阶段1：数据加载 -------------------------------------------------------
  cli_h1("数据加载阶段")
  
  # 1.1 加载表达矩阵
  cli_alert_info("正在读取矩阵文件: {args$input_matrix}")
  expr_dt <- fread(args$input_matrix)
  if (!"ID" %in% names(expr_dt)) {
    stop("矩阵文件必须包含 'ID' 列作为第一列", call. = FALSE)
  }
  expr_mat <- as.matrix(expr_dt, rownames = "ID")
  cli_alert_success("矩阵加载完成，维度: {nrow(expr_mat)} 行 x {ncol(expr_mat)} 列")
  
  # 1.2 加载分组信息
  group_files <- unlist(strsplit(args$group_files, ","))
  cli_alert_info("加载分组文件: {length(group_files)} 个")
  
  pheno_data <- rbindlist(lapply(group_files, function(f) {
    # 文件存在性检查
    if (!file.exists(f)) {
      stop(sprintf("分组文件不存在: %s", f), call. = FALSE)
    }
    
    # 数据加载
    dt <- fread(f)
    
    # 列名校验
    required_cols <- c("sample_id", "group")
    missing_cols <- setdiff(required_cols, names(dt))
    if (length(missing_cols) > 0) {
      stop(sprintf("分组文件 %s 缺少必要列: %s", 
                   basename(f), paste(missing_cols, collapse = ", ")), 
           call. = FALSE)
    }
    
    # 生成带批次前缀的样本ID
    gse_id <- gsub("_group\\.xls$", "", basename(f))
    dt[, sample_id := paste(gse_id, sample_id, sep = "_")]
    dt[, batch := gse_id]
    
    # 返回必要列
    dt[, .(sample_id, group, batch)]
  }))
  
  # 阶段2：数据校验 -------------------------------------------------------
  cli_h1("数据校验阶段")
  
  # 2.1 样本匹配性检查
  matrix_samples <- colnames(expr_mat)
  pheno_samples <- pheno_data$sample_id
  
  # 计算差异样本
  only_in_matrix <- setdiff(matrix_samples, pheno_samples)
  only_in_pheno <- setdiff(pheno_samples, matrix_samples)
  
  # 样本过滤
  matched_samples <- intersect(matrix_samples, pheno_samples)
  if (length(matched_samples) == 0) {
    stop("没有共同样本，请检查样本ID格式", call. = FALSE)
  }
  
  cli_alert_warning("样本差异统计:")
  cli_li("仅存在于矩阵: {length(only_in_matrix)} 个")
  cli_li("仅存在于分组文件: {length(only_in_pheno)} 个")
  cli_li("共同样本: {length(matched_samples)} 个")
  
  # 阶段3：数据预处理 -----------------------------------------------------
  cli_h1("数据预处理阶段")
  
  # 3.1 样本过滤
  cli_alert_info("过滤非共同样本...")
  expr_mat <- expr_mat[, matched_samples, drop = FALSE]
  pheno_data <- pheno_data[sample_id %in% matched_samples]
  
  # 3.2 样本排序
  cli_alert_info("对齐样本顺序...")
  setorder(pheno_data, sample_id)
  expr_mat <- expr_mat[, pheno_data$sample_id, drop = FALSE]
  
  # 3.3 数据转换
  cli_alert_info("检查数据范围...")
  if (max(expr_mat, na.rm = TRUE) > 100) {
    cli_alert_info("执行log2转换")
    expr_mat <- log2(expr_mat + 1)
    if (any(is.infinite(expr_mat))) {
      cli_alert_warning("发现无限值，重置为0")
      expr_mat[is.infinite(expr_mat)] <- 0
    }
  }
  
  # 阶段4：批次矫正 -------------------------------------------------------
  cli_h1("批次效应矫正阶段")
  
  # 4.1 创建模型矩阵
  mod <- model.matrix(~ group, data = pheno_data)
  
  # 4.2 ComBat矫正
  cli_alert_info("执行ComBat矫正...")
  corrected <- ComBat(
    dat = expr_mat,
    batch = pheno_data$batch,
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  
  # 阶段5：结果输出 -------------------------------------------------------
  cli_h1("结果保存阶段")
  
  # 5.1 保存矫正矩阵
  corrected_dt <- data.table(ID = rownames(corrected), corrected)
  fwrite(corrected_dt, args$output)
  cli_alert_success("矫正矩阵已保存: {args$output}")
  
  # 5.2 生成质控图
  cli_alert_info("生成质控图...")
  generate_pca <- function(mat, title) {
    pca <- prcomp(t(mat), scale. = TRUE)
    pca_dt <- data.table(
      PC1 = pca$x[,1],
      PC2 = pca$x[,2],
      Group = pheno_data$group,
      Batch = pheno_data$batch
    )
    
    ggplot(pca_dt, aes(PC1, PC2, color = Group, shape = Batch)) +
      geom_point(size = 3, alpha = 0.8) +
      labs(title = title) +
      theme_bw(base_size = 14)
  }
  
  p_before <- generate_pca(expr_mat, "Before Correction")
  ggsave(paste0(args$plot_prefix, "_before.png"), p_before, width = 8, height = 6)
  
  p_after <- generate_pca(corrected, "After Correction")
  ggsave(paste0(args$plot_prefix, "_after.png"), p_after, width = 8, height = 6)
  
  cli_alert_success("质控图已生成: {args$plot_prefix}_*.png")
  
}, error = function(e) {
  cli_alert_danger("处理失败: {e$message}")
  quit(status = 1)
})
