#!/usr/bin/env Rscript
# 文件名：limma_analysis.R
# 增强版：支持自定义logFC和FDR阈值

suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(ggplot2)
  library(cli)
  library(optparse)
  library(pheatmap)
})

option_list <- list(
  make_option(c("-m", "--matrix"), type = "character",
              help = "批次矫正后的表达矩阵文件路径"),
  make_option(c("-g", "--group_pattern"), type = "character",
              default = "*_group.xls",
              help = "分组文件路径模式 [默认: *_group.xls]"),
  make_option(c("-o", "--output_dir"), type = "character",
              default = "results",
              help = "输出目录 [默认: results]"),
  make_option(c("-f", "--fdr_threshold"), type = "numeric",
              default = 0.05,
              help = "FDR阈值 [默认: 0.05]"),
  make_option(c("-l", "--logfc_threshold"), type = "numeric",
              default = 1,
              help = "最小logFC绝对值 [默认: 1]")
)

parser <- OptionParser(usage = "%prog -m corrected_matrix.csv [options]", 
                       option_list = option_list)
args <- parse_args(parser)

perform_limma_analysis <- function() {
  tryCatch({
    # 参数校验
    if (args$fdr_threshold < 0 | args$fdr_threshold > 1) {
      stop("FDR阈值必须在0-1之间", call. = FALSE)
    }
    if (args$logfc_threshold < 0) {
      stop("logFC阈值必须≥0", call. = FALSE)
    }
    
    # 阶段1：数据加载与预处理 ------------------------------------------------
    cli_h1("数据加载阶段")
    
    # 读取表达矩阵
    cli_alert_info("正在读取表达矩阵: {args$matrix}")
    expr_dt <- fread(args$matrix)
    
    # 检查ID列
    if (!"ID" %in% colnames(expr_dt)) {
      stop("表达矩阵必须包含 'ID' 列作为第一列", call. = FALSE)
    }
    expr_mat <- as.matrix(expr_dt, rownames = "ID")
    
    # 检查样本名
    if (ncol(expr_mat) < 2) {
      stop("表达矩阵必须包含至少一个样本列", call. = FALSE)
    }
    cli_alert_success("矩阵加载完成，维度: {nrow(expr_mat)} 基因 x {ncol(expr_mat)} 样本")
    
    # 加载分组信息
    group_files <- Sys.glob(args$group_pattern)
    if (length(group_files) == 0) {
      stop("未找到分组文件，请检查模式: ", args$group_pattern, call. = FALSE)
    }
    
    cli_alert_info("加载分组文件: {length(group_files)} 个")
    pheno_data <- rbindlist(lapply(group_files, function(f) {
      dt <- fread(f)
      gse_id <- gsub("_group.xls", "", basename(f))
      dt[, sample_id := paste(gse_id, sample_id, sep = "_")]
      dt[, .(sample_id, group)]
    }))
    
    # 数据校验
    matrix_samples <- colnames(expr_mat)
    pheno_samples <- pheno_data$sample_id
    
    # 样本匹配检查
    if (!all(pheno_samples %in% matrix_samples)) {
      stop("存在未在矩阵中定义的分组样本", call. = FALSE)
    }
    
    # 阶段2：差异分析 ------------------------------------------------------
    cli_h1("差异分析阶段")
    
    # 构建设计矩阵
    design <- model.matrix(~ 0 + group, data = pheno_data)
    colnames(design) <- gsub("group", "", colnames(design))
    
    # 拟合线性模型
    cli_alert_info("拟合limma模型...")
    fit <- lmFit(expr_mat[, pheno_samples], design)
    
    # 设置对比矩阵（case vs control）
    contrast_matrix <- makeContrasts(
      contrasts = "case-control",
      levels = design
    )
    
    # 执行差异分析
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    # 提取结果
    full_results <- topTable(fit2, number = Inf, adjust.method = "BH")
    
    # 检查结果列
    required_cols <- c("logFC", "adj.P.Val")
    if (!all(required_cols %in% colnames(full_results))) {
      stop("差异分析结果缺少必要列: ", paste(required_cols, collapse = ", "), call. = FALSE)
    }
    
    # 阶段3：结果筛选与输出 ------------------------------------------------
    cli_h1("结果筛选阶段")
    
    # 应用阈值筛选
    sig_results <- full_results[
      abs(full_results$logFC) >= args$logfc_threshold & 
        full_results$adj.P.Val <= args$fdr_threshold,
    ]
    sig_results$gene=rownames(sig_results)
    cli_alert_success("筛选条件: |logFC| ≥ {args$logfc_threshold} & FDR ≤ {args$fdr_threshold}")
    cli_alert_info("总差异基因数: {nrow(sig_results)}")
    
    # 保存显著结果
    output_file <- file.path(args$output_dir, "differential_genes_significant.csv")
    fwrite(sig_results, output_file)
    cli_alert_success("显著结果已保存: {output_file}")
    
    # 生成火山图
    cli_alert_info("生成火山图...")
    volcano_plot <- ggplot(full_results, aes(logFC, -log10(adj.P.Val))) +
      geom_point(aes(color = adj.P.Val < args$fdr_threshold & abs(logFC) >= args$logfc_threshold), 
                 alpha = 0.6, size = 1.5) +
      geom_hline(yintercept = -log10(args$fdr_threshold), 
                 linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(-args$logfc_threshold, args$logfc_threshold),
                 linetype = "dashed", color = "blue") +
      scale_color_manual(values = c("gray", "red")) +
      labs(title = paste("Volcano Plot (FDR <", args$fdr_threshold, 
                         ", |logFC| >", args$logfc_threshold, ")"),
           x = "log2 Fold Change", 
           y = "-log10(FDR)") +
      theme_bw(base_size = 12) +
      theme(legend.position = "none")
    
    plot_file <- file.path(args$output_dir, "volcano_plot_enhanced.pdf")
    ggsave(plot_file, volcano_plot, width = 8, height = 6, dpi = 300)
    cli_alert_success("火山图已保存: {plot_file}")
    
    cli_h1("生成差异基因热图")
    
    # 提取显著基因的表达矩阵
    sig_genes <- rownames(sig_results)
    sig_expr_mat <- expr_mat[sig_genes, pheno_samples]
    
    # 样本分组注释
    annotation_col <- data.frame(
      Group = pheno_data$group,
      row.names = pheno_data$sample_id
    )
    
    # 绘制热图
    heatmap_file <- file.path(args$output_dir, "differential_genes_heatmap.pdf")
    pdf(heatmap_file)
    pheatmap(
      sig_expr_mat,
      scale = "row",
      annotation_col = annotation_col,
      show_rownames = FALSE,
      show_colnames = FALSE,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      main = paste("Heatmap of Differential Genes (FDR <", args$fdr_threshold, 
                   ", |logFC| >", args$logfc_threshold, ")"),
      silent = TRUE
    )
    dev.off()
    cli_alert_success("热图已保存: {heatmap_file}")
    
  }, error = function(e) {
    cli_alert_danger("分析失败: {e$message}")
    quit(status = 1)
  })
}

perform_limma_analysis()
