library(optparse)
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(janitor)
library(GEOquery)  # GEO数据下载与解析
library(limma)     # 差异表达分析
library(affy)      # Affymetrix芯片数据处理
# 定义参数
option_list <- list(
  make_option(c("-n", "--GSEnumber"), type = "character",help = "GSE编号(例如：GSE118370)"),
  make_option(c("-d", "--storage_dir"), type = "character",default=getwd(),help = "GEO矩阵和注释存储路径"),
  make_option(c("-a", "--GEO_variance_analysis"), action = "store_true", default = FALSE,help = "是否执行差异分析"),
  make_option(c("-e", "--enrich_analysis"), action = "store_true", default = FALSE,help = "是否执行富集分析")
  )
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

if (is.null(args$GSEnumber)) {
  print_help(parser)
  stop("必须指定--GSEnumber")
}
GSEnumber=args$GSEnumber
output_dir=args$storage_dir
setwd(output_dir)
series_matrix_path=paste0(GSEnumber,"_series_matrix.txt.gz")
f_soft=paste0(GSEnumber,"_family.soft.gz")


eSet <- getGEO(filename=series_matrix_path, destdir=output_dir,AnnotGPL = T,getGPL = F)
    # 提取临床信息
    pd <- pData(eSet)
    # 自动识别正常和肿瘤样本的分组脚本
    # 输入：样本信息表 pd（行名为样本ID，如GSM3325818）
    # 输出：分组向量 group 和样本信息增强后的 pd_with_group
    
    # --------------------- 方法1：基于标题特征识别 ---------------------
    # 通过标题（title）中的关键词判断样本类型
    pd$group <- ifelse(grepl("normal", pd$title, ignore.case = TRUE), "Normal", 
                       ifelse(grepl("adenocarcinoma|tumor", pd$title, ignore.case = TRUE), "Tumor", NA))
    
    # --------------------- 方法2：基于组织特征列识别（更准确）---------------------
    # 直接使用 characteristics_ch1 列的明确标注
    pd$group <- ifelse(grepl("tissue: normal tissue", pd$characteristics_ch1), "Normal",
                       ifelse(grepl("tissue: tumor tissue", pd$characteristics_ch1), "Tumor", NA))
    
    # --------------------- 结果验证 ---------------------
    # 检查分组结果
    table(pd$group, useNA = "always")
    
    # 输出带分组的样本信息（包含样本ID和分组）
    pd_with_group <- data.frame(
      sample_id = rownames(pd),
      group = pd$group,
      pd[, c("title", "characteristics_ch1")]  # 保留关键列用于人工复核
    )
    
    # 查看前6行
    head(pd_with_group)
    
    # --------------------- 与表达矩阵对齐 ---------------------
    # 确保表达矩阵列名与 pd 样本ID一致（重要！）
    exp_matrix <- exp_matrix[, rownames(pd)]  # exp_matrix 是你的表达矩阵
    
    # 创建分组因子向量（用于后续差异分析）
    group <- factor(pd$group, levels = c("Normal", "Tumor"))
    
    # 检查样本顺序是否一致
    identical(colnames(exp_matrix), rownames(pd))
    write.table(pd_with_group,paste0(GSEnumber,"_group.xls"),
                sep = "\t",na = "",row.names = F, col.names = T, quote = F)
    
    exp_ori0 <- exprs(eSet)
    identical(colnames(exp_ori0),rownames(pd))
    if( ! identical(rownames(pd),colnames(exp_ori0)) ) {
      exp_ori0 = exp_ori0[,match(rownames(pd),colnames(exp_ori0))] 
    }
    
    # 读入探针注释
    gpl <- parseGEO(fname = f_soft, GSElimits = NULL,AnnotGPL = F, 
                    getGPL = F)
    gpl_ori <- gpl@gpls[[1]]@dataTable@table # 获取平台注释表格
    # gpl_ori <- gpl@gpls[[1]]@dataTable %>% GEOquery::Table()
    # gpl_ori <- gpl@gpls$GPL16570@dataTable %>% GEOquery::Table()
    gpl_number <- gpl@header$platform_id # 获取平台编号
    colnames(gpl_ori)=tolower(colnames(gpl_ori))
    gpl<-gpl_ori[,c("id","gene symbol")]
    write.table(gpl, paste0(GSEnumber,"_anno_probeid2symbol.xls"),
                sep = "\t",na = "",row.names = F, col.names = T, quote = F)
    gpl$`gene symbol`<-data.frame(sapply(gpl$`gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]
    gpl$`gene symbol` <- trimws(as.character(gpl$`gene symbol`))
    exp<-as.data.frame(exp_ori0)
    exp$id<-rownames(exp)
    exp_symbol<-merge(exp,gpl,by="id")
    exp_symbol<-na.omit(exp_symbol)
    exp_matrix <- exp_symbol[, -c(1, ncol(exp_symbol))]
    exp_unique_matrix <- avereps(exp_matrix, ID = exp_symbol$`gene symbol`)
    exp_unique <- as.data.frame(exp_unique_matrix) %>%
      tibble::rownames_to_column(var = "gene_symbol")

    rownames(exp_unique) <- exp_unique$gene_symbol
    exp_matrix <- as.matrix(exp_unique[, -1])
write.csv(exp_matrix,paste0(GSEnumber,"_matrix.csv"))
#-----------------------------------------------------------------------#analysis
    #exp_matrix矩阵，pd_with_group[,c(1,2)]分组
if (args$GEO_variance_analysis==T){
      #判断矩阵是否log2转换-----------------------------------------------------------
  is_log2_transformed <- function(exp_matrix, 
                                  median_threshold = 30, 
                                  max_threshold = 1000,
                                  prob = 0.99) {
    # 核心判断逻辑 -------------------------------------------------------------
    matrix_values <- as.numeric(as.matrix(exp_matrix))
    q99 <- quantile(matrix_values, probs = prob, na.rm = TRUE)
    med <- median(matrix_values, na.rm = TRUE)
    has_negative <- any(matrix_values < 0)
    
    # 判断条件（可扩展）
    condition_logged <- (q99 < max_threshold) && 
      (med < median_threshold) && 
      (!has_negative)
    
    cat("Matrix diagnostic:\n",
        "  -", prob*100, "% quantile:", round(q99, 2), "\n",
        "  - Median:", round(med, 2), "\n",
        "  - Contains negative values:", has_negative, "\n",
        "  - Suggested log2 status:", ifelse(condition_logged, "Likely LOGGED", "Likely UNLOGGED"), "\n")
    return(condition_logged)
  }
  auto_log2_transform <- function(exp_matrix, 
                                  offset = 1, 
                                  verbose = TRUE) {
    # 检查是否为矩阵或数据框
    if (!is.matrix(exp_matrix) && !is.data.frame(exp_matrix)) {
      stop("输入必须是矩阵或数据框")
    }
    
    # 判断是否需要转换
    needs_log <- !is_log2_transformed(exp_matrix)  # 使用之前定义的判断函数
    
    if (needs_log) {
      if (verbose) cat("\n检测到数据可能需要 log2 转换...\n")
      
      # 处理非正数值（添加偏移量）
      min_value <- min(exp_matrix, na.rm = TRUE)
      if (min_value <= 0) {
        if (verbose) cat("检测到非正数值，添加偏移量", offset, "进行转换\n")
        exp_matrix <- exp_matrix + offset
      }
      
      # 执行 log2 转换
      if (verbose) cat("正在执行 log2 转换...\n")
      exp_matrix <- log2(exp_matrix)
    } else {
      if (verbose) cat("\n数据已符合 log2 转换特征，无需处理\n")
    }
    
    return(exp_matrix)
  }
  exp_matrix2=auto_log2_transform(exp_matrix)
  #---------------------------------------------分析
  group_list <- pd_with_group$group
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  contrast_matrix<-makeContrasts(Tumor-Normal,levels = design)
  ##PCA分析
  library(FactoMineR) 
  library(factoextra) 
  dat=as.data.frame(t(exp_matrix2))
  dat.pca <- PCA(dat, graph = FALSE)
  pca_plot <- fviz_pca_ind(dat.pca,
                           geom.ind = "point", # show points only (nbut not "text")
                           col.ind = group_list, # color by groups
                           palette = c("#00AFBB", "#E7B800"),
                           addEllipses = TRUE, # Concentration ellipses
                           legend.title = "Groups"
  )
  ggsave(plot = pca_plot,filename = paste0(GSEnumber,"_PCA.pdf"), dpi = 300)
  ## 差异分析
  fit <- lmFit(exp_matrix2,design)
  
  fit2 <- contrasts.fit(fit,contrast_matrix) 
  fit2 <- eBayes(fit2) 
  DiffEG<-topTable(fit2, coef=1, n=Inf) %>% na.omit()
  write.table(DiffEG,paste0(GSEnumber,"_DiffEG.xls"),
              sep = "\t",na = "",row.names = F, col.names = T, quote = F)
  library(ggplot2)
  library(ggrepel)  # 用于防止标签重叠
  
  # 添加基因名列（假设行名为基因）
  DiffEG$Gene <- rownames(DiffEG)
  
  # 定义显著性阈值
  padj_cutoff <- 0.05
  logFC_cutoff <- 1
  
  # 标记显著基因
  DiffEG <- DiffEG %>%
    mutate(
      Significance = case_when(
        adj.P.Val < padj_cutoff & logFC > logFC_cutoff ~ "Upregulated",
        adj.P.Val < padj_cutoff & logFC < -logFC_cutoff ~ "Downregulated",
        TRUE ~ "Not significant"
      ),
      log10_padj = -log10(adj.P.Val)
    )
  
  # 选择 top10 显著基因标注
  top_genes <- DiffEG %>%
    arrange(adj.P.Val) %>%
    head(10)
  
  # 创建火山图
  #颜色编码：
  #红色：显著上调 (adj.P.Val < 0.05 & logFC > 1)
  #蓝色：显著下调 (adj.P.Val < 0.05 & logFC < -1)
  #灰色：不显著基因
  #动态阈值线：
  #垂直虚线表示 log2FC 阈值 (±1)
  #水平虚线表示显著性阈值 (-log10(0.05))
  volcano_plot <- ggplot(DiffEG, aes(x = logFC, y = log10_padj, 
                                     color = Significance)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Downregulated" = "blue", 
                                  "Upregulated" = "red",
                                  "Not significant" = "grey60")) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", 
               color = "black", linewidth = 0.5) +
    geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), 
               linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_text_repel(data = top_genes, 
                    aes(label = Gene),
                    box.padding = 0.5,
                    max.overlaps = 20,
                    segment.color = "grey50") +
    labs(x = "log2 Fold Change",
         y = "-log10(Adjusted P-value)",
         title = "Differentially Expressed Genes") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")
  ggsave(paste0(GSEnumber,"_volcano_plot.pdf"), plot = volcano_plot, 
         width = 8, height = 6, dpi = 300)
  #绘制热图
  cg=rownames(DiffEG[DiffEG$adj.P.Val<0.05,])
  n=exp_matrix2[cg,]
  annotation_col=pd_with_group%>%select(group)
  library(pheatmap)
  if (length(cg)>30){
  pdf(paste0(GSEnumber,"_heatmap.pdf"))
  heatmap_plot <- pheatmap(n,
                           show_colnames =F,
                           show_rownames = F,
                           annotation_col=annotation_col,
                           scale = "row")
  dev.off()}else{
    pdf(paste0(GSEnumber,"_heatmap.pdf"))
    heatmap_plot <- pheatmap(n,
                             show_colnames =F,
                             annotation_col=annotation_col,
                             scale = "row")
    dev.off()
  }
}
###富集分析
if(args$enrich_analysis=T){
  if(exists("DiffEG")){
    library(clusterProfiler)
    library(ggthemes)
    library(org.Hs.eg.db)
    library(dplyr)
    library(ggplot2)
    library(stringr)
    library(enrichplot)
    s2e = bitr(DiffEG$Gene, 
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)#人类
    nrow(DiffEG)
    deg = inner_join(DiffEG,s2e,by=c("Gene"="SYMBOL"))
    nrow(deg)
    print(paste0("S2E过程的基因丢失:",nrow(DiffEG)-nrow(deg)))
    print(paste0("S2E过程的丢失占比:",(nrow(DiffEG)-nrow(deg))/nrow(DiffEG)))
    gene_diff = deg$ENTREZID
    ekk <- enrichKEGG(gene = gene_diff,organism = 'hsa')
    ekk <- setReadable(ekk,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
    ego <- enrichGO(gene = gene_diff,OrgDb= org.Hs.eg.db,
                    ont = "ALL",readable = TRUE)
    pdf(paste0(GSEnumber,"_GO_bar.pdf"),height = 15)
    barplot(ego, split = "ONTOLOGY") + 
      facet_grid(ONTOLOGY ~ ., space = "free_y",scales = "free_y") 
    dev.off()
    pdf(paste0(GSEnumber,"_KEGG_bar.pdf"))
    barplot(ekk)
    dev.off()
    pdf(paste0(GSEnumber,"_GO_dot.pdf"),height = 15)
    dotplot(ego, split = "ONTOLOGY") + 
      facet_grid(ONTOLOGY ~ ., space = "free_y",scales = "free_y") 
    dev.off()
    pdf(paste0(GSEnumber,"_KEGG_dot.pdf"))
    dotplot(ekk)
    dev.off()
write.table(as.data.frame(ekk), paste0(GSEnumber,"_KEGG.xls"),sep = "\t",na = "",row.names = F, col.names = T, quote = F)
write.table(as.data.frame(ego), paste0(GSEnumber,"_GO.xls"),sep = "\t",na = "",row.names = F, col.names = T, quote = F)
  }else{
    stop("进行富集分析前，需先进行差异分析。请检查-a是否为T")
  }
}
