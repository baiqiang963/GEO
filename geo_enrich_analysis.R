###富集分析
###默认仅适用于人类
suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
    library(cli)
    library(clusterProfiler)
    library(ggthemes)
    library(org.Hs.eg.db)#小鼠的话，把Hs改成Mm
    library(dplyr)
    library(ggplot2)
    library(stringr)
    library(enrichplot)})
option_list <- list(
  make_option(c("-d", "--data_dir"), type = "character",
help = "limma与wgcna的交集基因文件[limma_wgcna_merge_gene.csv]/来自limma或wgcna单独基因结果文件[differential_genes_significant.csv/WGCNA_ip_gene.csv]")
)
parser <- OptionParser(usage = "%prog -d batch_corrected_matrix.csv -t merged_trait_data.csv -o 输出文件夹路径", 
                       option_list = option_list)
args <- parse_args(parser)
if (is.null(args$data_dir)){
  stop("富集分析输入的矩阵为空，请检查shell脚本设置")
}
cli_h1("针对{basename(args$data_dir)}执行富集分析")
cli_alert_info("GO与KEGG富集分析，仅适用人类基因，如需分析其他物种须在geo_enrich_analysis.R中更改参考基因数据库")
setwd(dirname(args$data_dir))
DiffEG=fread(args$data_dir)
    s2e = bitr(DiffEG$gene, 
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)#人类
    nrow(DiffEG)
    deg = inner_join(DiffEG,s2e,by=c("gene"="SYMBOL"))
    nrow(deg)
    print(paste0("S2E过程的基因丢失:",nrow(DiffEG)-nrow(deg)))
    print(paste0("S2E过程的丢失占比:",(nrow(DiffEG)-nrow(deg))/nrow(DiffEG)))
    gene_diff = deg$ENTREZID
    ekk <- enrichKEGG(gene = gene_diff,organism = 'hsa')
    ekk <- setReadable(ekk,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
    ego <- enrichGO(gene = gene_diff,OrgDb= org.Hs.eg.db,
                    ont = "ALL",readable = TRUE)
    pdf(paste0("GO_bar.pdf"),height = 15)
    barplot(ego, split = "ONTOLOGY") + 
      facet_grid(ONTOLOGY ~ ., space = "free_y",scales = "free_y") 
    dev.off()
    pdf(paste0("KEGG_bar.pdf"))
    barplot(ekk)
    dev.off()
    pdf(paste0("GO_dot.pdf"),height = 15)
    dotplot(ego, split = "ONTOLOGY") + 
      facet_grid(ONTOLOGY ~ ., space = "free_y",scales = "free_y") 
    dev.off()
    pdf("KEGG_dot.pdf")
    dotplot(ekk)
    dev.off()
    write.table(as.data.frame(ekk), "KEGG.xls",sep = "\t",na = "",row.names = F, col.names = T, quote = F)
    write.table(as.data.frame(ego), "GO.xls",sep = "\t",na = "",row.names = F, col.names = T, quote = F)
    