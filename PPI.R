suppressPackageStartupMessages({
library(tidyverse)
  library(data.table)
library(clusterProfiler) 
library(org.Hs.eg.db)  #小鼠的话，把Hs改成Mm
library(STRINGdb)
library(igraph)
library(ggraph)
  library(cli)
  library(optparse)
  })
option_list <- list(
  make_option(c("-d", "--data_dir"), type = "character",
              help = "limma与wgcna的交集基因文件[limma_wgcna_merge_gene.csv]/来自limma或wgcna单独基因结果文件[differential_genes_significant.csv/WGCNA_ip_gene.csv]"),
  make_option(c("-s", "--STRINGdb"), type = "character",
              help = "STRINGdb R包11版人类蛋白参考文件本地路径。下载链接：https://stringdb-static.org/download/"),
  make_option(c("-t", "--score_threshold"), type = "character",default=400,
              help = "蛋白互作得分阈值默认400。低150，高700，极高900，越高可信度越强")
)
parser <- OptionParser(usage = "%prog -d batch_corrected_matrix.csv -t merged_trait_data.csv -o 输出文件夹路径", 
                       option_list = option_list)
args <- parse_args(parser)
args$data_dir="/home/baiqiang/GEO_STUDY/mul/limma_wgcna_merge_gene.csv"
args$STRINGdb="/home/baiqiang/GEO_STUDY/mul/"
if(is.null(args$STRINGdb)){
  stop("未发现蛋白参考文件，可从https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz下载")
}

cli_h1("针对{basename(args$data_dir)}执行PPI分析")
cli_alert_info("PPI分析，仅适用人类基因和蛋白，如需分析其他物种须分别更改PPI.R中的参考基因数据库和蛋白参考文件输入地址")
#9606是人类，小鼠是10090

string_db <- STRINGdb$new(version="11.5",species=9606,score_threshold=args$score_threshold,input_directory=args$STRINGdb)
cli_alert_success("创建STRINGdb对象成功")
setwd(dirname(args$data_dir))
gene=fread(args$data_dir)
cli_alert_info("矩阵基因数量：{nrow(gene)}")
gene <- gene$gene %>% bitr(fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Hs.eg.db", 
                      drop = T)
cli_alert_success("SYMBOL 2 ENTREZID对象成功")
cli_alert_info("S2E基因保留数量：{nrow(gene)}")
data_mapped <- gene %>% string_db$map(my_data_frame_id_col_names = "ENTREZID", 
                                      removeUnmappedRows = TRUE)
cli_alert_info("成功注释蛋白作用的基因数量：{length(unique(data_mapped$SYMBOL))}")
pdf("PPI_network_STRINGdb.pdf")
string_db$plot_network(data_mapped$STRING_id )
dev.off()
cli_alert_success("生成 STRINGdb R包自带的网络互作图：{dirname(args$data_dir)}/PPI_network_STRINGdb.pdf")
hit<-data_mapped$STRING_id
info <- string_db$get_interactions(hit)
fwrite(info,"PPI_interactions.csv",quote=F,row.names=F)
cli_alert_success("PPI结果生成：{dirname(args$data_dir)}/PPI_interactions.csv")

# 去除游离的互作关系后的network图
# 转换stringID为Symbol，只取前两列和最后一列
links <- info %>%
  mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
  mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
  dplyr::select(from, to , last_col()) %>% 
  dplyr::rename(weight = combined_score)


# 如果links数据框的一个link的from只出现过一次，同时to也只出现一次，则将其去除
links_2 <- links %>% mutate(from_c = count(., from)$n[match(from, count(., from)$from)]) %>%
  mutate(to_c = count(., to)$n[match(to, count(., to)$to)]) %>%
  filter(!(from_c == 1 & to_c == 1)) %>%
  dplyr::select(1,2,3)
# 新的节点数据
nodes_2 <- links_2 %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
# 创建网络图
net_2 <- igraph::graph_from_data_frame(d=links_2,vertices=nodes_2,directed = F)
# 添加必要的参数
igraph::V(net_2)$deg <- igraph::degree(net_2)
igraph::V(net_2)$size <- igraph::degree(net_2)/5
igraph::E(net_2)$width <- igraph::E(net_2)$weight/10
pdf("PPI_network_ggraph.pdf")
ggraph(net,layout = "linear", circular = TRUE)+
  geom_edge_arc(aes(edge_width=width), color = "lightblue", show.legend = F)+
  geom_node_point(aes(size=size), color="orange", alpha=0.7)+
  geom_node_text(aes(filter=deg>5, label=name), size = 5, repel = F)+
  scale_edge_width(range = c(0.2,1))+
  scale_size_continuous(range = c(1,10) )+
  guides(size=F)+
  theme_graph()
dev.off()
cli_alert_success("生成去除游离关系的网络互作图：{dirname(args$data_dir)}/PPI_network_ggraph.pdf")
cli_alert_success("PPI分析完成")
