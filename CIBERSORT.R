suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(cli)
  library(optparse)
  library(CIBERSORT)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
  library(ggplot2)
  library(tinyarray)
})
cli_h1("免疫浸润分析")
option_list <- list(
  make_option(c("-d", "--data_dir"), type = "character",
              help = "批次矫正后的矩阵[batch_corrected_matrix.csv]")
)
parser <- OptionParser(usage = "%prog -d batch_corrected_matrix.csv -t merged_trait_data.csv -o 输出文件夹路径", 
                       option_list = option_list)
args <- parse_args(parser)
expr_dt <- fread(args$data_dir)
expr_mat <- as.matrix(expr_dt, rownames = "ID")
# 伪计数
c <- 1
expr_mat<-2^expr_mat-c
cli_alert_success("矩阵log2计数还原为原始计数")
if (any(expr_mat < 0)){
  stop("还原后的矩阵存在负数，可能原先的矩阵计数并非log2转换而来，而是使用了Z-score标准化，请注意检查")
}
exp2 = as.data.frame(expr_mat)
exp2 = rownames_to_column(exp2)
setwd(dirname(args$data_dir))
write.table(exp2,file = "exp_CIBERSORT.txt",row.names = F,quote = F,sep = "\t")

lm22f = system.file("extdata", "LM22.txt", package = "CIBERSORT")
TME.results = cibersort(lm22f, 
                        "exp_CIBERSORT.txt" , 
                        perm = 1000, 
                        QN = T)
fwrite(TME.results,"TME_results.csv",row.names = T,quote = F)

re <- TME.results[,-(23:25)]
exp=expr_mat
cli_alert_info("生成免疫细胞丰度热图。。。")
pdf("immune_heatmap.pdf")
k <- apply(re,2,function(x) {sum(x == 0) < nrow(TME.results)/2})
table(k)
re2 <- as.data.frame(t(re[,k]))
Group = str_sub(colnames(exp),1,str_length(colnames(exp))-2)
table(Group)
an = data.frame(group = Group,
                row.names = colnames(exp))
pheatmap(re2,scale = "row",
         show_colnames = F,
         cluster_cols = F,
         annotation_col = an,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()
cli_alert_success("生成免疫细胞丰度热图成功：{dirname(args$data_dir)}/immune_heatmap.pdf")
cli_alert_info("生成免疫细胞柱状图。。。")
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
dat <- re %>% 
  as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  mutate(group = Group) %>% 
  gather(key = Cell_type,value = Proportion,-Sample,-group) %>% 
  arrange(group)

dat$Sample = factor(dat$Sample,ordered = T,levels = unique(dat$Sample)) #定横坐标顺序
# 先把group排序，然后将sample设为了因子，确定排序后的顺序为水平，所以两图的顺序是对应的。
dat2 = data.frame(a = 1:ncol(exp),
                  b = 1,
                  group = sort(Group)) 

p1 = ggplot(dat2,aes(x = a, y = b)) + 
  geom_tile(aes(fill = group)) + 
  scale_fill_manual(values = mypalette(22)[1:length(unique(Group))]) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank()) + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(fill = "Group")

p2 = ggplot(dat,aes(Sample, Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))
p3=p1 / p2 + plot_layout(heights = c(1,10),guides = "collect" ) &
  theme(legend.position = "bottom")
pdf("immune_bar.pdf")
p3
dev.off()
cli_alert_success("生成免疫细胞柱状图成功：{dirname(args$data_dir)}/immune_bar.pdf")
cli_alert_info("生成免疫细胞箱线图。。。")
pdf("immune_box22.pdf")
ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot() + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))
dev.off()
cli_alert_success("生成22种免疫细胞箱线图成功：{dirname(args$data_dir)}/immune_box22.pdf")
# 全是0的行去掉
k = colSums(re)>0;table(k)
re = re[,k]
pdf("immune_box_diff.pdf")
draw_boxplot(t(re),factor(Group),
             drop = T,
             color = mypalette(length(unique(Group))))+
  labs(x = "Cell Type", y = "Estimated Proportion") 
dev.off()
cli_alert_success("生成组间差异显著的免疫细胞箱线图成功：{dirname(args$data_dir)}/immune_box22.pdf")
