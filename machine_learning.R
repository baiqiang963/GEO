suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(cli)
  library(RColorBrewer)
  library(tidyverse)
  library(glmnet)
  library(randomForest)
  library(caret)
  library(ggplot2)
  library(cowplot)
  library(ggplotify)
  library(VennDiagram)
})
option_list <- list(
  make_option(c("-d", "--data_dir"), type = "character",
              help = "批次矫正后的矩阵文件存储路径 [batch_corrected_matrix.csv]"),
  make_option(c("-t", "--trait_data"), type = "character",
              help = "合并后的临床信息路径 [merged_trait_data.csv]"),
  make_option(c("-l", "--hub_data_limma"), 
              type = "character",
              help = "来自limma分析的hub基因矩阵结果文件 [differential_genes_significant.csv]"),
  make_option(c("-w", "--hub_data_wgcna"), 
              type = "character",
              help = "来自WGCNA分析的hub基因矩阵结果文件 [WGCNA_ip_gene.csv]"),
  make_option(c("-o", "--out_dir"), type = "character",help = "输出路径")
)
parser <- OptionParser(usage = "%prog -d batch_corrected_matrix.csv -t merged_trait_data.csv -o 输出文件夹路径", 
                       option_list = option_list)
args <- parse_args(parser)

if(is.null(args$data_dir)){
  stop("未找到矩阵文件，检查输入")
}
if(is.null(args$trait_data)){
  stop("未找到临床信息/表型文件，检查输入")
}
if(is.null(args$out_dir)){
  stop("未定义输出路径")
}

setwd(args$out_dir)
cli_h1("机器学习筛选基因的分析模式判断")
if (!is.null(args$hub_data_limma) && !is.null(args$hub_data_wgcna)) {
  cli_alert("发现hub_gene输入来自limma和WGCNA，两者取交集操作后，执行机器学习筛选")
  limma=fread(args$hub_data_limma)[,"gene"]
  cli_alert_success("读取limma差异基因数量：{nrow(limma)}")
  wgcna=fread(args$hub_data_wgcna)[,"gene"]
  cli_alert_success("读取wgcna差异基因数量：{nrow(wgcna)}")
  hub_gene=merge(limma,wgcna,by="gene")
  cli_alert_success("取交集成功，保留基因：{nrow(hub_gene)}")
  fwrite(hub_gene,"limma_wgcna_merge_gene.csv",quote = F,row.names = F)
  cli_alert_success("Limma与wgcna交集结果：{args$out_dir}/limma_wgcna_merge_gene.csv")
  #pdf("limma_wgcna_vn.pdf")
  venn.diagram(
    x = list(limma = limma$gene, WGCNA = wgcna$gene),
    filename="limma_wgcna_vn.pdf",
    fill = c("lightblue", "lightgreen"), # 设置颜色
    alpha = 0.5, # 设置透明度
    cat.cex = 1.2, # 设置类别标签大小
    cex = 1.5, # 设置数字大小
    cat.pos = c(-30, 30), # 设置类别标签位置
    cat.dist = c(0.05, 0.05), # 设置类别标签距离
    cat.fontface = "bold", # 设置类别标签字体
    fontface = "bold", # 设置数字字体
    margin = 0.1 # 设置边距
  )
  #dev.off()
  cli_alert_success("Limma与wgcna交集韦恩图生成：{args$out_dir}/limma_wgcna_vn.pdf")
} else if (!is.null(args$hub_data_limma)) {
  cli_alert("发现hub_gene输入来自limma，执行机器学习筛选")
  hub_gene=fread(args$hub_data_limma)[,"gene"]
  cli_alert_success("读取limma差异基因数量：{nrow(hub_gene)}")
} else if (!is.null(args$hub_data_wgcna)) {
  cli_alert("发现hub_gene输入来自WGCNA，执行机器学习筛选")
  hub_gene=fread(args$hub_data_wgcna)[,"gene"]
  cli_alert_success("读取wgcna差异基因数量：{nrow(hub_gene)}")
} else {
  stop("请检查 -l 和 -w，至少输入一个hub gene矩阵")
}
colnames(hub_gene) <- "symbol"
train_data <- read.csv(args$data_dir, row.names = 1, check.names = F)  # 行名为全部基因名，每列为样本名
cli_alert("训练集：{args$data_dir}")
group <- read.csv(args$trait_data)
cli_alert("分组信息：{args$trait_data}")
dat <- train_data[rownames(train_data) %in% hub_gene$symbol, ] %>%
  t() %>%
  as.data.frame() # 整理后行为样本名，列为基因名
dat$sample_id <- rownames(dat)
dat <- merge(dat, group, var = "sample_id") # 筛选后的表达矩阵与分组信息表合并
dat <- column_to_rownames(dat, var = "sample_id") %>% as.data.frame() # 列名转为行名

table(dat$group)
dat$group <- factor(dat$group, levels = c('case', 'control'))
cli_h1("执行LASSO")
cc=381
set.seed(cc) # 设置种子
cli_alert_info("种子:{cc}")
res.lasso <- cv.glmnet(x = as.matrix(dat[-ncol(dat)]), 
                       y = dat$group, 
                       family = "binomial", 
                       type.measure = "default",
                       nfolds = 10)
s="lambda.min"
coef.min <- coef(res.lasso, s = s)
cli_alert_info("lambda.min——是交叉验证过程中使得平均偏差（或其他指定的损失函数）最小的lambda值。换
句话说，它是根据交叉验证结果，直接选择出的最优lambda值，因为它在平均意义上给出了最好的预测性能")
cli_alert_info("lambda.1se——lambda.1se是lambda.min之后，使得平均偏差增加不超过一倍标准误差的最大的lambda值。
这个选择标准旨在在模型的预测性能和复杂度之间找到一个折中。通过选择一个稍微更大的lambda值（即更强的正则化），l
ambda.1se可以帮助避免过拟合，在接受一些预测性能的损失的同时来换取更好的稳健性")
cli_alert_info("目前使用的LASSO选择标准：{s}")

# 找出那些回归系数没有被惩罚为0的
active.min <- which(coef.min@i != 0)
# 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i + 1]
lasso_geneids <- lasso_geneids[-1] %>% as.data.frame()
colnames(lasso_geneids) <- 'gene'

fwrite(lasso_geneids, file = 'lasso_gene.csv',row.names = F,quote = F)
cli_alert_success("LASSO筛选成功，hub基因结果保存至：{args$out_dir}/lasso_gene.csv")
pdf("LASSO.pdf",width = 15,height = 8)
par(mfrow = c(1, 2))
plot(res.lasso$glmnet.fit, 
     xvar = 'lambda', 
     label = TRUE, 
     las = 1)
abline(v=log(c(res.lasso$lambda.min, res.lasso$lambda.1se)), lty="dashed")
plot(res.lasso, 
     las =1)
dev.off()
cli_alert_success("LASSO结果可视化：{args$out_dir}/LASSO.pdf")
cli_alert_success("LASSO筛选完成")
#-----------------------------------------------------------------------
cli_h1("执行SVM-RFE")
aa=21
set.seed(aa) # 设置种子
cli_alert_info("种子:{aa}")
number = 10
control <- rfeControl(functions = caretFuncs, method = "cv", number = number) 
cli_alert_info("交叉验证次数：{number}")
# 执行SVM-RFE算法
num <- ncol(dat)-1
results <- rfe(x = dat[, 1:num],
               y = dat$group, 
               sizes = c(1:num), 
               rfeControl = control,
               method = "svmRadial"
)

## 结果分析
svmrfe_result <- data.frame(gene = predictors(results))

fwrite(svmrfe_result, file = 'svm_rfe_gene.csv',row.names = F,quote = F)
cli_alert_success("SVM-RFE筛选成功，hub基因结果保存至：{args$out_dir}/svm_rfe_gene.csv")
# SVM-RFE结果简单可视化
p1 <- plot(results, type=c("o"),
           xgap.axis = 1)
p1 <- as.ggplot(plot_grid(p1))+
  labs(title="SVM_RFE_analyse", x="", y = "",size=25) +
  # theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=25),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pdf("SVM-RFE.pdf")
p1
dev.off()
cli_alert_success("SVM-RFE结果可视化：{args$out_dir}/SVM-RFE.pdf")
cli_alert_success("SVM-RFE筛选完成")
#-----------------------------------------------------------------------
cli_h1("执行随机森林")
x<-dat[,!colnames(dat)%in%"group"]
y<-factor(dat$group)
rf<-randomForest(x, y,data=dat,ntree=1000)
pdf("random_forest.pdf", width = 6, height = 6)
plot(rf, main = "Random forest", lwd = 2)
dev.off()
# 找出误差最小的点
optionTrees <- which.min(rf$err.rate[, 1])
rf2 <- randomForest(x, y,data=dat,ntree= optionTrees)
# 基因重要性
importance <- importance(x = rf2)
# 绘图
pdf("geneImportance.pdf", width = 6.2, height = 6)
varImpPlot(rf2, main = "")
dev.off()
cli_alert_success("随机森林结果可视化：{args$out_dir}/geneImportance.pdf")
rfGenes <- importance[order(importance[, "MeanDecreaseGini"], decreasing = TRUE), ]
# 重要性评分大于2的基因
rfGenes <- names(rfGenes[rfGenes>2])
fwrite(data.frame(gene = rfGenes), file = "randomforest_Gene.csv", quote = F, row.names = F)
cli_alert_success("随机森林筛选成功，hub基因结果保存至：{args$out_dir}/randomforest_Gene.csv")
cli_h1("所有机器学习筛选基因执行完成")
