suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
  library(optparse)
  library(cli)
  library(RColorBrewer)
})
cli_h1("开启多线程运行")
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# 命令行参数配置 ---------------------------------------------------------------
option_list <- list(
  make_option(c("-e", "--expr_matrix"), type = "character",
              help = "表达矩阵文件路径（CSV格式，行名为基因，列名为样本）"),
  make_option(c("-t", "--trait_data"), type = "character",
              help = "表型数据文件路径（CSV格式，需包含sample_id和group列）"),
  make_option(c("-o", "--output_dir"), type = "character",
              default = "wgcna_results",
              help = "输出目录 [默认: wgcna_results]"),
  make_option("--var_threshold", type = "numeric", default = 0.25,
              help = "基因方差阈值 [默认: 0.25]"),
  make_option("--sample_cut_height", type = "numeric",
              help = "样本聚类切割高度"),
  make_option("--min_module_size", type = "integer", default = 200,
              help = "一步法网络构建时每个模块的最小模块基因数 [默认: 200]"),
  make_option("--network_type", type = "character", default = "signed",
              help = "网络类型：signed/unsigned [默认: signed]"),
  make_option("--MEDissThres", type = "numeric", default = 0.25,
              help = "切割高度0.25，即模块之间相关性达到(1-0.25) [默认: 0.25]")
)

parser <- OptionParser(usage = "%prog -e expr.csv -t trait.csv [options]", 
                       option_list = option_list)
args <- parse_args(parser)
setwd(args$output_dir)

# 自定义函数：自动选择软阈值 -----------------------------------------------------
auto_select_power <- function(sft, Rsquared_cut = 0.8, slope_cut = -1,
                              min_power = 5) {
  fit <- sft$fitIndices
  
  valid1 <- which(fit$SFT.R.sq >= Rsquared_cut & 
                    fit$slope < 0 & 
                    fit$Power >= min_power)
  if (length(valid1) > 0) return(fit$Power[valid1[1]])
  
  valid2 <- which(fit$SFT.R.sq >= (Rsquared_cut - 0.1) & 
                    fit$slope < 0 & 
                    fit$Power >= min_power)
  if (length(valid2) > 0) return(fit$Power[valid2[1]])
  
  slope_diff <- abs(fit$slope - slope_cut)
  candidate <- which(slope_diff <= 0.2 & fit$Power >= min_power)
  if (length(candidate) > 0) return(fit$Power[candidate[which.max(fit$SFT.R.sq[candidate])]])
  
  stop("无法自动选择Power值")
}
tryCatch({
fpkm <- read.csv(args$expr_matrix) 
rownames(fpkm) <- fpkm[,1]
fpkm<- fpkm[,-1]
WGCNA_matrix = t(fpkm[order(apply(fpkm,1,mad), decreasing = T)[1:round(nrow(fpkm)*0.3,0)],]) #按mad进行排序，取前30%的基因进行后续分析
datExpr0 = as.data.frame(WGCNA_matrix) 
datExpr_all = as.data.frame(t(fpkm))
cli_h1("数据质控")
gsg_all = goodSamplesGenes(datExpr_all, verbose = 3)
gsg_all$allOK
if (!gsg_all$allOK){
if (sum(!gsg_all$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr_all)[!gsg_all$goodGenes], collapse = ", ")));
if (sum(!gsg_all$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr_all)[!gsg_all$goodSamples], collapse = ", ")));
    datExpr_all = datExpr_all[gsg_all$goodSamples, gsg_all$goodGenes]
}
gsg_all = goodSamplesGenes(datExpr_all, verbose = 3)
cli_alert_success("goodSamplesGenes质控完成：{gsg_all$allOK}")
sampleTree = hclust(dist(datExpr0), method = "average")
if(is.null(args$sample_cut_height)){
pdf("Sample_clustering_to_detect_outliers.pdf")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
cli_warn("Sample_clustering质控完成：未设置sample_cut_height")
cli_warn("本次流程运行完成后，查看：{args$output_dir}/Sample_clustering_to_detect_outliers.pdf,确认是否存在需要剔除的离群值")
}
if(!is.null(args$sample_cut_height)){
  pdf("Sample_clustering_to_detect_outliers.pdf")
  par(cex = 0.6)
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  abline(h = args$sample_cut_height, col = "red")
  dev.off()
  clust = cutreeStatic(sampleTree, cutHeight = args$sample_cut_height, minSize = 10)
  table(clust)
  keepSamples = (clust==1)
  datExpr = datExpr0[keepSamples, ]
  cli_alert_success("Sample_clustering质控完成，已设置sample_cut_height：{args$sample_cut_height}")
  cli_alert_success("阈值线质控图已生成：{args$output_dir}/Sample_clustering_to_detect_outliers.pdf")
}
datTraits <- read.csv(args$trait_data,header=TRUE)
datTraits=datTraits[datTraits$sample_id%in%rownames(datExpr),]
rownames(datTraits) <- datTraits[,1]
datTraits<-datTraits[,c("group"),drop=F]
datTraits$group[datTraits$group=="case"]=1
datTraits$group[datTraits$group=="control"]=0
dim(datTraits)
traitColors = numbers2colors(as.numeric(as.matrix(datTraits)), signed = FALSE)
sampleTree = hclust(dist(datExpr), method = "average")
pdf("Sample_dendrogram_and_trait_heatmap_after_rm_outliers.pdf")
plotDendroAndColors(sampleTree, traitColors, groupLabels = names(datTraits),main = "Sample dendrogram and trait heatmap")
dev.off()
cli_alert_success("最终质控图已生成，请以此判断是否存在离群值：{args$output_dir}/Sample_dendrogram_and_trait_heatmap_after_rm_outliers.pdf")
save(datExpr, datTraits, file = "WGCNA0.3-dataInput.RData")

cli_h1("开始自动选择阈值")
powers = c(1:20) 
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
power=auto_select_power(sft)
fwrite(sft$fitIndices,"power_table.csv",quote = F,row.names = F)
cli_alert_success("自动选择阈值{power}")
# 设置PDF输出参数
pdf("soft-thresholding_power.pdf", width = 9, height = 5, paper = "special")

# 优化图形参数
par(
  mfrow = c(1, 2),  # 1行2列布局
  mar = c(4.5, 5, 2.5, 1),  # 边距调整：下、左、上、右
  mgp = c(3, 0.8, 0),       # 坐标轴标签位置
  cex.main = 1.2,           # 主标题缩放
  cex.lab = 1.1,            # 坐标轴标签缩放
  cex.axis = 0.9             # 坐标轴刻度缩放
)

# 绘制第一个图：无标度拓扑拟合
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     type = "n",  # 初始化画布不绘制点
     xlab = "Soft Threshold (power)",
     ylab = expression(paste("Scale Free Topology Fit (", R^2, ")")),
     main = "Scale Independence",
     ylim = c(0, 1))  # 固定y轴范围

# 添加半透明网格线
grid(col = "gray90", lty = "dotted")
points(sft$fitIndices[,1], 
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       pch = 21, 
       bg = scales::alpha("red", 0.6), 
       col = "darkred", 
       cex = 1.2)
text(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, 
     cex = 0.8, 
     pos = 3,  # 文字显示在点上侧
     col = "darkred")
abline(h = sft$fitIndices$SFT.R.sq[power], col = "blue", lwd = 2, lty = "dashed")

# 绘制第二个图：平均连接度
plot(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     type = "n",
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Network Connectivity",
     log = "y")  # y轴对数变换

grid(col = "gray90", lty = "dotted")
points(sft$fitIndices[,1], 
       sft$fitIndices[,5],
       pch = 22, 
       bg = scales::alpha("darkgreen", 0.6), 
       col = "darkgreen", 
       cex = 1.2)
text(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     labels = powers, 
     cex = 0.8, 
     pos = 3, 
     col = "darkgreen")

# 关闭图形设备
invisible(dev.off())
cli_h1("开始进行网络构建")
cli_text("deepSplit 参数调整划分模块的敏感度，值越大，越敏感，得到的模块就越多，默认是2")
cli_text("minModuleSize 参数设置最小模块的基因数，值越小，小的模块就会被保留下来")
cli_text("目前minModuleSize 参数：--min_module_size {args$min_module_size}")
cli_text("mergeCutHeight设置合并相似性模块的距离，值越小，就越不容易被合并，保留下来的模块就越多")
cli_text("目前mergeCutHeight参数：--MEDissThres {args$MEDissThres}")
net = blockwiseModules(
  datExpr,
  power = power,
  maxBlockSize = 6000,
  TOMType = "signed", minModuleSize = args$min_module_size,
  reassignThreshold = 0, mergeCutHeight = args$MEDissThres,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = F, 
  verbose = 3
)
table(net$colors) 

mergedColors = labels2colors(net$colors)
pdf("Cluster_Dendrogram.pdf")
plotDendroAndColors(net$dendrograms[[1]], 
                    mergedColors[net$blockGenes[[1]]], 
                    "Module colors", 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, file = "networkConstruction-auto.RData")
cli_alert_success("网络构建完成")
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
table(moduleColors)
pdf("Module-trait associations.pdf",width = 8, height=10)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = "trait", yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()
cli_alert_success("模块识别")
modNames = substring(names(MEs), 3)
traitNames = "case"
#计算 MM
MMcor= as.data.frame(cor(datExpr, MEs, use = "p"))
MMP = as.data.frame(corPvalueStudent(as.matrix(MMcor),
                                     nSamples))
names(MMcor) = paste("MM", modNames, sep="")
names(MMP) = paste("p.MM", modNames, sep="")
cli_alert_success("计算 MM完成")
#计算 GS
datcli=datTraits
GScor = as.data.frame(cor(datExpr,datcli, use = "p"))
GSP = as.data.frame(corPvalueStudent(as.matrix(GScor),
                                     nSamples))
names(GScor) = paste("GS.", traitNames, sep="");
names(GSP) = paste("p.GS.", traitNames, sep="");
MMGS<-cbind(MMcor, MMP,GScor, GSP)
write.table(MMGS, "MMGS", sep = "\t", quote = F)
write.table(MMGS, "MMGS.txt", sep = "\t", quote = F)
module=substring(rownames(moduleTraitPvalue)[which.min(moduleTraitPvalue[, "group"])], 3)
column = match(module, modNames)
moduleGenes = moduleColors==module
trait="case"
traitColumn=match(trait,traitNames)
cli_alert_success("计算 GS完成")
pdf("verboseScatterplot.pdf", width = 8, height = 6)
par(mfrow = c(1,1), mar = c(5,5,4,2))

# 计算参考线数值
h_value <- mean(abs(GScor[moduleGenes, 1]))
v_value <- mean(abs(MMcor[moduleGenes, column]))
verboseScatterplot(abs(MMcor[moduleGenes, column]),
                   abs(GScor[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", trait),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, 
                   cex.lab = 1.2, 
                   cex.axis = 1.2, 
                   col = module)
abline(h = h_value, v = v_value, col = "red", lty = 2)
# 添加数值标签
text(x = par("usr")[2]*0.95,  # 右侧95%位置
     y = h_value, 
     labels = sprintf("Avg GS = %.2f", h_value),
     col = "red", pos = 3, cex = 0.9, font = 2)

text(x = v_value, 
     y = par("usr")[4]*0.95,  # 顶部95%位置
     labels = sprintf("Avg MM = %.2f", v_value),
     col = "red", pos = 4, cex = 0.9, font = 2)

dev.off()
cli_alert_success("最重要的模块基因散点图已生成：{args$output_dir}/verboseScatterplot.pdf")
MMGS_module=MMGS
colnames(MMGS_module)=gsub("MM","",colnames(MMGS_module))
MMGS_module$gene=rownames(MMGS_module)
MMGS_module=MMGS_module[,c(module,"GS.case","gene")]
MMGS_module=MMGS_module[abs(MMGS_module$GS.case)>h_value,]
MMGS_module=MMGS_module[abs(MMGS_module[,1])>v_value,]
fwrite(MMGS_module,"WGCNA_ip_gene.csv",quote = F,row.names = F)
cli_alert_success("MM与GS平均线之上的重要基因已筛选完成：{nrow(MMGS_module)}个基因")
cli_alert_success("结果生成：{args$output_dir}/WGCNA_ip_gene.csv")
}, error = function(e) {
  cli_alert_danger("处理失败: {e$message}")
  quit(status = 1)
})
