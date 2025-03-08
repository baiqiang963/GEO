# GEO数据分析管道
这是一套可自定义参数和步骤的高度自动化的对多个GEO数据进行合并分析的管道流程。
![](https://github.com/baiqiang963/GEO/blob/main/images/pipeline.png)
## 第一次运行前的环境准备
### R>=4.2.0
### 将PPI分析所需的蛋白参考文件下载到脚本所在目录中(9606是人类的蛋白参考文件代号)
```
#bash
cd your_code_dir_path
wget https://stringdb-static.org/download/protein.aliases.v11.5/9606.protein.aliases.v11.5.txt.gz
wget https://stringdb-static.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz
wget https://stringdb-static.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz
```
### 依赖R包缺失检查
```
#bash
module load R-4.2.0
Rscript check.R
```
## 📋 脚本功能
**1.GEO数据预处理（探针ID转Gene symbol，自动检测并修正log2计数，多个GEO数据合并并生成对应分组信息）**
**2.limma差异分析和WGCNA分析筛选差异基因**
**3.机器学习筛选特征：LASSO、SVM、随机森林**
**4.GO和KEGG分析**
**5.PPI分析**
**6.免疫浸润分析**
## 🚀 快速开始示例
**step0：使用以下脚本做好每个GEO数据集的性状分组case/control**
调动你的R语言
```
#bash
module load R-4.2.0
```
**step1：使用以下脚本做好每个GEO数据集的性状分组case/control**
geo_group_analysis.R -n 数据集的GSE编号 -d matrix和soft文件所在的目录(要求为同一目录)
geo_group_analysis2.R -i 由geo_group_analysis.R生成潜在分组指标表格 -c 指示性状的列/可以被用来分组的列的列名 -f 模糊匹配指示性状的列/可以被用来分组的列中case性状的关键词
```
#bash
module load R-4.2.0
Rscript geo_group_analysis.R -n GSE31821 -d /home/baiqiang/GEO_STUDY/mul
Rscript geo_group_analysis.R -n GSE41177 -d /home/baiqiang/GEO_STUDY/mul
Rscript eo_group_analysis.R -n GSE79768 -d /home/baiqiang/GEO_STUDY/mul
Rscript geo_group_analysis2.R -i /home/baiqiang/GEO_STUDY/mul/GSE31821_group.xls -c characteristics_ch1 -f "patient disease status: atrial fibrillation patient" --overwrite
Rscript geo_group_analysis2.R -i /home/baiqiang/GEO_STUDY/mul/GSE79768_group.xls -c characteristics_ch1.3 -f "condition: Atrial Fibrillation" --overwrite
Rscript geo_group_analysis2.R -i /home/baiqiang/GEO_STUDY/mul/GSE41177_group.xls -c title -f "AF_*" --overwrite
```
**step2：geo_pipeline分析(第一次运行，未剔除WGCNA分析中的离群样本)**
这次运行会对三个数据集进行预处理以及合并处理，生成合并后表型分组数据。随后执行limma差异分析和wgcna分析。
但是无法提前得知样本剪切高度，需要等这次流程中WGCNA分析执行完毕以后。查看Sample_clustering_to_detect_outliers.pdf确认是否存在离群样本。
在Sample_clustering_to_detect_outliers.pdf生成并发现离群样本后，可以提前终止流程继续往下执行，从而提高效率。
如果未发现离群样本，step2可以执行到底，无需进行step3.
--logfc 0.5 筛选标准：limma分析的差异基因logfc绝对值 默认1
--fdr 0.05 筛选标准：limma分析的差异基因fdr值 默认0.05
```
#bash
sh GEO_multiple.sh -d /home/baiqiang/GEO_STUDY/mul --logfc 0.5 --fdr 0.05 GSE41177 GSE31821 GSE79768
```
**step3：geo_pipeline分析(第二次运行，剔除WGCNA分析中的离群样本)**
第二次执行时，会自动沿用第一次执行时的预处理结果，从limma分析开始跑。
在第一次运行中Sample_clustering_to_detect_outliers.pdf中发现离群样本后，在第二次运行时使用--sample_cut_height 设置剪切高度进行进行离群样本剔除。
```
#bash
module load R-4.2.0
sh GEO_multiple.sh -d /home/baiqiang/GEO_STUDY/mul --logfc 0.5 --fdr 0.05 --sample_cut_height 80 GSE41177 GSE31821 GSE79768
```

##自定义分析流程
**修改GEO_multiple.sh中的参数可以改变分析流程：**
*默认情况下，机器学习筛选特征是针对limma分析和WGCNA分析的交集基因来筛选特征。*
```
# 定义机器学习脚本路径
ML_SCRIPT="$(dirname "$0")/machine_learning.R"

# 检查脚本是否存在
if [[ -f "$ML_SCRIPT" ]]; then
  echo "发现机器学习脚本: $ML_SCRIPT，开始执行..."
  
  # 执行机器学习
  Rscript "$ML_SCRIPT" \
          -d "$CORRECTED_FILE" \
          -t "$TRAIT_DATA" \
          -l "${DATA_HOME}/differential_genes_significant.csv" \
          -w "${DATA_HOME}/WGCNA_ip_gene.csv" \
          -o "$DATA_HOME"
  
  # 检查执行结果
  if [[ $? -eq 0 ]]; then
    echo "机器学习执行成功"
  else
    echo "警告：机器学习执行失败，请检查日志"
  fi
else
  echo "警告：未找到机器学习脚本: $ML_SCRIPT，跳过执行"
fi
```
如果针对limma分析或WGCNA分析的基因进行筛选，可以仅保留-l和-w中的一个。
```
# 定义机器学习脚本路径
ML_SCRIPT="$(dirname "$0")/machine_learning.R"

# 检查脚本是否存在
if [[ -f "$ML_SCRIPT" ]]; then
  echo "发现机器学习脚本: $ML_SCRIPT，开始执行..."
  
  # 执行机器学习
  Rscript "$ML_SCRIPT" \
          -d "$CORRECTED_FILE" \
          -t "$TRAIT_DATA" \
          -l "${DATA_HOME}/differential_genes_significant.csv" \
		  #-w "${DATA_HOME}/WGCNA_ip_gene.csv" \
          -o "$DATA_HOME"
  
  # 检查执行结果
  if [[ $? -eq 0 ]]; then
    echo "机器学习执行成功"
  else
    echo "警告：机器学习执行失败，请检查日志"
  fi
else
  echo "警告：未找到机器学习脚本: $ML_SCRIPT，跳过执行"
fi
```
去除-w参数后，那么机器学习仅针对limma分析的结果执行。
*同理：默认情况下，富集分析是针对limma分析和WGCNA分析的交集基因进行的。*
```
# 定义富集分析脚本路径
ENRICH_SCRIPT="$(dirname "$0")/geo_enrich_analysis.R"

# 检查脚本是否存在
if [[ -f "$ENRICH_SCRIPT" ]]; then
  echo "发现富集分析脚本: $ENRICH_SCRIPT，开始执行..."
  
  # 执行富集分析
  Rscript "$ENRICH_SCRIPT" \
          -d "${DATA_HOME}/limma_wgcna_merge_gene.csv" \
  
  # 检查执行结果
  if [[ $? -eq 0 ]]; then
    echo "富集分析执行成功"
  else
    echo "警告：富集分析执行失败，请检查日志"
  fi
else
  echo "警告：未找到富集分析脚本: $ENRICH_SCRIPT，跳过执行"
fi
```
可以改-d 使其针对limma分析或WGCNA分析的结果进行富集分析
```
#针对limma分析
-d "${DATA_HOME}/differential_genes_significant.csv"
#针对WGCNA分析
-d "${DATA_HOME}/WGCNA_ip_gene.csv"
```
*同理：默认情况下，PPI是针对limma分析和WGCNA分析的交集基因进行的。*
```
# 定义PPI分析脚本路径
PPI_SCRIPT="$(dirname "$0")/PPI.R"

# 检查脚本是否存在
if [[ -f "$PPI_SCRIPT" ]]; then
  echo "发现富集分析脚本: $PPI_SCRIPT，开始执行..."
  
  # 执行PPI分析
  Rscript "$PPI_SCRIPT" \
          -d "${DATA_HOME}/limma_wgcna_merge_gene.csv" \
          -s "$(dirname "$0")" \
          -t 400
  # 检查执行结果
  if [[ $? -eq 0 ]]; then
    echo "PPI分析执行成功"
  else
    echo "警告：PPI分析执行失败，请检查日志"
  fi
else
  echo "警告：未找到富集分析脚本: $PPI_SCRIPT，跳过执行"
fi
```
可以改-d 使其针对limma分析或WGCNA分析的结果进行PPI分析
```
#针对limma分析
-d "${DATA_HOME}/differential_genes_significant.csv"
#针对WGCNA分析
-d "${DATA_HOME}/WGCNA_ip_gene.csv"
```
## ⚠️ 注意事项
1. **富集分析需要在联网条件下进行**：
2. **默认情况下，富集分析和PPI分析参考基因和蛋白基于人类物种**
3. **limma分析和WGCNA分析的交集基因结果limma_wgcna_merge_gene.csv在机器学习脚本的开头代码中生成，如有需要可以手动取交集，生成一个命名为limma_wgcna_merge_gene.csv的无行名，列名为gene，列内容为gene symbol的table**

## 📜 许可证
本项目采用 [GPL-3.0 license](LICENSE)，商业使用需额外授权


