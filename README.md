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
## 🚀 快速开始

### 单个GEO处理命令
```bash
Rscript geo_analysis.R \
  --GSEnumber GSE118370 \    # 必填参数
  --storage_dir ./results \  # 指定输出目录
  --GEO_variance_analysis \ # 启用差异分析
  --enrich_analysis         # 启用富集分析
```

#### 📌 参数详解
| 参数 | 缩写 | 类型 | 默认值 | 功能 |
|------|------|------|--------|------|
| `--GSEnumber` | `-n` | string | 必填 | GEO 数据集编号 (如 GSE118370) |
| `--storage_dir` | `-d` | path | 必填 | GEO矩阵输入和结果存储目录 |
| `--GEO_variance_analysis` | `-a` | flag | FALSE | 启用差异表达分析 |
| `--enrich_analysis` | `-e` | flag | FALSE | 启用富集分析 (需配合 -a) |

### 多个GEO处理命令（高级流程）
```bash
Rscript geo_analysis.R \
  --GSEnumber GSE118370 \    # 必填参数
  --storage_dir ./results \  # 指定输出目录
  --GEO_variance_analysis \ # 启用差异分析
  --enrich_analysis         # 启用富集分析
```

### 🧪 单个GEO使用示例
#### 数据准备
cd ./your_path
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118370/matrix/GSE118370_series_matrix.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118370/soft/GSE118370_family.soft.gz
##### 如果使用此数据进行科学研究，请引用：
Xu L, Lu C, Huang Y, Zhou J et al. SPINK1 promotes cell growth and metastasis of lung adenocarcinoma and acts as a novel prognostic biomarker. BMB Rep 2018 Dec;51(12):648-653.
#### 📂 输入文件结构（./path_of_folder_including_GSE118370_files）
```
./your_path/
├── GSE118370_series_matrix.txt.gz
├── GSE118370_family.soft.gz
```
#### 案例 1：数据处理流程
```bash
# 预处理
Rscript geo_analysis.R -n GSE118370 -d ./path_of_folder_including_GSE118370_files
```

#### 案例 2：完整分析流程
```bash
# 预处理+差异分析+富集分析
Rscript geo_analysis.R -n GSE118370 -d ./path_of_folder_including_GSE118370_files -a -e
```

## 📂 输出文件结构
```
./path_of_folder_including_GSE118370_files
├── GSE118370_matrix.csv      # 标准化表达矩阵（log2）
├── GSE118370_group.xls       # 样本分组信息
├── GSE118370_DiffEG.xls      # 差异基因表 (需 -a 参数)
├── GSE118370_PCA.pdf         # PCA 可视化
├── volcano_plot.pdf          # 火山图
├── GO_enrichment/           # 富集分析结果 (需 -e 参数)
│   ├── GSE118370_GO_barplot.pdf
│   └── GSE118370_GO_dotplot.pdf
│   └── GSE118370_GO.xls
└── KEGG_enrichment/
    ├── GSE118370_KEGG_barplot.pdf
    └── GSE118370_KEGG_dotplot.pdf
    └── GSE118370_KEGG.xls
```

### 🧪 多个GEO使用示例
#### 数据准备
#### 📂 输入文件结构（）
#### 案例 1：分析流程
```bash


```
## 📂 输出文件结构
```

```

## ⚠️ 注意事项
1. **富集分析需要在联网条件下进行**：
2. **默认情况下，富集分析和PPI分析参考基因和蛋白基于人类物种**


## 📜 许可证
本项目采用 [GPL-3.0 license](LICENSE)，商业使用需额外授权


