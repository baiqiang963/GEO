# 列出所有需要检查的包（去重后）
required_packages <- unique(c(
  "optparse", "GEOquery", "dplyr", "cli", "readr", "stringr", "data.table", 
  "tidyr", "tidyverse", "janitor", "limma", "affy", "ggplot2", "pheatmap", 
  "sva", "RColorBrewer", "glmnet", "randomForest", "caret", "cowplot", 
  "ggplotify", "VennDiagram", "WGCNA", "clusterProfiler", "ggthemes", 
  "org.Hs.eg.db", "enrichplot", "STRINGdb", "igraph", "ggraph","CIBERSORT"
))

# 初始化一个空向量来存储未安装的包
not_installed <- c()

# 遍历所有包，检查是否已安装
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    not_installed <- c(not_installed, pkg)
  }
}

# 输出未安装的包的名称
if (length(not_installed) > 0) {
  cat("以下包未安装:\n")
  cat(paste(not_installed, collapse = "\n"), "\n")
} else {
  cat("所有包均已安装。\n")
}
