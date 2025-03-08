#!/bin/bash

# R脚本路径（默认与Shell脚本同目录）
R_SCRIPT="$(dirname "$0")/GEO_detect.R"

# 使用帮助说明
usage() {
    echo "使用方法: $0 -d 数据目录 [--logfc 阈值] [--fdr 阈值] [--sample_cut_height 高度] GSE列表..."
    echo "示例：$0 -d /path/to/data --logfc 1 --fdr 0.05 --sample_cut_height 120 GSE123 GSE456"
    exit 1
}

# 检查R脚本是否存在
if [[ ! -f "$R_SCRIPT" ]]; then
    echo "错误：未找到R脚本 GEO_detect.R"
    echo "请确保 GEO_detect.R 与脚本 $0 位于同一目录下"
    exit 1
fi

# 默认值
LOGFC="1"
FDR="0.05"
SAMPLE_CUT_HEIGHT=""

# 解析参数
DATA_HOME=""  # 数据目录，必须手动指定
while [[ $# -gt 0 ]]; do
    case "$1" in
        -d)
            DATA_HOME="$2"
            shift 2
            ;;
        --logfc)
            LOGFC="$2"
            shift 2
            ;;
        --fdr)
            FDR="$2"
            shift 2
            ;;
        --sample_cut_height)
            SAMPLE_CUT_HEIGHT="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        -*)
            echo "错误：未知选项 $1"
            usage
            ;;
        *)
            break
            ;;
    esac
done

# 检查数据目录是否指定
if [[ -z "$DATA_HOME" ]]; then
    echo "错误：必须指定数据目录"
    usage
fi

# 检查数据目录是否存在
if [[ ! -d "$DATA_HOME" ]]; then
    echo "错误：数据目录不存在: $DATA_HOME"
    exit 1
fi

# 检查是否指定了GSE编号
if [[ $# -eq 0 ]]; then
    echo "错误：需要指定GSE编号"
    usage
fi

# 定义关键文件路径
CORRECTED_FILE="${DATA_HOME}/batch_corrected_matrix.csv"
TRAIT_DATA="${DATA_HOME}/merged_trait_data.csv"

# 检查是否跳过数据处理
SKIP_PROCESSING=false
if [[ -f "$CORRECTED_FILE" && -f "$TRAIT_DATA" ]]; then
    echo "检测到已存在的校正矩阵和表型数据，跳过数据处理步骤。"
    SKIP_PROCESSING=true
fi

if [ "$SKIP_PROCESSING" = false ]; then
    # 预处理主循环
    for GSE_ID in "$@"; do
        matrix_file="${DATA_HOME}/${GSE_ID}_series_matrix.txt.gz"
        soft_file="${DATA_HOME}/${GSE_ID}_family.soft.gz"
        
        # 严格模式检查文件
        check_failed=0
        [ ! -f "$matrix_file" ] && echo "错误：缺失矩阵文件 $matrix_file" && check_failed=1
        [ ! -f "$soft_file" ] && echo "错误：缺失注释文件 $soft_file" && check_failed=1
        [ $check_failed -eq 1 ] && continue

        # 执行预处理并记录耗时
        echo "▨▨▶ 预处理启动 $GSE_ID [$(date +%T)]"
        start_time=$SECONDS
        
        if Rscript "$R_SCRIPT" -n "$GSE_ID" -d "$DATA_HOME"; then
            elapsed=$((SECONDS - start_time))
            echo "◉ 预处理完成 [耗时: ${elapsed}s]"
        else
            echo "⊗ 预处理失败，检查日志"
        fi
        echo "───────────────────────────────────"
    done

    # 收集成功处理的GSE列表
    SUCCESS_GSE=()
    for GSE_ID in "$@"; do
        grouped_file="${DATA_HOME}/${GSE_ID}_group.xls"
        matrix_file="${DATA_HOME}/${GSE_ID}_matrix.csv"
        
        if [[ -f "$grouped_file" && -f "$matrix_file" ]]; then
            SUCCESS_GSE+=("$GSE_ID")
        else
            echo "警告：$GSE_ID 预处理结果不完整，跳过合并"
        fi
    done

    # 执行矩阵合并
    if [[ ${#SUCCESS_GSE[@]} -gt 0 ]]; then
        echo "▨▨▨▶ 开始合并矩阵 [总数: ${#SUCCESS_GSE[@]}]"
        GSE_LIST=$(IFS=','; echo "${SUCCESS_GSE[*]}")
        MERGE_SCRIPT="$(dirname "$0")/geo_matrix_merger.R"
        if [[ -f "$MERGE_SCRIPT" ]]; then
            Rscript "$MERGE_SCRIPT" \
                -g "$GSE_LIST" \
                -d "$DATA_HOME" \
                -o "${DATA_HOME}/merged_matrix_$(date +%Y%m%d).csv"
                
            merged_file="${DATA_HOME}/merged_matrix_$(date +%Y%m%d).csv"
            if [[ -f "$merged_file" ]]; then
                echo "◉ 合并完成，文件大小: $(du -h $merged_file | cut -f1)"
            else
                echo "⊗ 合并失败，请检查R脚本输出"
            fi
        else
            echo "错误：未找到合并脚本 geo_matrix_merger.R"
        fi
    else
        echo "警告：没有可合并的有效数据集"
    fi

    # 执行批次效应校正
    CORRECTION_SCRIPT="$(dirname "$0")/geo_batch_correction.R"
    if [[ -f "$CORRECTION_SCRIPT" && -f "$merged_file" ]]; then
        echo "▨▨▨▶ 开始批次效应矫正"
        GROUP_FILES=()
        for GSE_ID in "${SUCCESS_GSE[@]}"; do
            GROUP_FILE="${DATA_HOME}/${GSE_ID}_group.xls"
            [[ -f "$GROUP_FILE" ]] && GROUP_FILES+=("$GROUP_FILE")
        done
        
        if [[ ${#GROUP_FILES[@]} -ge 1 ]]; then
            GROUP_LIST=$(IFS=','; echo "${GROUP_FILES[*]}")
            Rscript "$CORRECTION_SCRIPT" \
                 -i "$merged_file" \
                -g "$GROUP_LIST" \
                -o "$CORRECTED_FILE" \
                -p "${DATA_HOME}/batch_correction_qc"
            
            if [[ -f "$CORRECTED_FILE" ]]; then
                echo "◉ 批次矫正完成，输出文件: $CORRECTED_FILE"
            else
                echo "⊗ 批次矫正失败"
            fi
        else
            echo "⊗ 没有有效的分组文件，跳过批次矫正"
        fi
    fi
    # 合并表型数据
    if [[ ${#SUCCESS_GSE[@]} -gt 0 ]]; then
        echo "▨▨▨▶ 合并表型数据"
        GSE_LIST=$(IFS=','; echo "${SUCCESS_GSE[*]}")
        MERGE_PHENO_SCRIPT="$(dirname "$0")/merge_pheno_data.R"
        Rscript "$MERGE_PHENO_SCRIPT" \
            -i "$DATA_HOME" \
            -o "$TRAIT_DATA" \
            -g "$GSE_LIST" \
            --id_col "sample_id" \
            --group_col "group"
        
        if [[ -f "$TRAIT_DATA" ]]; then
            echo "◉ 表型数据合并完成: $(wc -l < "$TRAIT_DATA") 行"
        else
            echo "⊗ 表型数据合并失败"
        fi
    fi
fi
# 执行差异分析
DIFF_SCRIPT="$(dirname "$0")/limma_analysis.R"
if [[ -f "$DIFF_SCRIPT" && -f "$CORRECTED_FILE" ]]; then
    echo "▨▨▨▶ 开始差异分析"
    Rscript "$DIFF_SCRIPT" \
        -m "$CORRECTED_FILE" \
        -g "${DATA_HOME}/*_group.xls" \
        -o "$DATA_HOME" \
        -f "$FDR" \
        -l "$LOGFC"
        
    if [[ -f "${DATA_HOME}/differential_genes_significant.csv" ]]; then
        echo "◉ 差异分析完成，结果目录: $DATA_HOME"
    else
        echo "⊗ 差异分析失败"
    fi
fi
# 执行WGCNA分析
WGCNA_SCRIPT="$(dirname "$0")/geo_wgcna_analysis.R"
if [[ -f "$WGCNA_SCRIPT" && -f "$CORRECTED_FILE" && -f "$TRAIT_DATA" ]]; then
    echo "▨▨▨▶ 开始WGCNA分析（高变异基因筛选）"
    WGCNA_PARAMS=()
    [[ -n "$SAMPLE_CUT_HEIGHT" ]] && WGCNA_PARAMS+=("--sample_cut_height" "$SAMPLE_CUT_HEIGHT")
    
    Rscript "$WGCNA_SCRIPT" \
        -e "$CORRECTED_FILE" \
        -t "$TRAIT_DATA" \
        -o "$DATA_HOME" \
        "${WGCNA_PARAMS[@]}" \
    
    if [[ -f "${DATA_HOME}/gene_module_membership.csv" ]]; then
        echo "◉ WGCNA分析完成，结果目录: $DATA_HOME"
    else
        echo "⊗ WGCNA分析失败"
    fi
else
    echo "警告：缺少必要文件，跳过WGCNA分析"
fi
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
echo "▨▨▨ 所有分析流程执行完毕 ▨▨▨"
