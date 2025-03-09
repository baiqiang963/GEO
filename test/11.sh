###手动分组
Rscript /home/baiqiang/GEO_STUDY/mul/geo_group_analysis.R -n GSE31821 -d /home/baiqiang/GEO_STUDY/mul
Rscript /home/baiqiang/GEO_STUDY/mul/geo_group_analysis.R -n GSE41177 -d /home/baiqiang/GEO_STUDY/mul
Rscript /home/baiqiang/GEO_STUDY/mul/geo_group_analysis.R -n GSE79768 -d /home/baiqiang/GEO_STUDY/mul
Rscript /home/baiqiang/GEO_STUDY/mul/geo_group_analysis2.R -i /home/baiqiang/GEO_STUDY/mul/GSE31821_group.xls -c characteristics_ch1 -f "patient disease status: atrial fibrillation patient" --overwrite
Rscript /home/baiqiang/GEO_STUDY/mul/geo_group_analysis2.R -i /home/baiqiang/GEO_STUDY/mul/GSE79768_group.xls -c characteristics_ch1.3 -f "condition: Atrial Fibrillation" --overwrite
Rscript /home/baiqiang/GEO_STUDY/mul/geo_group_analysis2.R -i /home/baiqiang/GEO_STUDY/mul/GSE41177_group.xls -c title -f "AF_*" --overwrite
###预处理（id转symbol->log2计数检测并转换->矩阵合并）
sh /home/baiqiang/GEO_STUDY/mul/GEO_multiple.sh -d /home/baiqiang/GEO_STUDY/mul --logfc 0.5 --fdr 0.05 --sample_cut_height 80 GSE41177 GSE31821 GSE79768