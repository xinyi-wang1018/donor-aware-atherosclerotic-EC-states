GSE247238 Step 4：参考状态映射与患者层面pseudobulk验证
============================================================

研究对象：第二独立队列 GSE247238。
参考来源：53样本主队列Seurat对象。

核心设计：
1. 从参考对象自动提取 C7、C3、C5、C1/C2、独立C4 五类特征。
2. C1与C2合并为晚期重塑区；C4作为独立拒绝/旁支类。
3. 九个候选基因和SULF1不参与分类器构建，避免循环验证。
4. 细胞映射只用于状态归属；正式统计以患者pseudobulk为单位。
5. primary pseudobulk仅纳入 high/moderate 映射细胞。

最重要文件：
- /root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk/02_cell_mapping/05_cell_level_state_mapping.csv.gz
  每个EC细胞的映射状态、相关性、margin和置信度。
- /root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk/02_cell_mapping/08_state_counts_proportions_by_patient_condition.csv
  每位患者Control/Plaque中的状态数量和比例。
- /root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk/02_cell_mapping/09_state_proportion_paired_statistics.csv
  CAR2-CAR10状态比例的患者配对统计。
- /root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk/03_patient_pseudobulk/07_patient_state_candidate_SULF1_logCPM.csv
  patient×state候选基因和SULF1的pseudobulk表达。
- /root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk/03_patient_pseudobulk/11_candidate_SULF1_state_contrast_edgeR.csv
  C3/C5/晚期状态之间的患者配对edgeR结果。
- /root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk/03_patient_pseudobulk/13_candidate_SULF1_Plaque_vs_Control_within_state_edgeR.csv
  同一映射状态内Plaque-Control患者配对edgeR结果。
- /root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk/03_patient_pseudobulk/15_candidate_vs_SULF1_patient_state_precedence.csv
  每位患者中候选状态顺序是否早于SULF1。
- /root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk/03_patient_pseudobulk/16_candidate_preSULF1_state_order_summary.csv
  九候选的最终“早于SULF1状态顺序支持”汇总。
- /root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk/03_patient_pseudobulk/README_patient_pseudobulk_and_preSULF1.txt
  自动摘要和解释边界。

图片：
- /root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk/plots/01_UMAP_mapped_states.png
- /root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk/plots/02_state_composition_by_sample.png
- /root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk/plots/03_candidate_SULF1_state_heatmap.png
- /root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk/plots/04_mapping_confidence_margin.png
- /root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk/plots/05_candidate_SULF1_patient_state_profiles.png

重要解释：
- “早于SULF1”仅表示在映射状态顺序中更早出现，不等于真实纵向时间或因果关系。
- C4不作为C3后的第二条动脉分支。
- low_reject细胞不进入主要pseudobulk分析。
