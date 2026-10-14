GSE247238 Step 5：内皮身份与污染敏感性验证
================================================

目的：
1. 检查C3-like细胞是否保留内皮身份；
2. 排查LOX、CTHRC1等信号是否主要来自明显成纤维、平滑肌/周细胞或免疫共表达风险细胞；
3. 在三套逐步严格的细胞集合中重复患者级pseudobulk；
4. 确认ANGPT2、LOX、CTHRC1和SULF1的状态方向是否保留。

三套敏感性集合：
- S1_HM_ALL：high + moderate映射主轴EC；
- S2_HIGH：仅high映射主轴EC；
- S3_STRICT_EC：high + 两种映射方法一致 + 严格内皮身份QC。

严格过滤不使用九个候选或SULF1，因此不存在用目标基因定义过滤后再验证目标基因的循环分析。

最重要结果：
1. /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/01_cell_identity_QC/05_sensitivity_set_state_counts.csv
2. /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/01_cell_identity_QC/07_core_candidate_positive_cell_identity_QC.csv
3. /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/01_cell_identity_QC/08_C3_candidate_identity_score_correlations.csv
4. /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/03_statistics/01_state_contrast_QC_by_sensitivity_set.csv
5. /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/03_statistics/02_candidate_SULF1_sensitivity_edgeR_long.csv
6. /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/03_statistics/04_patient_direction_summary_by_sensitivity_set.csv
7. /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/03_statistics/05_candidate_robustness_summary.csv
8. /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/03_statistics/06_core_ANGPT2_LOX_CTHRC1_SULF1_summary.csv
9. /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/03_statistics/README_statistics_and_interpretation.txt

关键图：
- /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/plots/01_state_identity_risk_scores.png
- /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/plots/02_state_cell_retention_across_filters.png
- /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/plots/03_C3_core_candidate_nonEC_risk.png
- /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/plots/04_UMAP_strict_EC_filter.png（若h5ad含UMAP）
- /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/plots/05_C3_candidate_vs_EC_identity_margin.png
- /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/plots/06_candidate_edgeR_logFC_sensitivity_heatmap.png
- /root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity/plots/07_core_gene_patient_differences_across_filters.png

解释规则：
- ANGPT2/LOX/CTHRC1在S3_STRICT_EC中C3_vs_C7仍为正方向：支持早期C3信号不完全由明显非EC共表达风险驱动；
- SULF1在S3_STRICT_EC中Late_vs_C3仍为正方向：支持其晚期锚点定位；
- 严格过滤后方向消失：需降级为队列依赖或潜在污染/转化相关信号；
- 严格集合患者数不足：不能判定失败，应报告为证据不足。
