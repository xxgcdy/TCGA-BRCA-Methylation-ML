# 图表索引 / Figure guide

| 图 | 回答的问题 | 主要结论 | 注意事项 |
|---|---|---|---|
| `cohort_qc_flow.png` | 数据从下载到QC保留了多少？ | 890文件最终形成869个QC后样本 | 文件数不等于独立生物学样本 |
| `pca_pc1_pc2.png` | QC后样本的全局甲基化结构和离群情况如何？ | 正常组织较集中，肿瘤更分散；PC1/PC2解释27.26%/9.92%方差 | 不能单独排除批次、纯度或细胞组成影响 |
| `dmp_delta_beta_concordance.png` | 全队列和配对效应是否一致？ | Δβ高度一致，r=0.9924 | 高相关不消除混杂风险 |
| `paired_dmp_volcano.png` | 哪些CpG同时显著且大效应？ | 存在大量方向明确的大效应位点 | 点数巨大，不能把显著等同因果 |
| `top50_replicated_cpg_heatmap.png` | 可重复CpG能否区分样本？ | Top位点呈清晰肿瘤/正常模式 | 图用于描述，不是独立泛化证据 |
| `GO_top20.png` / `KEGG_top20.png` | 差异位点邻近基因富集什么功能？ | 黏附/信号等通路富集 | 富集不代表激活或抑制 |
| `internal_models_corrected.png` | 修正后的四模型内部测试表现如何？ | 分组CV XGBoost最高，Transformer仍接近 | 单次划分且仅16个正常测试样本 |
| `model_metrics_patient_bootstrap_ci.png` | 初始模型比较的不确定性有多大？ | AUC较稳，特异度区间较宽 | 对应初始模型版本；bootstrap不能替代外部验证 |
| `transformer_architecture.png` | 初始Transformer如何工作？ | 10-CpG patch、64维、2层 | 该图记录初始顺序位置实现；修正版已消除跨染色体patch |
| `masked_pretraining_ablation.png` | 预训练收益是否依赖标签比例和位置编码？ | 基因组位置设置小幅受益，其他设置不稳定 | 不支持“预训练普遍更优”的结论 |
| `external_validation_gse52865.png` | 冻结模型能否迁移到独立队列？ | 排序能力较强，但阈值校准发生漂移 | 57例外部样本不足以证明临床泛化 |
| `external_calibration_gse52865.png` | 外部队列中的冻结概率是否仍校准？ | 树模型校准较好，Transformer概率明显漂移 | 仅57例；五分位曲线和ECE是描述性估计 |
