# 图表索引 / Figure guide

| 图 | 回答的问题 | 主要结论 | 注意事项 |
|---|---|---|---|
| `cohort_qc_flow.png` | 数据从下载到QC保留了多少？ | 890文件最终形成869个QC后样本 | 文件数不等于独立生物学样本 |
| `dmp_delta_beta_concordance.png` | 全队列和配对效应是否一致？ | Δβ高度一致，r=0.9924 | 高相关不消除混杂风险 |
| `paired_dmp_volcano.png` | 哪些CpG同时显著且大效应？ | 存在大量方向明确的大效应位点 | 点数巨大，不能把显著等同因果 |
| `top50_replicated_cpg_heatmap.png` | 可重复CpG能否区分样本？ | Top位点呈清晰肿瘤/正常模式 | 图用于描述，不是独立泛化证据 |
| `GO_top20.png` / `KEGG_top20.png` | 差异位点邻近基因富集什么功能？ | 黏附/信号等通路富集 | 富集不代表激活或抑制 |
| `test_roc_all_models.png` | 四模型排序如何？ | RF/LASSO最高，Transformer接近 | 单次划分且正常样本少 |
| `model_metrics_patient_bootstrap_ci.png` | 排名的不确定性有多大？ | AUC较稳，特异度区间较宽 | bootstrap不能替代外部验证 |
| `test_confusion_matrices.png` | 错误发生在哪一类？ | 四模型灵敏度相同，特异度有差异 | 1个正常样本改变特异度6.25个百分点 |
| `transformer_architecture.png` | 当前Transformer实际怎么工作？ | 10-CpG patch、64维、2层 | 顺序位置而非基因组距离，patch可能跨染色体 |
