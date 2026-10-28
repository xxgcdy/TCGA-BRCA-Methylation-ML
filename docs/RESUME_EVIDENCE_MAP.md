# 简历证据映射 / Resume evidence map

| 可写入简历的陈述 | 可核查证据 | 不应扩大的说法 |
|---|---|---|
| 构建TCGA-BRCA HM450端到端分析流程 | `scripts/01`–`08`、`scripts/run_pipeline.sh` | 不称临床生产系统 |
| 完成869样本、390,644 CpG的QC | `qc_filter_summary.json`、QC流程图 | 不声称从原始IDAT完成全套QC |
| 发现25,975个跨全队列/配对一致的大效应CpG | `dmp_combined_summary.json` | 不声称因果位点 |
| 获得2,986个高置信DMR | `dmr_overlap_summary.json` | 不称严格共识或功能验证DMR |
| 以正确背景完成GO/KEGG偏倚校正富集 | enrichment表、KEGG缓存脚本 | 不称通路被激活 |
| 实现患者级划分与训练集内特征选择 | evidence checks、训练特征表 | 不称完全嵌套CV已完成 |
| 比较LASSO/RF/XGBoost/Transformer | 匿名预测、重算指标、ROC/CI图 | 不据单次划分宣称某算法普遍最优 |
| 实现78k参数CpG-patch Transformer | 训练脚本、架构图、summary | 不称甲基化大模型或已预训练模型 |
| RF测试AUC 0.9947，Transformer 0.9799 | `model_metrics_recomputed.tsv` | 必须同时说明134个测试样本、16个正常和CI |
