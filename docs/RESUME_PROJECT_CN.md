# 中文简历项目描述

## A. 一行版

基于TCGA-BRCA 869例HM450组织甲基化样本，构建泄漏控制的DMP/DMR—功能富集—机器学习/Transformer端到端流程，并以匿名预测、患者级bootstrap和独立测试完成可复核评估。

## B. 三条版

- 完成890个GDC文件MD5校验、技术重复处理与HM450 QC，获得869个样本和390,644个常染色体CpG；结合90对患者配对分析识别25,975个方向一致的大效应CpG及2,986个高置信DMR。
- 采用患者级train/validation/test划分与训练集内limma筛选Top 1,000 CpG，避免配对样本跨集合和特征选择泄漏；比较LASSO、随机森林、XGBoost与78k参数CpG-patch Transformer。
- 在134例独立测试样本上，RF/LASSO ROC-AUC分别为0.9947/0.9926，Transformer为0.9799；公开匿名预测、混淆矩阵及患者级bootstrap置信区间，并制定独立GEO外部验证方案。

## C. AI for Science强调版

- 将甲基化探针按基因组顺序组织为10-CpG patch，设计64维、2层Transformer监督基线，并系统审计其顺序位置编码、跨染色体patch与小样本单种子限制。
- 建立传统统计、生物学解释和深度学习的统一比较框架；结果显示强组织甲基化信号下简单模型优于Transformer，提出以masked methylation预训练、少标签微调和跨队列泛化作为下一阶段可证伪问题。
- 通过证据校验脚本复核患者隔离、特征唯一性、四模型测试样本一致性和指标重算，使每条简历结论可追溯到代码与机器可读结果。
