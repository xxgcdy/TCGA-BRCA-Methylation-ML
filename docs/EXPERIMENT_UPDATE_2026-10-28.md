# 可靠性修复、Masked预训练与外部验证（2026-10-28）

## 1. 实验边界

- 所有模型继续使用仅由TCGA训练集筛选的1,000个CpG。
- 患者级训练/验证/测试划分保持不变（604/131/134）。
- 模型选择和阈值只能使用训练集或验证集；测试集和GSE52865标签不参与选择。
- 服务器任务遵守最多8个CPU核心限制；预训练在本地RTX 5090 D上运行。
- 本地GPU环境使用Python 3.10、PyTorch 2.11.0+cu128；RTX 5090 D的CUDA矩阵运算在实验前已实际验证。

## 2. XGBoost可靠性修复

旧模型出现`best_iteration=0`，测试AUC仅0.8904，不能作为可信XGBoost基线。修正版在训练集内执行确定性的5折患者分组交叉验证，以log-loss选择参数与轮数：`max_depth=2`、`eta=0.02`、382轮。训练CV AUC为0.9982，冻结测试AUC为0.9958；测试敏感度0.9831、特异度0.9375。

脚本：`scripts/06_model_training/run_xgboost_grouped_cv.R`

## 3. 染色体安全Transformer

每条染色体分别按hg19位置排序并切分10-CpG patch。总计111个patch，22个染色体末端patch使用padding，跨染色体patch为0。五个随机种子的平均测试AUC为0.9870 ± 0.0052；概率集成测试AUC为0.9873，特异度1.0000。

脚本：`scripts/07_transformer/train_chromosome_safe_transformer.py`

## 4. Masked Methylation Pretraining

预训练目标为30%随机patch的β值重建，且只使用TCGA训练样本。下游比较：

- 位置编码：none、learned sequence + chromosome、chromosome + hg19 coordinate Fourier features；
- 初始化：随机初始化、masked预训练；
- 标签比例：10%、25%、50%、100%；
- 三个随机种子，共72次微调。

基因组位置编码下，预训练相对随机初始化的平均测试AUC增益为：10%标签+0.0009、25%+0.0032、50%+0.0019、100%+0.0028。其他位置设置下增益不稳定并出现负值。因此本结果证明流程已实现，但没有证明当前小规模预训练普遍优于随机初始化。

脚本：`scripts/09_pretraining/run_masked_pretraining_experiments.py`

## 5. GSE52865独立外部验证

GSE52865包含40个乳腺肿瘤和17个正常乳腺组织HM450样本。下载文件为228,875,072字节，SHA256为`8cc1dded61438a1acd72769d0d5c6069450ee0b2d8e0de6504648b4137f7423b`。由于NCBI FTP数据通道发生静默重启，采用2,048个小分块、独立双副本比较和30个不一致块的第三次多数裁决，最终完整gzip流和CRC通过。冻结1,000个CpG全部命中、缺失率0。

| 模型 | ROC-AUC | Accuracy | Sensitivity | Specificity |
|---|---:|---:|---:|---:|
| LASSO | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| Random Forest | 1.0000 | 0.9649 | 1.0000 | 0.8824 |
| XGBoost grouped CV | 1.0000 | 0.9298 | 1.0000 | 0.7647 |
| Genomic Transformer, random | 0.9779 | 0.8947 | 0.9750 | 0.7059 |
| Genomic Transformer, pretrained | 0.9941 | 0.7368 | 0.6250 | 1.0000 |

ROC-AUC反映排序能力；accuracy/sensitivity/specificity使用TCGA验证集冻结阈值。预训练Transformer的高AUC与较低敏感度并不矛盾，而是概率分布跨队列漂移的表现。不能用GSE52865标签重新选阈值后再把结果称为独立测试。

## 6. 主要结论

1. XGBoost故障已修复，可靠基线与LASSO/RF一致强。
2. 跨染色体patch已消除，多种子结果比单次训练更可靠。
3. Masked预训练已在GPU完整实现，但当前收益小且依赖位置编码。
4. GSE52865验证支持排序泛化，同时暴露了阈值和概率校准的跨队列漂移。
5. 下一步的研究价值主要来自扩大无标签多队列语料、采用标签无关token词表及校准迁移，而不是继续在单一TCGA测试集追求更高小数位AUC。

## English summary

XGBoost was repaired with training-only patient-grouped cross-validation (test AUC 0.9958), and chromosome-safe patching removed all cross-chromosome patches (five-seed ensemble AUC 0.9873). Masked methylation pretraining was implemented on an RTX 5090 D across three positional encodings and four label fractions. Gains were small and configuration-dependent. Frozen external evaluation on GSE52865 showed strong ranking transfer, but several models—especially the pretrained Transformer—exhibited threshold/calibration shift. No external-label retuning was performed.
