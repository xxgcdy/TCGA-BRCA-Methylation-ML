#!/usr/bin/env python3
"""Plot the corrected frozen-test comparison without sample-level exports."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
base = pd.read_csv(ROOT / "results/model_metrics/baseline_metrics.tsv", sep="\t")
base = base.query("split == 'test' and model in ['LASSO', 'RandomForest']").copy()
xgb = pd.read_csv(ROOT / "results/models/xgboost_grouped_cv/metrics.tsv", sep="\t").query("split == 'test'")
transformer = pd.read_csv(ROOT / "results/transformer_chromosome_safe/metrics_all_seeds.tsv", sep="\t").query("split == 'test' and seed == 'ensemble'")
transformer = transformer.assign(model="ChromosomeSafeTransformer")
data = pd.concat([base, xgb, transformer], ignore_index=True)
labels = ["LASSO", "Random forest", "XGBoost\ngrouped CV", "Transformer\n5-seed ensemble"]
colors = ["#4C78A8", "#59A14F", "#F28E2B", "#B279A2"]

fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.7))
x = range(len(data))
axes[0].bar(x, data.roc_auc, color=colors)
axes[0].set_xticks(list(x), labels); axes[0].set_ylim(0.95, 1.002)
axes[0].set_ylabel("Frozen test ROC-AUC"); axes[0].grid(axis="y", alpha=0.22)
for i, value in enumerate(data.roc_auc): axes[0].text(i, value + 0.001, f"{value:.4f}", ha="center", fontsize=9)

width = 0.24
for offset, metric, color in zip([-width, 0, width], ["accuracy", "sensitivity", "specificity"], ["#4C78A8", "#E45756", "#72B7B2"]):
    axes[1].bar([i + offset for i in x], data[metric], width, label=metric.capitalize(), color=color)
axes[1].set_xticks(list(x), labels); axes[1].set_ylim(0.90, 1.01)
axes[1].set_ylabel("Validation-threshold performance"); axes[1].grid(axis="y", alpha=0.22)
axes[1].legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.03))
fig.suptitle("Corrected held-out TCGA test comparison (n=134)", fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(ROOT / "results/figures/internal_models_corrected.png", dpi=220, bbox_inches="tight")
plt.close(fig)
