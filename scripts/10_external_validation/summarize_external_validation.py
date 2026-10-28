#!/usr/bin/env python3
"""Bootstrap and visualize frozen GSE52865 external-validation results."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results/external_validation/GSE52865"
FIGURES = ROOT / "results/figures"
FIGURES.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20261028)

baseline = pd.read_csv(RESULTS / "baseline_external_predictions.tsv.gz", sep="\t")
deep = pd.read_csv(RESULTS / "pretraining_ablation_external_predictions.tsv.gz", sep="\t")
deep["model"] = deep["position_mode"] + "_" + deep["initialization"]
predictions = pd.concat([
    baseline[["model", "sample_id", "label", "probability", "threshold"]],
    deep[["model", "sample_id", "label", "probability", "threshold"]],
], ignore_index=True)


def point_metrics(group):
    y = group.label.to_numpy(int)
    p = group.probability.to_numpy(float)
    pred = p >= float(group.threshold.iloc[0])
    return {
        "roc_auc": roc_auc_score(y, p),
        "accuracy": np.mean(pred == y),
        "sensitivity": np.mean(pred[y == 1] == 1),
        "specificity": np.mean(pred[y == 0] == 0),
    }


rows = []
for model, group in predictions.groupby("model", sort=False):
    group = group.reset_index(drop=True)
    normal = np.flatnonzero(group.label.to_numpy() == 0)
    tumor = np.flatnonzero(group.label.to_numpy() == 1)
    boot = []
    for _ in range(2000):
        selected = np.r_[RNG.choice(normal, len(normal), replace=True), RNG.choice(tumor, len(tumor), replace=True)]
        boot.append(point_metrics(group.iloc[selected]))
    point = point_metrics(group)
    row = {"model": model, "n": len(group), "n_tumor": len(tumor), "n_normal": len(normal), **point}
    for metric in point:
        values = [x[metric] for x in boot]
        row[f"{metric}_ci_low"] = np.quantile(values, 0.025)
        row[f"{metric}_ci_high"] = np.quantile(values, 0.975)
    rows.append(row)

summary = pd.DataFrame(rows)
summary.to_csv(RESULTS / "external_metrics_bootstrap_ci.tsv", sep="\t", index=False)

selected = ["LASSO", "RandomForest", "XGBoostGroupedCV", "genomic_random", "genomic_pretrained"]
labels = ["LASSO", "Random forest", "XGBoost (grouped CV)", "Genomic Transformer\nrandom init.", "Genomic Transformer\nmasked pretrained"]
short_labels = ["LASSO", "Random\nforest", "XGBoost", "Genomic\nrandom", "Genomic\npretrained"]
plot = summary.set_index("model").loc[selected]
colors = ["#4C78A8", "#59A14F", "#F28E2B", "#9C755F", "#B279A2"]

fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.3), gridspec_kw={"width_ratios": [1.0, 1.55]})
y = np.arange(len(plot))
auc = plot.roc_auc.to_numpy()
err = np.vstack([auc - plot.roc_auc_ci_low.to_numpy(), plot.roc_auc_ci_high.to_numpy() - auc])
axes[0].barh(y, auc, color=colors, alpha=0.9)
axes[0].errorbar(auc, y, xerr=err, fmt="none", ecolor="#222222", capsize=3, lw=1.1)
axes[0].set_yticks(y, labels); axes[0].invert_yaxis(); axes[0].set_xlim(0.85, 1.01)
axes[0].set_xlabel("ROC-AUC (stratified bootstrap 95% CI)")
axes[0].grid(axis="x", alpha=0.22)

metrics = ["accuracy", "sensitivity", "specificity"]
x = np.arange(len(plot)); width = 0.23
metric_colors = ["#4C78A8", "#E45756", "#72B7B2"]
for offset, metric, color in zip([-width, 0, width], metrics, metric_colors):
    axes[1].bar(x + offset, plot[metric], width, label=metric.capitalize(), color=color)
axes[1].set_xticks(x, short_labels); axes[1].set_ylim(0.55, 1.03)
axes[1].set_ylabel("Frozen-threshold performance")
axes[1].grid(axis="y", alpha=0.22)
axes[1].legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.02))
axes[1].text(0.5, -0.18, "Thresholds fixed on TCGA validation; no GSE52865 retuning",
             transform=axes[1].transAxes, ha="center", fontsize=8.5, color="#555555")

fig.suptitle("Independent external validation: GSE52865 (40 tumor, 17 normal)", fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0.05, 1, 0.94])
fig.savefig(FIGURES / "external_validation_gse52865.png", dpi=220, bbox_inches="tight")
plt.close(fig)
print(summary.to_string(index=False))
