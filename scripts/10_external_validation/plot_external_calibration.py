#!/usr/bin/env python3
"""Summarize frozen-probability calibration on the GSE52865 cohort."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, log_loss

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results/external_validation/GSE52865"
FIGURES = ROOT / "results/figures"
FIGURES.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20261029)

baseline = pd.read_csv(RESULTS / "baseline_external_predictions.tsv.gz", sep="\t")
deep = pd.read_csv(RESULTS / "pretraining_ablation_external_predictions.tsv.gz", sep="\t")
deep["model"] = deep["position_mode"] + "_" + deep["initialization"]
predictions = pd.concat(
    [
        baseline[["model", "label", "probability"]],
        deep[["model", "label", "probability"]],
    ],
    ignore_index=True,
)

selected = ["LASSO", "RandomForest", "XGBoostGroupedCV", "genomic_random", "genomic_pretrained"]
labels = ["LASSO", "Random forest", "XGBoost", "Genomic\nrandom", "Genomic\npretrained"]
colors = ["#4C78A8", "#59A14F", "#F28E2B", "#9C755F", "#B279A2"]


def calibration_bins(y, probability, n_bins=5):
    frame = pd.DataFrame({"label": y, "probability": probability})
    frame["bin"] = pd.qcut(frame["probability"], q=n_bins, labels=False, duplicates="drop")
    return frame.groupby("bin", observed=True).agg(
        mean_probability=("probability", "mean"),
        observed_fraction=("label", "mean"),
        n=("label", "size"),
    ).reset_index(drop=True)


metric_rows = []
curves = {}
for model in selected:
    frame = predictions[predictions.model == model].copy()
    y = frame.label.to_numpy(int)
    probability = frame.probability.to_numpy(float)
    curve = calibration_bins(y, probability)
    curves[model] = curve
    ece = np.sum(curve.n / len(frame) * np.abs(curve.mean_probability - curve.observed_fraction))
    metric_rows.append(
        {
            "model": model,
            "n": len(frame),
            "prevalence": y.mean(),
            "mean_predicted_probability": probability.mean(),
            "brier_score": brier_score_loss(y, probability),
            "log_loss": log_loss(y, np.clip(probability, 1e-7, 1 - 1e-7)),
            "ece_5_quantile_bins": ece,
        }
    )

metrics = pd.DataFrame(metric_rows)
metrics.to_csv(
    RESULTS / "external_calibration_metrics.tsv",
    sep="\t",
    index=False,
    float_format="%.6f",
)

fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.4), gridspec_kw={"width_ratios": [1.45, 1.0]})

# Probability distributions by observed class. Jitter is deterministic.
for model_index, (model, color) in enumerate(zip(selected, colors)):
    frame = predictions[predictions.model == model]
    for class_index, truth in enumerate([0, 1]):
        values = frame.loc[frame.label == truth, "probability"].to_numpy()
        position = model_index + (-0.16 if truth == 0 else 0.16)
        bp = axes[0].boxplot(values, positions=[position], widths=0.25, patch_artist=True,
                             showfliers=False, manage_ticks=False)
        bp["boxes"][0].set(facecolor="#72B7B2" if truth == 0 else "#E45756", alpha=0.55)
        for key in ["whiskers", "caps", "medians"]:
            for artist in bp[key]: artist.set(color="#333333", linewidth=1)
        jitter = RNG.uniform(-0.055, 0.055, size=len(values))
        axes[0].scatter(position + jitter, values, s=13, alpha=0.55,
                        color="#2A9D8F" if truth == 0 else "#D62728", edgecolors="none")
axes[0].set_xticks(np.arange(len(labels)), labels)
axes[0].set_ylabel("Frozen predicted tumor probability")
axes[0].set_ylim(-0.04, 1.04)
axes[0].grid(axis="y", alpha=0.22)
axes[0].set_title("Probability distributions by observed class")
axes[0].scatter([], [], color="#2A9D8F", label="Normal")
axes[0].scatter([], [], color="#D62728", label="Tumor")
axes[0].legend(frameon=False, loc="upper left")

axes[1].plot([0, 1], [0, 1], linestyle="--", color="#777777", label="Ideal calibration")
for model, label, color in zip(selected, labels, colors):
    curve = curves[model]
    axes[1].plot(curve.mean_probability, curve.observed_fraction, marker="o", lw=1.6,
                 color=color, label=label.replace("\n", " "))
axes[1].set_xlabel("Mean predicted probability")
axes[1].set_ylabel("Observed tumor fraction")
axes[1].set_xlim(-0.03, 1.03); axes[1].set_ylim(-0.03, 1.03)
axes[1].grid(alpha=0.22)
axes[1].set_title("Five-bin reliability curves")
axes[1].legend(frameon=False, fontsize=8, loc="upper left")

fig.suptitle("GSE52865 frozen-probability calibration (40 tumor, 17 normal)",
             fontsize=14, fontweight="bold")
fig.text(0.5, 0.015, "External labels are used for evaluation only; no threshold or probability recalibration was fitted.",
         ha="center", fontsize=9, color="#555555")
fig.tight_layout(rect=[0, 0.04, 1, 0.94])
fig.savefig(FIGURES / "external_calibration_gse52865.png", dpi=220, bbox_inches="tight")
plt.close(fig)

print(metrics.to_string(index=False))
