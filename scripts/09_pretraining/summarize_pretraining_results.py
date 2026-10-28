#!/usr/bin/env python3
"""Summarize matched-seed pretraining gains and create an experiment figure."""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
RES = ROOT / "results/masked_pretraining"
FIG = ROOT / "results/figures"
FIG.mkdir(parents=True, exist_ok=True)

metrics = pd.read_csv(RES / "finetune_metrics_all_runs.tsv", sep="\t").query("split == 'test'")
history = pd.read_csv(RES / "pretraining_history.tsv", sep="\t")

pivot = metrics.pivot_table(index=["position_mode", "label_fraction", "seed"],
                            columns="initialization", values="roc_auc").reset_index()
pivot["pretraining_gain_auc"] = pivot["pretrained"] - pivot["random"]
gain = pivot.groupby(["position_mode", "label_fraction"])["pretraining_gain_auc"].agg(
    n_seeds="count", mean_gain="mean", sd_gain="std", min_gain="min", max_gain="max").reset_index()
gain.to_csv(RES / "pretraining_gain_paired.tsv", sep="\t", index=False)

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
for mode, frame in history.groupby("position_mode"):
    axes[0].plot(frame["epoch"], frame["masked_mse"], label=mode)
axes[0].set(xlabel="Pretraining epoch", ylabel="Masked reconstruction MSE",
            title="Masked methylation pretraining")
axes[0].grid(alpha=0.25); axes[0].legend(title="Position encoding")

colors = {"random": "#6B7280", "pretrained": "#0072B2"}
markers = {"none": "o", "sequence": "s", "genomic": "^"}
summary = metrics.groupby(["position_mode", "initialization", "label_fraction"])["roc_auc"].agg(["mean", "std"]).reset_index()
for (mode, init), frame in summary.groupby(["position_mode", "initialization"]):
    axes[1].errorbar(frame["label_fraction"] * 100, frame["mean"], yerr=frame["std"],
                     label=f"{mode} / {init}", color=colors[init], marker=markers[mode],
                     linestyle="-" if init == "pretrained" else "--", capsize=3, alpha=0.9)
axes[1].set(xlabel="Labeled training patients used (%)", ylabel="Held-out TCGA test ROC-AUC",
            title="Low-label fine-tuning (mean ± SD; 3 seeds)")
axes[1].set_xticks([10, 25, 50, 100]); axes[1].grid(alpha=0.25)
axes[1].legend(fontsize=7, ncol=2)
fig.suptitle("Position encoding and masked-pretraining ablation", fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "masked_pretraining_ablation.png", dpi=180, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(gain.to_string(index=False))
