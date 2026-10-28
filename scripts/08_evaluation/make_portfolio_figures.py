#!/usr/bin/env python3
"""Generate portfolio figures from real QC, DMP, matrix, and model outputs."""

import json
import os
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
OUT = ROOT / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

COLORS = {"tumor": "#D55E00", "normal": "#0072B2", "accent": "#009E73", "grey": "#6B7280"}


def save(fig, name):
    fig.savefig(OUT / name, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# Cohort and probe flow from saved QC summaries.
prep = json.loads((ROOT / "results/qc/sample_preparation_summary.json").read_text())
qc = json.loads((ROOT / "results/qc/qc_filter_summary.json").read_text())
fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
sample_labels = ["GDC files", "Unique biological\nsamples", "QC-passed\nsamples"]
sample_values = [prep["input_files"], prep["unique_samples"], qc["retained_samples"]]
axes[0].bar(sample_labels, sample_values, color=["#9CA3AF", "#56B4E9", COLORS["accent"]])
for i, value in enumerate(sample_values):
    axes[0].text(i, value + 10, f"{value:,}", ha="center", fontweight="bold")
axes[0].set_ylim(0, 980); axes[0].set_ylabel("Count"); axes[0].set_title("Sample flow")
probe_labels = ["HM450 rows", "All-chromosome\nQC probes", "Autosomal\nanalysis probes"]
probe_values = [qc["input_probes"], qc["retained_probes_all_chromosomes"], qc["retained_probes_autosomal"]]
axes[1].bar(probe_labels, probe_values, color=["#9CA3AF", "#56B4E9", COLORS["accent"]])
for i, value in enumerate(probe_values):
    axes[1].text(i, value + 6000, f"{value:,}", ha="center", fontweight="bold")
axes[1].set_ylim(0, 540000); axes[1].set_ylabel("CpG probes"); axes[1].set_title("Probe filtering")
fig.suptitle("TCGA-BRCA HM450 cohort construction", fontweight="bold", y=0.98)
fig.text(0.27, 0.035, "890 files → 888 unique samples → 869 QC-passed samples",
         ha="center", fontsize=9, color="#374151")
fig.text(0.27, 0.012, "2 technical replicates removed; final cohort: 773 tumor / 96 normal",
         ha="center", fontsize=8.5, color="#374151")
fig.subplots_adjust(left=0.08, right=0.98, top=0.82, bottom=0.22, wspace=0.20)
save(fig, "cohort_qc_flow.png")


# Full-versus-paired delta-beta concordance.
full = pd.read_csv(ROOT / "results/dmp/dmp_full_cohort.tsv.gz", sep="\t",
                   usecols=["probe_id", "adj.P.Val", "delta_beta"])
paired = pd.read_csv(ROOT / "results/dmp/dmp_paired_difference.tsv.gz", sep="\t",
                     usecols=["probe_id", "adj.P.Val", "mean_delta_beta"])
dmp = full.merge(paired, on="probe_id", suffixes=("_full", "_paired"), validate="one_to_one")
r = np.corrcoef(dmp["delta_beta"], dmp["mean_delta_beta"])[0, 1]
common = ((dmp["adj.P.Val_full"] < 0.05) & (dmp["delta_beta"].abs() >= 0.20) &
          (dmp["adj.P.Val_paired"] < 0.05) & (dmp["mean_delta_beta"].abs() >= 0.20) &
          (np.sign(dmp["delta_beta"]) == np.sign(dmp["mean_delta_beta"])))
fig, ax = plt.subplots(figsize=(6.2, 5.5))
hb = ax.hexbin(dmp["delta_beta"], dmp["mean_delta_beta"], gridsize=100, bins="log", cmap="viridis", mincnt=1)
lim = 0.9
ax.plot([-lim, lim], [-lim, lim], "--", color="white", linewidth=1.2)
ax.set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel="Full cohort Δβ (tumor − normal)",
       ylabel="Paired mean Δβ (tumor − normal)", title="Effect-size concordance")
ax.text(0.03, 0.96, f"Pearson r = {r:.4f}\nCommon large-effect CpGs = {common.sum():,}\nAll CpGs = {len(dmp):,}",
        transform=ax.transAxes, va="top", bbox=dict(boxstyle="round", facecolor="white", alpha=0.9))
fig.colorbar(hb, ax=ax, label="log10 bin count")
save(fig, "dmp_delta_beta_concordance.png")


# Heatmap of 50 strongest replicated large-effect CpGs.
dmp_common = dmp.loc[common].copy()
dmp_common["replicated_effect"] = np.minimum(dmp_common["delta_beta"].abs(), dmp_common["mean_delta_beta"].abs())
top_ids = dmp_common.nlargest(50, "replicated_effect")["probe_id"].tolist()
probes = pd.read_csv(ROOT / "processed/probe_order_qc_autosomal.tsv.gz", sep="\t", usecols=["probe_id"])
indices = probes.reset_index().set_index("probe_id").loc[top_ids, "index"].to_numpy()
samples = pd.read_csv(ROOT / "processed/sample_order_qc_autosomal.tsv", sep="\t", usecols=["sample_type"])
with h5py.File(ROOT / "processed/beta_qc_autosomal.h5", "r") as h:
    read_order = np.argsort(indices)
    restore = np.argsort(read_order)
    beta = h["beta"][:, indices[read_order]][:, restore].astype(np.float32)
order = np.argsort(samples["sample_type"].eq("Primary Tumor").to_numpy(), kind="stable")
beta = beta[order].T
row_med = np.nanmedian(beta, axis=1, keepdims=True)
beta = np.where(np.isnan(beta), row_med, beta)
z = (beta - beta.mean(axis=1, keepdims=True)) / np.maximum(beta.std(axis=1, keepdims=True), 1e-6)
z = np.clip(z, -3, 3)
fig, ax = plt.subplots(figsize=(12, 7.2))
im = ax.imshow(z, aspect="auto", interpolation="nearest", cmap="RdBu_r", vmin=-3, vmax=3)
n_normal = int((~samples["sample_type"].eq("Primary Tumor")).sum())
ax.axvline(n_normal - 0.5, color="black", linewidth=1)
ax.set(xlabel=f"Samples ordered by tissue type (normal n={n_normal}; tumor n={len(samples)-n_normal})",
       ylabel="Top replicated CpGs", title="Top 50 replicated large-effect CpGs")
ax.set_yticks(np.arange(50), top_ids, fontsize=5.5); ax.set_xticks([])
fig.colorbar(im, ax=ax, label="Per-CpG beta z-score", fraction=0.02)
save(fig, "top50_replicated_cpg_heatmap.png")


# Model estimates with patient-bootstrap intervals.
ci = pd.read_csv(ROOT / "results/evaluation/patient_bootstrap_95ci.tsv", sep="\t")
show_metrics = ["roc_auc", "sensitivity", "specificity"]
models = ci["model"].drop_duplicates().tolist()
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8), sharey=True)
for ax, metric in zip(axes, show_metrics):
    zc = ci[ci["metric"] == metric].set_index("model").loc[models]
    x = np.arange(len(models))
    ax.errorbar(x, zc["estimate"], yerr=[zc["estimate"] - zc["ci_low"], zc["ci_high"] - zc["estimate"]],
                fmt="o", capsize=4, color=COLORS["normal"])
    ax.set_xticks(x, [m.replace("CpGPatchTransformer", "Transformer") for m in models], rotation=35, ha="right")
    ax.set_title(metric.replace("_", " ").title()); ax.grid(axis="y", alpha=0.25)
    ax.set_ylim(0.55, 1.02)
axes[0].set_ylabel("Estimate and 95% patient-bootstrap CI")
fig.suptitle("Held-out test performance uncertainty (2,000 patient bootstraps)", fontweight="bold")
save(fig, "model_metrics_patient_bootstrap_ci.png")


# Architecture schematic derived from the implemented tensor flow.
fig, ax = plt.subplots(figsize=(15.5, 4.6)); ax.axis("off")
boxes = [
    ("β values", "[B, 1000]"),
    ("Genomic\nordering", "chr + hg19 position"),
    ("Fixed CpG\npatches", "[B, 100, 10]"),
    ("Linear projection\n+ LN + GELU", "[B, 100, 64]"),
    ("Chromosome embedding\n+ sequence position", "+ [CLS] token"),
    ("2-layer\nTransformer", "4 heads; FFN = 128"),
    ("MLP classifier", "[B] tumor logit"),
]
box_width = 0.122
xs = np.linspace(0.025, 0.855, len(boxes))
for i, ((title, subtitle), x0) in enumerate(zip(boxes, xs)):
    ax.add_patch(plt.Rectangle((x0, 0.34), box_width, 0.36, facecolor="#E8F1FA", edgecolor="#2563EB", linewidth=1.2))
    ax.text(x0 + box_width / 2, 0.56, title, ha="center", va="center",
            fontsize=8.2, fontweight="bold", linespacing=1.15)
    ax.text(x0 + box_width / 2, 0.405, subtitle, ha="center", va="center",
            fontsize=7.1, color="#374151")
    if i < len(boxes) - 1:
        ax.annotate("", xy=(xs[i + 1] - 0.004, 0.52), xytext=(x0 + box_width + 0.004, 0.52),
                    arrowprops=dict(arrowstyle="->", color="#4B5563", linewidth=1.1))
ax.text(0.5, 0.92, "Implemented CpG-patch Transformer (78,017 trainable parameters)",
        ha="center", fontsize=13, fontweight="bold")
ax.text(0.5, 0.14, "Baseline limitations", ha="center", fontsize=9.2,
        color="#991B1B", fontweight="bold")
ax.text(0.5, 0.085, "Sequence positions do not encode genomic distance  •  Fixed patches may cross chromosome boundaries",
        ha="center", fontsize=8.7, color="#991B1B")
save(fig, "transformer_architecture.png")

print("Portfolio figures written to", OUT)
