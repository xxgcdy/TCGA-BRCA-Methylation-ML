#!/usr/bin/env python3
import json
from datetime import datetime, timezone
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse.linalg import svds

PROJECT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
PROCESSED = PROJECT / "processed"
QC = PROJECT / "results" / "qc"

matrix = np.load(PROCESSED / "beta_matrix_float32.npy", mmap_mode="r")
sample_qc = pd.read_csv(QC / "sample_filter_decisions.tsv", sep="\t")
probe = pd.read_csv(QC / "probe_filter_decisions.tsv.gz", sep="\t", low_memory=False)
probe_stats = pd.read_csv(QC / "probe_qc_unfiltered.tsv.gz", sep="\t")
if not np.array_equal(probe["matrix_column"].to_numpy(), probe_stats["matrix_column"].to_numpy()):
    raise RuntimeError("Probe QC table order mismatch")
probe["sd_beta"] = probe_stats["sd_beta"].to_numpy()

sample_indices = np.flatnonzero(sample_qc["qc_keep"].to_numpy(dtype=bool))
eligible = (
    probe["qc_keep_autosomal"].to_numpy(dtype=bool)
    & (probe["missing_rate_after_sample_qc"].to_numpy() <= 0.05)
)
eligible_idx = np.flatnonzero(eligible)
variability = probe.loc[eligible, "sd_beta"].to_numpy()
n_features = min(20000, len(eligible_idx))
top_local = np.argpartition(variability, -n_features)[-n_features:]
probe_indices = eligible_idx[top_local[np.argsort(variability[top_local])[::-1]]]

x = np.asarray(matrix[np.ix_(sample_indices, probe_indices)], dtype=np.float32)
medians = np.nanmedian(x, axis=0)
missing = np.isnan(x)
x[missing] = np.take(medians, np.nonzero(missing)[1])
means = x.mean(axis=0, dtype=np.float64)
sd = x.std(axis=0, ddof=1, dtype=np.float64)
sd[sd == 0] = 1
x = ((x - means) / sd).astype(np.float32)

k = 10
u, singular, vt = svds(x, k=k, which="LM", return_singular_vectors=True)
order = np.argsort(singular)[::-1]
singular = singular[order]
u = u[:, order]
scores = u * singular
total_variance = np.var(x, axis=0, ddof=1).sum()
explained = (singular * singular / (x.shape[0] - 1)) / total_variance

kept = sample_qc.loc[sample_indices].reset_index(drop=True).copy()
for i in range(k):
    kept[f"PC{i+1}"] = scores[:, i]

outlier = np.zeros(len(kept), dtype=bool)
max_abs_robust_z = np.zeros(len(kept), dtype=float)
for sample_type, idx in kept.groupby("sample_type").groups.items():
    block = scores[np.asarray(list(idx)), :]
    center = np.median(block, axis=0)
    mad = np.median(np.abs(block - center), axis=0)
    scale = 1.4826 * mad
    scale[scale == 0] = np.inf
    rz = np.abs((block - center) / scale)
    maxz = rz.max(axis=1)
    max_abs_robust_z[np.asarray(list(idx))] = maxz
    outlier[np.asarray(list(idx))] = maxz > 6

kept["max_abs_robust_z_pc1_pc10_within_type"] = max_abs_robust_z
kept["pca_outlier_flag_gt6"] = outlier
kept.to_csv(QC / "pca_sample_scores.tsv", sep="\t", index=False)
pd.DataFrame({
    "rank": np.arange(1, n_features + 1),
    "matrix_column": probe.loc[probe_indices, "matrix_column"].to_numpy(),
    "probe_id": probe.loc[probe_indices, "probe_id"].to_numpy(),
    "sd_beta_unfiltered": probe.loc[probe_indices, "sd_beta"].to_numpy(),
}).to_csv(QC / "pca_probe_set_top20000.tsv.gz", sep="\t", index=False, compression="gzip")

colors = {"Primary Tumor": "#d62728", "Solid Tissue Normal": "#1f77b4"}
fig, ax = plt.subplots(figsize=(8, 6))
for sample_type, group in kept.groupby("sample_type"):
    ax.scatter(group["PC1"], group["PC2"], s=18, alpha=0.7, label=sample_type, c=colors[sample_type])
flagged = kept[kept["pca_outlier_flag_gt6"]]
if len(flagged):
    ax.scatter(flagged["PC1"], flagged["PC2"], s=70, facecolors="none", edgecolors="black", linewidths=1, label="PCA outlier flag")
ax.set_xlabel(f"PC1 ({explained[0]*100:.2f}%)")
ax.set_ylabel(f"PC2 ({explained[1]*100:.2f}%)")
ax.set_title("TCGA-BRCA HM450 QC PCA")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(QC / "pca_pc1_pc2.png", dpi=180)
fig.savefig(QC / "pca_pc1_pc2.pdf")
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for sample_type, group in sample_qc.groupby("sample_type"):
    axes[0].hist(group["missing_rate"], bins=30, alpha=0.55, label=sample_type, color=colors[sample_type])
axes[0].axvline(0.20, color="black", linestyle="--", label="QC cutoff")
axes[0].set_xlabel("Sample missing rate")
axes[0].set_ylabel("Samples")
axes[0].legend(frameon=False, fontsize=8)
axes[1].hist(probe["missing_rate_after_sample_qc"], bins=50, color="#666666")
axes[1].axvline(0.20, color="black", linestyle="--")
axes[1].set_xlabel("Probe missing rate after sample QC")
axes[1].set_ylabel("Probes")
fig.tight_layout()
fig.savefig(QC / "missingness_distributions.png", dpi=180)
fig.savefig(QC / "missingness_distributions.pdf")
plt.close(fig)

summary = {
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "samples": int(len(kept)),
    "eligible_autosomal_probes_missing_le_0.05": int(len(eligible_idx)),
    "pca_probes": int(n_features),
    "imputed_values_in_pca_matrix": int(missing.sum()),
    "pc_explained_variance_ratio": {f"PC{i+1}": float(explained[i]) for i in range(k)},
    "pca_outlier_flags_within_sample_type_max_abs_robust_z_gt6": int(outlier.sum()),
}
(QC / "pca_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2))
if outlier.any():
    print(kept.loc[outlier, ["sample_submitter_id", "sample_type", "PC1", "PC2", "max_abs_robust_z_pc1_pc10_within_type"]].to_string(index=False))

