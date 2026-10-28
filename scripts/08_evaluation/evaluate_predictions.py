#!/usr/bin/env python3
"""Recompute held-out metrics and patient-bootstrap confidence intervals.

This script uses saved per-sample probabilities. Bootstrap sampling is performed
at the patient level, keeping paired tumor/normal samples from a patient together.
"""

import hashlib
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
SEED = 20261028
N_BOOTSTRAP = int(os.environ.get("N_BOOTSTRAP", "2000"))


def roc_auc(y: np.ndarray, p: np.ndarray) -> float:
    y = np.asarray(y, dtype=int)
    ranks = pd.Series(p).rank(method="average").to_numpy()
    n_pos, n_neg = int(y.sum()), int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    return float((ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def average_precision(y: np.ndarray, p: np.ndarray) -> float:
    order = np.argsort(-np.asarray(p), kind="mergesort")
    yy = np.asarray(y, dtype=int)[order]
    return float(((np.cumsum(yy) / np.arange(1, len(yy) + 1)) * yy).sum() / yy.sum())


def metrics(y: np.ndarray, p: np.ndarray, threshold: float) -> dict:
    y = np.asarray(y, dtype=int)
    pred = (np.asarray(p) >= threshold).astype(int)
    tp = int(((pred == 1) & (y == 1)).sum())
    tn = int(((pred == 0) & (y == 0)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum())
    return {
        "roc_auc": roc_auc(y, p),
        "average_precision": average_precision(y, p),
        "accuracy": float((pred == y).mean()),
        "sensitivity": tp / (tp + fn),
        "specificity": tn / (tn + fp),
        "tp": tp, "tn": tn, "fp": fp, "fn": fn,
    }


def patient_bootstrap(frame: pd.DataFrame, threshold: float, rng: np.random.Generator) -> pd.DataFrame:
    cases = frame["case_id"].drop_duplicates().to_numpy()
    grouped = {case: g for case, g in frame.groupby("case_id", sort=False)}
    rows = []
    for _ in range(N_BOOTSTRAP):
        sampled = rng.choice(cases, size=len(cases), replace=True)
        boot = pd.concat([grouped[case] for case in sampled], ignore_index=True)
        if boot["label"].nunique() < 2:
            continue
        rows.append(metrics(boot["label"].to_numpy(), boot["probability"].to_numpy(), threshold))
    return pd.DataFrame(rows)


def anon(value: str) -> str:
    return hashlib.sha256(str(value).encode("utf-8")).hexdigest()[:16]


out = ROOT / "results" / "evaluation"
fig_dir = ROOT / "results" / "figures"
out.mkdir(parents=True, exist_ok=True)
fig_dir.mkdir(parents=True, exist_ok=True)

baseline = pd.read_csv(ROOT / "results/models/baseline_predictions.tsv.gz", sep="\t")
transformer = pd.read_csv(ROOT / "results/transformer/transformer_predictions.tsv.gz", sep="\t")
predictions = pd.concat([baseline, transformer], ignore_index=True)
base_metrics = pd.read_csv(ROOT / "results/models/baseline_metrics.tsv", sep="\t")
tf_metrics = pd.read_csv(ROOT / "results/transformer/transformer_metrics.tsv", sep="\t")
reported = pd.concat([base_metrics, tf_metrics.assign(model="CpGPatchTransformer")], ignore_index=True)

recomputed, ci_rows = [], []
rng = np.random.default_rng(SEED)
test_frames = {}
for model, frame in predictions[predictions["split"] == "test"].groupby("model", sort=False):
    threshold = float(reported.loc[(reported["model"] == model) & (reported["split"] == "test"), "threshold"].iloc[0])
    point = metrics(frame["label"].to_numpy(), frame["probability"].to_numpy(), threshold)
    recomputed.append({"model": model, "split": "test", "n": len(frame), "threshold": threshold, **point})
    boot = patient_bootstrap(frame, threshold, rng)
    for metric in ["roc_auc", "average_precision", "accuracy", "sensitivity", "specificity"]:
        low, high = np.nanpercentile(boot[metric], [2.5, 97.5])
        ci_rows.append({"model": model, "metric": metric, "estimate": point[metric],
                        "ci_low": float(low), "ci_high": float(high),
                        "bootstrap_unit": "patient", "n_bootstrap": len(boot)})
    test_frames[model] = (frame, threshold, point)

recomputed_df = pd.DataFrame(recomputed)
ci_df = pd.DataFrame(ci_rows)
recomputed_df.to_csv(out / "model_metrics_recomputed.tsv", sep="\t", index=False)
ci_df.to_csv(out / "patient_bootstrap_95ci.tsv", sep="\t", index=False)

public_predictions = predictions.copy()
public_predictions["sample_id"] = public_predictions["sample_id"].map(anon)
public_predictions["case_id"] = public_predictions["case_id"].map(anon)
public_predictions.to_csv(out / "anonymized_predictions.tsv.gz", sep="\t", index=False, compression="gzip")

models = list(test_frames)
fig, axes = plt.subplots(1, len(models), figsize=(3.2 * len(models), 3.2), constrained_layout=True)
for ax, model in zip(axes, models):
    _, _, m = test_frames[model]
    mat = np.array([[m["tn"], m["fp"]], [m["fn"], m["tp"]]])
    ax.imshow(mat, cmap="Blues")
    for i in range(2):
        for j in range(2):
            ax.text(j, i, str(mat[i, j]), ha="center", va="center", fontsize=12)
    ax.set_title(model)
    ax.set_xticks([0, 1], ["Normal", "Tumor"])
    ax.set_yticks([0, 1], ["Normal", "Tumor"])
    ax.set_xlabel("Predicted")
    if ax is axes[0]:
        ax.set_ylabel("Observed")
fig.suptitle("Held-out test confusion matrices (validation-selected thresholds)")
fig.savefig(fig_dir / "test_confusion_matrices.png", dpi=180)
plt.close(fig)

summary = {
    "seed": SEED,
    "bootstrap_unit": "patient",
    "n_bootstrap_requested": N_BOOTSTRAP,
    "models": models,
    "metrics_recomputed_from_saved_probabilities": True,
}
(out / "evaluation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(recomputed_df.to_string(index=False))
print(ci_df.to_string(index=False))
