#!/usr/bin/env python3
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
QC = PROJECT / "results" / "qc"

sample = pd.read_csv(QC / "sample_qc_unfiltered.tsv", sep="\t")
probe = pd.read_csv(QC / "probe_qc_unfiltered.tsv.gz", sep="\t")

quantiles = [0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1]
sample_quantiles = sample["missing_rate"].quantile(quantiles).to_dict()
probe_quantiles = probe["missing_rate"].quantile(quantiles).to_dict()

median = float(sample["missing_rate"].median())
mad = float(np.median(np.abs(sample["missing_rate"] - median)))
robust_upper = median + 3 * 1.4826 * mad

sample_thresholds = {}
for threshold in [0.10, 0.15, 0.20, 0.25, 0.30]:
    sample_thresholds[str(threshold)] = int((sample["missing_rate"] > threshold).sum())

probe_thresholds = {}
for threshold in [0, 0.01, 0.05, 0.10, 0.20, 0.30, 0.50]:
    keep = probe["missing_rate"] <= threshold
    probe_thresholds[str(threshold)] = {
        "retained": int(keep.sum()),
        "removed": int((~keep).sum()),
        "retained_percent": float(100 * keep.mean()),
    }

by_type = (
    sample.groupby("sample_type")["missing_rate"]
    .agg(["count", "min", "median", "mean", "max"])
    .reset_index()
)
by_type.to_csv(QC / "sample_missingness_by_type.tsv", sep="\t", index=False)

outliers = sample[sample["missing_rate"] > robust_upper].sort_values("missing_rate", ascending=False)
outliers.to_csv(QC / "sample_missingness_robust_outliers.tsv", sep="\t", index=False)

summary = {
    "samples": int(len(sample)),
    "probes": int(len(probe)),
    "sample_missing_rate_quantiles": {str(k): float(v) for k, v in sample_quantiles.items()},
    "probe_missing_rate_quantiles": {str(k): float(v) for k, v in probe_quantiles.items()},
    "sample_missing_rate_median": median,
    "sample_missing_rate_mad": mad,
    "sample_robust_upper_cutoff_median_plus_3_scaled_mad": float(robust_upper),
    "robust_outlier_samples": int(len(outliers)),
    "samples_above_fixed_threshold": sample_thresholds,
    "probe_retention_by_max_missing_rate": probe_thresholds,
    "probes_all_missing": int((probe["missing_rate"] == 1).sum()),
    "probes_zero_variance_among_observed": int((probe["sd_beta"] == 0).sum()),
}
(QC / "missingness_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2))
print("\nBY_SAMPLE_TYPE")
print(by_type.to_string(index=False))
print("\nTOP_SAMPLE_OUTLIERS")
print(outliers[["sample_submitter_id", "sample_type", "missing_rate", "mean_beta"]].head(20).to_string(index=False))

