#!/usr/bin/env python3
"""Validate split integrity and consistency of saved project evidence."""

import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
OUT = ROOT / "results" / "evaluation"
OUT.mkdir(parents=True, exist_ok=True)

checks = {}

split = pd.read_csv(ROOT / "metadata/sample_split.tsv", sep="\t")
case_split_counts = split.groupby("case_id")["split"].nunique()
checks["patient_overlap_between_splits"] = int((case_split_counts > 1).sum())
checks["split_sample_counts"] = split.groupby(["split", "sample_type"]).size().rename("n").reset_index().to_dict("records")
checks["n_samples"] = int(len(split))
checks["n_cases"] = int(split["case_id"].nunique())

features = pd.read_csv(ROOT / "results/models/training_only_top1000_features.tsv.gz", sep="\t")
checks["selected_features"] = int(len(features))
checks["selected_feature_ids_unique"] = bool(features["probe_id"].is_unique)
checks["selected_p_values_nondecreasing"] = bool(features["P.Value"].is_monotonic_increasing)

baseline = pd.read_csv(ROOT / "results/models/baseline_predictions.tsv.gz", sep="\t")
transformer = pd.read_csv(ROOT / "results/transformer/transformer_predictions.tsv.gz", sep="\t")
all_pred = pd.concat([baseline, transformer], ignore_index=True)
test = all_pred[all_pred["split"] == "test"]
test_sets = test.groupby("model")["sample_id"].apply(set)
reference = test_sets.iloc[0]
checks["test_sample_sets_identical_across_models"] = bool(all(s == reference for s in test_sets))
label_consistency = test.groupby("sample_id")["label"].nunique()
checks["test_labels_consistent_across_models"] = bool((label_consistency == 1).all())
checks["test_samples_per_model"] = {str(k): int(v) for k, v in test.groupby("model").size().items()}

with h5py.File(ROOT / "processed/transformer/tcga_brca_hm450_top1000.h5", "r") as h:
    beta = h["beta"][:]
    split_code = h["split_code"][:]
    checks["transformer_matrix_shape"] = list(beta.shape)
    checks["transformer_nan_after_imputation"] = int(np.isnan(beta).sum())
    checks["transformer_split_counts"] = {
        "train": int((split_code == 0).sum()),
        "validation": int((split_code == 1).sum()),
        "test": int((split_code == 2).sum()),
    }
    checks["transformer_feature_selection_attribute"] = str(h.attrs["feature_selection"])

required_truths = [
    checks["patient_overlap_between_splits"] == 0,
    checks["selected_features"] == 1000,
    checks["selected_feature_ids_unique"],
    checks["selected_p_values_nondecreasing"],
    checks["test_sample_sets_identical_across_models"],
    checks["test_labels_consistent_across_models"],
    checks["transformer_nan_after_imputation"] == 0,
    checks["transformer_matrix_shape"] == [869, 1000],
]
checks["all_required_checks_passed"] = bool(all(required_truths))

(OUT / "project_evidence_checks.json").write_text(json.dumps(checks, indent=2), encoding="utf-8")
print(json.dumps(checks, indent=2))
if not checks["all_required_checks_passed"]:
    raise SystemExit("One or more required integrity checks failed")
