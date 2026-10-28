#!/usr/bin/env python3
import json
from datetime import datetime, timezone
import os
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
PROCESSED = PROJECT / "processed"
META = PROJECT / "metadata"
SEED = 20261027

samples = pd.read_csv(PROCESSED / "samples_after_qc.tsv", sep="\t")


def case_stratum(group: pd.DataFrame) -> str:
    types = set(group["sample_type"])
    if {"Primary Tumor", "Solid Tissue Normal"}.issubset(types):
        return "paired_tumor_normal"
    if "Primary Tumor" in types:
        return "tumor_only"
    if "Solid Tissue Normal" in types:
        return "normal_only"
    raise ValueError(f"Unexpected sample types: {types}")


cases = (
    samples.groupby("case_id", sort=True)
    .apply(lambda g: pd.Series({
        "case_submitter_id": g["case_submitter_id"].iloc[0],
        "stratum": case_stratum(g),
        "n_samples": len(g),
        "n_tumor": int((g["sample_type"] == "Primary Tumor").sum()),
        "n_normal": int((g["sample_type"] == "Solid Tissue Normal").sum()),
    }), include_groups=False)
    .reset_index()
)

rng = np.random.default_rng(SEED)
assignments = []
for stratum, group in cases.groupby("stratum", sort=True):
    ids = group["case_id"].to_numpy(copy=True)
    rng.shuffle(ids)
    n = len(ids)
    n_train = int(np.floor(0.70 * n))
    n_val = int(np.floor(0.15 * n))
    if n >= 3:
        n_train = max(n_train, 1)
        n_val = max(n_val, 1)
    n_test = n - n_train - n_val
    if n >= 3 and n_test < 1:
        n_test = 1
        n_train -= 1
    for split, subset in (
        ("train", ids[:n_train]),
        ("validation", ids[n_train:n_train + n_val]),
        ("test", ids[n_train + n_val:]),
    ):
        assignments.extend({"case_id": case_id, "split": split} for case_id in subset)

case_split = cases.merge(pd.DataFrame(assignments), on="case_id", how="left", validate="one_to_one")
sample_split = samples.merge(case_split[["case_id", "stratum", "split"]], on="case_id", how="left", validate="many_to_one")

if sample_split["split"].isna().any():
    raise RuntimeError("Unassigned samples")
sets = {s: set(case_split.loc[case_split["split"] == s, "case_id"]) for s in ["train", "validation", "test"]}
if (sets["train"] & sets["validation"]) or (sets["train"] & sets["test"]) or (sets["validation"] & sets["test"]):
    raise RuntimeError("Patient leakage across splits")

case_split.to_csv(META / "patient_split.tsv", sep="\t", index=False)
sample_split.to_csv(META / "sample_split.tsv", sep="\t", index=False)
for split in ["train", "validation", "test"]:
    sample_split[sample_split["split"] == split].to_csv(META / f"{split}_samples.tsv", sep="\t", index=False)

table = (
    sample_split.groupby(["split", "sample_type"]).size().rename("samples").reset_index()
)
case_table = (
    case_split.groupby(["split", "stratum"]).size().rename("cases").reset_index()
)
table.to_csv(META / "split_sample_counts.tsv", sep="\t", index=False)
case_table.to_csv(META / "split_case_counts.tsv", sep="\t", index=False)

summary = {
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "seed": SEED,
    "strategy": "patient-level stratified by paired_tumor_normal, tumor_only, normal_only",
    "cases": int(len(case_split)),
    "samples": int(len(sample_split)),
    "patient_overlap_between_splits": 0,
    "split_sample_counts": table.to_dict(orient="records"),
    "split_case_counts": case_table.to_dict(orient="records"),
}
(META / "split_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2))

