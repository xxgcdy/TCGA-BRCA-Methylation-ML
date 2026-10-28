#!/usr/bin/env python3
import csv
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
SHEET = PROJECT / "metadata" / "sample_sheet_tcga_brca_hm450.tsv"
RAW = PROJECT / "raw_beta"
QC = PROJECT / "results" / "qc"
QC.mkdir(parents=True, exist_ok=True)

samples = pd.read_csv(SHEET, sep="\t", dtype=str, encoding="utf-8-sig")
samples["size"] = pd.to_numeric(samples["size"], errors="raise")
samples["relative_path"] = samples["id"] + "/" + samples["filename"]
samples["absolute_path"] = samples["relative_path"].map(lambda p: str(RAW / p))
samples["label"] = samples["sample_type"].map({"Primary Tumor": 1, "Solid Tissue Normal": 0})


def read_beta(path: str) -> pd.Series:
    frame = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["probe_id", "beta"],
        dtype={"probe_id": "string", "beta": "float64"},
        na_values=["NA"],
    )
    if frame["probe_id"].duplicated().any():
        raise ValueError(f"Duplicate probe IDs in {path}")
    return frame.set_index("probe_id")["beta"]


selected_ids = set(samples["id"])
reports = []
excluded = []
duplicate_groups = samples.groupby("sample_id", sort=True).filter(lambda x: len(x) > 1).groupby("sample_id")

for sample_id, group in duplicate_groups:
    loaded = {}
    file_metrics = {}
    for row in group.itertuples(index=False):
        beta = read_beta(row.absolute_path)
        loaded[row.id] = beta
        valid = beta.dropna()
        file_metrics[row.id] = {
            "file_id": row.id,
            "filename": row.filename,
            "aliquot_submitter_id": row.aliquot_submitter_id,
            "n_probes": int(beta.size),
            "n_missing": int(beta.isna().sum()),
            "missing_rate": float(beta.isna().mean()),
            "mean_beta": float(valid.mean()),
            "median_beta": float(valid.median()),
            "out_of_range": int(((valid < 0) | (valid > 1)).sum()),
        }

    ids = list(group["id"])
    if len(ids) != 2:
        raise ValueError(f"Expected exactly two files for duplicate sample {sample_id}")
    left, right = loaded[ids[0]].align(loaded[ids[1]], join="inner")
    mask = left.notna() & right.notna()
    delta = left[mask] - right[mask]
    pearson = float(left[mask].corr(right[mask], method="pearson"))
    spearman = float(left[mask].corr(right[mask], method="spearman"))
    rmse = float(np.sqrt(np.mean(np.square(delta.to_numpy()))))
    mad = float(np.mean(np.abs(delta.to_numpy())))

    ranked = sorted(ids, key=lambda fid: (file_metrics[fid]["n_missing"], -int(samples.loc[samples["id"] == fid, "size"].iloc[0]), fid))
    keep_id = ranked[0]
    drop_id = ranked[1]
    selected_ids.discard(drop_id)
    excluded.append({
        "sample_id": sample_id,
        "kept_file_id": keep_id,
        "excluded_file_id": drop_id,
        "reason": "technical replicate; retained file with fewer missing beta values",
    })
    for fid in ids:
        metric = file_metrics[fid]
        metric.update({
            "sample_id": sample_id,
            "case_submitter_id": group["case_submitter_id"].iloc[0],
            "paired_file_id": right.name if False else ids[1] if fid == ids[0] else ids[0],
            "overlap_nonmissing": int(mask.sum()),
            "pearson_r": pearson,
            "spearman_r": spearman,
            "mean_absolute_difference": mad,
            "rmse": rmse,
            "decision": "KEEP" if fid == keep_id else "EXCLUDE",
        })
        reports.append(metric)

final = samples[samples["id"].isin(selected_ids)].copy().sort_values(
    ["sample_type", "case_submitter_id", "sample_submitter_id", "id"]
)
if len(final) != 888 or final["sample_id"].nunique() != 888:
    raise RuntimeError(f"Unexpected final cohort size: rows={len(final)}, samples={final['sample_id'].nunique()}")

pd.DataFrame(reports).to_csv(QC / "technical_replicate_report.tsv", sep="\t", index=False)
pd.DataFrame(excluded).to_csv(PROJECT / "metadata" / "excluded_technical_replicates.tsv", sep="\t", index=False)
final.to_csv(PROJECT / "metadata" / "final_sample_sheet.tsv", sep="\t", index=False)

summary = {
    "input_files": int(len(samples)),
    "duplicate_biological_samples": int(samples.groupby("sample_id").size().gt(1).sum()),
    "excluded_technical_replicate_files": int(len(excluded)),
    "final_files": int(len(final)),
    "unique_samples": int(final["sample_id"].nunique()),
    "unique_cases": int(final["case_id"].nunique()),
    "primary_tumor_samples": int((final["sample_type"] == "Primary Tumor").sum()),
    "solid_tissue_normal_samples": int((final["sample_type"] == "Solid Tissue Normal").sum()),
    "cases_with_both_primary_tumor_and_normal": int(
        final.groupby("case_id")["sample_type"].agg(set).map(
            lambda s: {"Primary Tumor", "Solid Tissue Normal"}.issubset(s)
        ).sum()
    ),
}
(QC / "sample_preparation_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2))

