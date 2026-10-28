#!/usr/bin/env python3
"""Extract frozen TCGA model features and labels from GSE52865 series matrix."""

import csv
import gzip
import hashlib
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
SOURCE = ROOT / "external/GSE52865/GSE52865_series_matrix.txt.gz"
FEATURES = ROOT / "results/model_metrics/training_only_top1000_features.tsv.gz"
if not FEATURES.exists():
    FEATURES = ROOT / "results/models/training_only_top1000_features.tsv.gz"
OUT = ROOT / "external/GSE52865/processed"
OUT.mkdir(parents=True, exist_ok=True)

feature_table = pd.read_csv(FEATURES, sep="\t")
wanted = feature_table["probe_id"].astype(str).tolist()
wanted_set = set(wanted)

sha256 = hashlib.sha256()
with SOURCE.open("rb") as raw:
    for block in iter(lambda: raw.read(8 * 1024 * 1024), b""):
        sha256.update(block)

metadata = {}
values = {}
sample_ids = None
in_table = False
with gzip.open(SOURCE, "rt", encoding="utf-8", errors="replace", newline="") as handle:
    for line in handle:
        line = line.rstrip("\r\n")
        if not in_table:
            if line.startswith("!Sample_"):
                row = next(csv.reader([line], delimiter="\t", quotechar='"'))
                metadata.setdefault(row[0], []).append(row[1:])
            elif line == "!series_matrix_table_begin":
                in_table = True
            continue
        if sample_ids is None:
            row = next(csv.reader([line], delimiter="\t", quotechar='"'))
            if row[0] != "ID_REF": raise RuntimeError("Unexpected table header")
            sample_ids = row[1:]
            continue
        if line == "!series_matrix_table_end": break
        probe = line.split("\t", 1)[0].strip('"')
        if probe not in wanted_set: continue
        row = next(csv.reader([line], delimiter="\t", quotechar='"'))
        values[probe] = np.asarray([float(x) if x not in ("", "NA", "null") else np.nan for x in row[1:]], dtype=np.float32)

missing_features = [probe for probe in wanted if probe not in values]
matrix = np.stack([values.get(probe, np.full(len(sample_ids), np.nan, dtype=np.float32)) for probe in wanted], axis=1)

source = metadata["!Sample_source_name_ch1"][0]
titles = metadata["!Sample_title"][0]
if len(sample_ids) != len(source): raise RuntimeError("Metadata length mismatch")
labels = np.asarray([1 if "tumour" in x.lower() or "tumor" in x.lower() else 0 for x in source], dtype=np.int8)
sample_meta = pd.DataFrame({"sample_id": sample_ids, "title": titles, "source_name": source, "label": labels})
sample_meta.to_csv(OUT / "sample_metadata.tsv", sep="\t", index=False)
pd.DataFrame({"probe_id": wanted, "missing_from_external": [x in missing_features for x in wanted]}).to_csv(
    OUT / "feature_coverage.tsv", sep="\t", index=False)
with h5py.File(OUT / "gse52865_top1000.h5", "w") as h:
    h.create_dataset("beta", data=matrix, compression="gzip")
    h.create_dataset("label", data=labels)
    h.create_dataset("probe_id", data=np.asarray(wanted, dtype="S"))
    h.create_dataset("sample_id", data=np.asarray(sample_ids, dtype="S"))

summary = {
    "accession": "GSE52865", "source_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52865/matrix/GSE52865_series_matrix.txt.gz",
    "source_sha256": sha256.hexdigest(), "n_samples": len(sample_ids), "n_tumor": int(labels.sum()),
    "n_normal": int((labels == 0).sum()), "requested_features": len(wanted),
    "features_found": len(values), "features_missing": len(missing_features),
    "matrix_missing_rate": float(np.isnan(matrix).mean()), "beta_min": float(np.nanmin(matrix)),
    "beta_max": float(np.nanmax(matrix)),
}
(OUT / "preparation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
