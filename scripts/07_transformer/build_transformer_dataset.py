#!/usr/bin/env python3
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
OUT = ROOT / "processed" / "transformer"
OUT.mkdir(parents=True, exist_ok=True)

samples = pd.read_csv(ROOT / "processed/sample_order_qc_autosomal.tsv", sep="\t")
splits = pd.read_csv(ROOT / "metadata/sample_split.tsv", sep="\t")
samples = samples[["sample_id", "case_id", "sample_type"]].merge(
    splits[["sample_id", "split", "label"]], on="sample_id", how="left", validate="one_to_one"
)
features = pd.read_csv(ROOT / "results/models/training_only_top1000_features.tsv.gz", sep="\t")
probes = pd.read_csv(ROOT / "processed/probe_order_qc_autosomal.tsv.gz", sep="\t")

# R probe_idx is one-based; Python HDF5 indices are zero-based.
features["probe_index_zero_based"] = features["probe_idx"].astype(int) - 1
ann_cols = [c for c in ["probe_id", "chr", "pos", "UCSC_RefGene_Name", "Relation_to_Island"] if c in probes.columns]
features = features.merge(probes[ann_cols], on="probe_id", how="left", validate="one_to_one")
features["chr_num"] = features["chr"].str.replace("chr", "", regex=False).astype(int)
features = features.sort_values(["chr_num", "pos", "probe_id"]).reset_index(drop=True)

with h5py.File(ROOT / "processed/beta_qc_autosomal.h5", "r") as src:
    ds = src["beta"]
    if ds.shape != (len(samples), len(probes)):
        raise RuntimeError(f"Unexpected beta shape {ds.shape}; expected {(len(samples), len(probes))}")
    requested = features["probe_index_zero_based"].to_numpy()
    read_order = np.argsort(requested)
    restore_genomic_order = np.argsort(read_order)
    x = ds[:, requested[read_order]].astype(np.float32)[:, restore_genomic_order]

missing_mask = np.isnan(x)
train_mask = samples["split"].eq("train").to_numpy()
train_median = np.nanmedian(x[train_mask], axis=0).astype(np.float32)
inds = np.where(missing_mask)
x[inds] = train_median[inds[1]]
train_mean = x[train_mask].mean(axis=0).astype(np.float32)
train_std = x[train_mask].std(axis=0).astype(np.float32)
train_std[train_std < 1e-6] = 1.0

split_code = samples["split"].map({"train": 0, "validation": 1, "test": 2}).astype(np.int8).to_numpy()
out_h5 = OUT / "tcga_brca_hm450_top1000.h5"
with h5py.File(out_h5, "w") as h:
    h.create_dataset("beta", data=x, compression="lzf")
    h.create_dataset("missing_mask", data=missing_mask, compression="lzf")
    h.create_dataset("label", data=samples["label"].astype(np.int8).to_numpy())
    h.create_dataset("split_code", data=split_code)
    h.create_dataset("chromosome", data=features["chr_num"].astype(np.int8).to_numpy())
    h.create_dataset("position", data=features["pos"].astype(np.int32).to_numpy())
    h.create_dataset("train_median", data=train_median)
    h.create_dataset("train_mean", data=train_mean)
    h.create_dataset("train_std", data=train_std)
    h.attrs["split_code_map"] = json.dumps({"train": 0, "validation": 1, "test": 2})
    h.attrs["feature_selection"] = "training-only limma top 1000; genome-coordinate token order"

samples.to_csv(OUT / "sample_metadata.tsv", sep="\t", index=False)
features.to_csv(OUT / "token_metadata.tsv.gz", sep="\t", index=False, compression="gzip")
summary = {
    "n_samples": int(len(samples)),
    "n_tokens": int(len(features)),
    "shape_samples_by_tokens": list(x.shape),
    "n_train": int((split_code == 0).sum()),
    "n_validation": int((split_code == 1).sum()),
    "n_test": int((split_code == 2).sum()),
    "missing_values_before_imputation": int(missing_mask.sum()),
    "missing_rate": float(missing_mask.mean()),
    "imputation": "per-token median learned from training split only",
    "token_order": "chromosome then hg19 position",
    "output_h5": str(out_h5),
}
(OUT / "dataset_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))

