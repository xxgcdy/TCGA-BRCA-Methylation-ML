#!/usr/bin/env python3
import json
from datetime import datetime, timezone
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

PROJECT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
PROCESSED = PROJECT / "processed"
matrix = np.load(PROCESSED / "beta_matrix_float32.npy", mmap_mode="r")
sample_keep = np.load(PROCESSED / "sample_keep_qc.npy")
probe_keep = np.load(PROCESSED / "probe_keep_qc_autosomal.npy")
sample_indices = np.flatnonzero(sample_keep)
probe_indices = np.flatnonzero(probe_keep)
samples = pd.read_csv(PROCESSED / "sample_order.tsv", sep="\t").iloc[sample_indices].reset_index(drop=True)
probes = pd.read_csv(PROCESSED / "probes_after_qc_autosomal.tsv.gz", sep="\t", low_memory=False)

if len(samples) != len(sample_indices) or len(probes) != len(probe_indices):
    raise RuntimeError("QC mask and metadata size mismatch")

path = PROCESSED / "beta_qc_autosomal.h5"
missing_total = 0
with h5py.File(path, "w") as handle:
    dataset = handle.create_dataset(
        "beta",
        shape=(len(sample_indices), len(probe_indices)),
        dtype="float32",
        chunks=(8, 4096),
        compression="lzf",
        shuffle=True,
        fillvalue=np.nan,
    )
    for start in range(0, len(sample_indices), 8):
        end = min(start + 8, len(sample_indices))
        block = np.asarray(matrix[np.ix_(sample_indices[start:end], probe_indices)], dtype=np.float32)
        missing_total += int(np.isnan(block).sum())
        dataset[start:end, :] = block
        if end % 80 == 0 or end == len(sample_indices):
            print(f"progress={end}/{len(sample_indices)}", flush=True)
    handle.create_dataset("source_sample_matrix_row", data=sample_indices.astype(np.int32))
    handle.create_dataset("source_probe_matrix_column", data=probe_indices.astype(np.int32))
    handle.attrs["orientation"] = "samples_by_probes"
    handle.attrs["sample_metadata"] = "sample_order_qc_autosomal.tsv"
    handle.attrs["probe_metadata"] = "probe_order_qc_autosomal.tsv.gz"

samples.to_csv(PROCESSED / "sample_order_qc_autosomal.tsv", sep="\t", index=False)
probes.to_csv(PROCESSED / "probe_order_qc_autosomal.tsv.gz", sep="\t", index=False, compression="gzip")
summary = {
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "path": str(path),
    "shape_samples_by_probes": [int(len(sample_indices)), int(len(probe_indices))],
    "dtype": "float32",
    "compression": "lzf",
    "file_bytes": int(path.stat().st_size),
    "missing_values": missing_total,
    "missing_rate": float(missing_total / (len(sample_indices) * len(probe_indices))),
}
(PROCESSED / "beta_qc_autosomal_metadata.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2))

