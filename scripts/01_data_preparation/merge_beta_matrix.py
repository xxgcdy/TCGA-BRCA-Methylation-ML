#!/usr/bin/env python3
import json
import time
from datetime import datetime, timezone
import os
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
SHEET = PROJECT / "metadata" / "final_sample_sheet.tsv"
OUT = PROJECT / "processed"
QC = PROJECT / "results" / "qc"
LOG = PROJECT / "logs" / "merge_beta_matrix.progress.log"
OUT.mkdir(parents=True, exist_ok=True)
QC.mkdir(parents=True, exist_ok=True)


def log(message: str) -> None:
    line = f"{datetime.now().astimezone().isoformat()}\t{message}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def read_file(path: str):
    frame = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["probe_id", "beta"],
        dtype={"probe_id": "string", "beta": "float32"},
        na_values=["NA"],
    )
    return frame["probe_id"].to_numpy(dtype=str), frame["beta"].to_numpy(dtype=np.float32)


samples = pd.read_csv(SHEET, sep="\t", dtype=str)
samples["label"] = pd.to_numeric(samples["label"], errors="raise").astype("int8")
n_samples = len(samples)
first_ids, first_beta = read_file(samples.iloc[0]["absolute_path"])
n_probes = len(first_ids)
if len(np.unique(first_ids)) != n_probes:
    raise RuntimeError("Duplicate probe IDs in first file")

matrix_path = OUT / "beta_matrix_float32.npy"
matrix = np.lib.format.open_memmap(matrix_path, mode="w+", dtype=np.float32, shape=(n_samples, n_probes))
valid_count = np.zeros(n_probes, dtype=np.uint16)
sum_beta = np.zeros(n_probes, dtype=np.float64)
sum_sq_beta = np.zeros(n_probes, dtype=np.float64)
sample_qc = []
started = time.time()
LOG.write_text("", encoding="utf-8")
log(f"start samples={n_samples} probes={n_probes}")

for index, row in samples.iterrows():
    if index == 0:
        probe_ids, beta = first_ids, first_beta
    else:
        probe_ids, beta = read_file(row["absolute_path"])
    if len(probe_ids) != n_probes or not np.array_equal(probe_ids, first_ids):
        raise RuntimeError(f"Probe order mismatch: {row['id']} {row['filename']}")
    matrix[index, :] = beta
    valid = ~np.isnan(beta)
    values = beta[valid].astype(np.float64)
    valid_count[valid] += 1
    sum_beta[valid] += values
    sum_sq_beta[valid] += values * values
    if np.any((values < 0) | (values > 1)):
        raise RuntimeError(f"Out-of-range beta values: {row['id']}")
    sample_qc.append({
        "matrix_row": index,
        "file_id": row["id"],
        "case_id": row["case_id"],
        "case_submitter_id": row["case_submitter_id"],
        "sample_id": row["sample_id"],
        "sample_submitter_id": row["sample_submitter_id"],
        "sample_type": row["sample_type"],
        "label": int(row["label"]),
        "n_probes": n_probes,
        "n_missing": int((~valid).sum()),
        "missing_rate": float((~valid).mean()),
        "mean_beta": float(np.mean(values)),
        "sd_beta": float(np.std(values, ddof=1)),
        "q01_beta": float(np.quantile(values, 0.01)),
        "q25_beta": float(np.quantile(values, 0.25)),
        "median_beta": float(np.quantile(values, 0.50)),
        "q75_beta": float(np.quantile(values, 0.75)),
        "q99_beta": float(np.quantile(values, 0.99)),
    })
    if (index + 1) % 25 == 0 or index + 1 == n_samples:
        elapsed = time.time() - started
        rate = (index + 1) / elapsed
        eta = (n_samples - index - 1) / rate if rate else 0
        matrix.flush()
        log(f"progress={index + 1}/{n_samples} rate={rate:.2f}_samples_per_s eta_seconds={eta:.0f}")

matrix.flush()
del matrix

missing_count = n_samples - valid_count.astype(np.int32)
mean = np.divide(sum_beta, valid_count, out=np.full(n_probes, np.nan), where=valid_count > 0)
numerator = sum_sq_beta - np.divide(sum_beta * sum_beta, valid_count, out=np.zeros(n_probes), where=valid_count > 0)
variance = np.divide(numerator, valid_count - 1, out=np.full(n_probes, np.nan), where=valid_count > 1)
variance = np.maximum(variance, 0, where=~np.isnan(variance), out=variance)

probe_qc = pd.DataFrame({
    "matrix_column": np.arange(n_probes, dtype=np.int32),
    "probe_id": first_ids,
    "n_missing": missing_count,
    "missing_rate": missing_count / n_samples,
    "mean_beta": mean,
    "sd_beta": np.sqrt(variance),
})
probe_qc.to_csv(QC / "probe_qc_unfiltered.tsv.gz", sep="\t", index=False, compression="gzip")
pd.DataFrame(sample_qc).to_csv(QC / "sample_qc_unfiltered.tsv", sep="\t", index=False)
pd.DataFrame({"matrix_column": np.arange(n_probes), "probe_id": first_ids}).to_csv(
    OUT / "probe_order.tsv.gz", sep="\t", index=False, compression="gzip"
)
samples.assign(matrix_row=np.arange(n_samples)).to_csv(OUT / "sample_order.tsv", sep="\t", index=False)

summary = {
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "matrix_path": str(matrix_path),
    "matrix_shape_samples_by_probes": [n_samples, n_probes],
    "dtype": "float32",
    "bytes": int(matrix_path.stat().st_size),
    "total_missing_values": int(missing_count.sum()),
    "overall_missing_rate": float(missing_count.sum() / (n_samples * n_probes)),
    "samples_with_probe_order_mismatch": 0,
    "elapsed_seconds": float(time.time() - started),
}
(OUT / "beta_matrix_metadata.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
log("complete " + json.dumps(summary, separators=(",", ":")))

