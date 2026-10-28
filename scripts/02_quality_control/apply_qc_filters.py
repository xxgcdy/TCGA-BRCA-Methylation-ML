#!/usr/bin/env python3
import json
from datetime import datetime, timezone
import os
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
PROCESSED = PROJECT / "processed"
QC = PROJECT / "results" / "qc"

matrix = np.load(PROCESSED / "beta_matrix_float32.npy", mmap_mode="r")
sample_qc = pd.read_csv(QC / "sample_qc_unfiltered.tsv", sep="\t")
probe_order = pd.read_csv(PROCESSED / "probe_order.tsv.gz", sep="\t")
annotation = pd.read_csv(PROJECT / "annotation" / "hm450_annotation.tsv.gz", sep="\t", low_memory=False)

sample_qc["qc_keep"] = sample_qc["missing_rate"] <= 0.20
sample_qc["qc_reason"] = np.where(sample_qc["qc_keep"], "PASS", "sample_missing_rate_gt_0.20")
sample_qc.to_csv(QC / "sample_filter_decisions.tsv", sep="\t", index=False)
sample_keep = sample_qc["qc_keep"].to_numpy(dtype=bool)
sample_indices = np.flatnonzero(sample_keep)

n_kept_samples = int(sample_keep.sum())
n_probes = matrix.shape[1]
missing_count = np.zeros(n_probes, dtype=np.int32)
chunk = 25000
for start in range(0, n_probes, chunk):
    end = min(start + chunk, n_probes)
    missing_count[start:end] = np.isnan(matrix[sample_indices, start:end]).sum(axis=0)

probe = probe_order.copy()
probe["n_missing_after_sample_qc"] = missing_count
probe["missing_rate_after_sample_qc"] = missing_count / n_kept_samples
anno_cols = [
    "probe_id", "chr", "pos", "strand", "Name", "Type", "Probe_rs", "Probe_maf",
    "CpG_rs", "CpG_maf", "SBE_rs", "SBE_maf", "Islands_Name", "Relation_to_Island",
    "Random_Loci", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group",
    "DMR", "Enhancer", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DHS"
]
probe = probe.merge(annotation[anno_cols], on="probe_id", how="left", validate="one_to_one")

probe["annotated"] = probe["chr"].notna()
probe["missing_pass"] = probe["missing_rate_after_sample_qc"] <= 0.20
probe["random_locus"] = probe["Random_Loci"].astype(str).str.upper().eq("TRUE")
probe["common_cpg_snp"] = pd.to_numeric(probe["CpG_maf"], errors="coerce").ge(0.05).fillna(False)
probe["common_sbe_snp"] = pd.to_numeric(probe["SBE_maf"], errors="coerce").ge(0.05).fillna(False)
probe["autosomal"] = probe["chr"].isin([f"chr{i}" for i in range(1, 23)])
probe["qc_keep_all_chromosomes"] = (
    probe["annotated"]
    & probe["missing_pass"]
    & ~probe["random_locus"]
    & ~probe["common_cpg_snp"]
    & ~probe["common_sbe_snp"]
)
probe["qc_keep_autosomal"] = probe["qc_keep_all_chromosomes"] & probe["autosomal"]

probe.to_csv(QC / "probe_filter_decisions.tsv.gz", sep="\t", index=False, compression="gzip")
np.save(PROCESSED / "sample_keep_qc.npy", sample_keep)
np.save(PROCESSED / "probe_keep_qc_all_chromosomes.npy", probe["qc_keep_all_chromosomes"].to_numpy(dtype=bool))
np.save(PROCESSED / "probe_keep_qc_autosomal.npy", probe["qc_keep_autosomal"].to_numpy(dtype=bool))

kept_samples = sample_qc[sample_qc["qc_keep"]].copy()
kept_samples.to_csv(PROCESSED / "samples_after_qc.tsv", sep="\t", index=False)
probe.loc[probe["qc_keep_all_chromosomes"], ["matrix_column"] + anno_cols].to_csv(
    PROCESSED / "probes_after_qc_all_chromosomes.tsv.gz", sep="\t", index=False, compression="gzip"
)
probe.loc[probe["qc_keep_autosomal"], ["matrix_column"] + anno_cols].to_csv(
    PROCESSED / "probes_after_qc_autosomal.tsv.gz", sep="\t", index=False, compression="gzip"
)

summary = {
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "input_samples": int(len(sample_qc)),
    "retained_samples": n_kept_samples,
    "excluded_samples": int((~sample_keep).sum()),
    "retained_primary_tumor": int(((sample_qc["sample_type"] == "Primary Tumor") & sample_qc["qc_keep"]).sum()),
    "retained_solid_tissue_normal": int(((sample_qc["sample_type"] == "Solid Tissue Normal") & sample_qc["qc_keep"]).sum()),
    "input_probes": int(len(probe)),
    "unannotated_probes": int((~probe["annotated"]).sum()),
    "probes_failing_missing_rate": int((~probe["missing_pass"]).sum()),
    "random_loci": int(probe["random_locus"].sum()),
    "common_cpg_snp": int(probe["common_cpg_snp"].sum()),
    "common_sbe_snp": int(probe["common_sbe_snp"].sum()),
    "retained_probes_all_chromosomes": int(probe["qc_keep_all_chromosomes"].sum()),
    "retained_probes_autosomal": int(probe["qc_keep_autosomal"].sum()),
}
(QC / "qc_filter_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2))

