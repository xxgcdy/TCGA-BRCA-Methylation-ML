# Data dictionary

The repository does not contain TCGA sample-level beta matrices. The following files are expected under `METHYLATION_ROOT` when reproducing the workflow.

| Path | Description |
|---|---|
| `raw_beta/` | GDC SeSAMe Level-3 HM450 beta files, one directory per GDC file ID |
| `metadata/final_sample_sheet.tsv` | One retained aliquot per biological sample after technical-replicate handling |
| `metadata/sample_split.tsv` | Patient-level train/validation/test assignment and binary label |
| `annotation/hm450_annotation.tsv.gz` | HM450 hg19 annotation exported from the Bioconductor annotation package |
| `processed/beta_matrix_float32.npy` | Initial sample-by-probe beta matrix |
| `processed/beta_qc_autosomal.h5` | QC-filtered sample-by-autosomal-CpG beta matrix |
| `processed/probe_order_qc_autosomal.tsv.gz` | Probe order, chromosome, position, gene and island annotation |
| `processed/sample_order_qc_autosomal.tsv` | Sample order corresponding to the HDF5 matrix |

Important model outputs and summary tables committed to this repository are documented in the README and `docs/analysis_report.md`.
