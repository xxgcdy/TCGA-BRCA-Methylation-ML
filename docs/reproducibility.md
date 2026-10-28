# Reproducibility guide

## Data access

The analysis used open-access TCGA-BRCA SeSAMe Methylation Beta Estimation files from the NCI Genomic Data Commons. Generate a current GDC manifest and sample metadata through the GDC API or portal. Sample-level data are deliberately excluded from this repository.

## Project root

Run scripts from the repository root or set an explicit project directory:

```bash
export METHYLATION_ROOT=/path/to/TCGA-BRCA-Methylation-ML
```

All published scripts read this environment variable. If it is unset, they use the current working directory.

## Suggested execution order

1. `scripts/01_data_preparation/prepare_samples.py`
2. `scripts/01_data_preparation/merge_beta_matrix.py`
3. `scripts/01_data_preparation/export_hm450_annotation.R`
4. `scripts/02_quality_control/apply_qc_filters.py`
5. `scripts/02_quality_control/run_qc_pca.py`
6. `scripts/01_data_preparation/make_patient_splits.py`
7. `scripts/02_quality_control/build_qc_hdf5.py`
8. `scripts/03_dmp/run_dmp_limma.R` and `run_paired_dmp_limma.R`
9. `scripts/03_dmp/summarize_dmp.R`
10. `scripts/04_dmr/run_dmr_full.R` and `run_dmr_paired.R`
11. `scripts/04_dmr/summarize_dmr_overlap.R`
12. `scripts/05_enrichment/run_missmethyl_enrichment.R`
13. Cache the official KEGG human pathway links under `annotation/kegg/`, then run `scripts/05_enrichment/run_kegg_local.R`
14. `scripts/06_model_training/run_baseline_models.R`
15. `scripts/07_transformer/build_transformer_dataset.py`
16. `scripts/07_transformer/train_methylation_transformer.py`
17. `scripts/make_final_figures.R`

## Leakage controls

- Splitting is performed at the patient level.
- The top 1,000 model CpGs are selected with limma using training samples only.
- Imputation statistics are estimated from training samples only.
- The validation split is used for classification thresholds and Transformer early stopping.
- The test split is evaluated only after model development.

## Compute

The reported run used CPU execution capped at eight threads. The compact Transformer does not require a GPU. Raw data plus intermediate matrices required approximately 16 GB in this implementation; selected GitHub outputs are much smaller.

## KEGG mapping and offline reproducibility

The official human gene-pathway links and pathway names were downloaded from the KEGG REST API and cached locally before enrichment. The cached mapping was supplied to `missMethyl::gsameth`, preserving HM450 probe-number bias correction while avoiding repeated network requests. Because KEGG pathway definitions are updated over time, regenerate the cache when reproducing the analysis with a newer pathway release.
