#!/usr/bin/env bash
set -euo pipefail

ROOT="${METHYLATION_ROOT:-$(pwd)}"
STAGE="${1:-all}"
export METHYLATION_ROOT="$ROOT"

run() { echo "+ $*"; "$@"; }

case "$STAGE" in
  prepare)
    run python scripts/01_data_preparation/verify_download.py
    run python scripts/01_data_preparation/prepare_samples.py
    run python scripts/01_data_preparation/merge_beta_matrix.py
    run Rscript scripts/01_data_preparation/export_hm450_annotation.R
    ;;
  qc)
    run python scripts/02_quality_control/apply_qc_filters.py
    run python scripts/02_quality_control/run_qc_pca.py
    run python scripts/01_data_preparation/make_patient_splits.py
    run python scripts/02_quality_control/build_qc_hdf5.py
    ;;
  biology)
    run Rscript scripts/03_dmp/run_dmp_limma.R
    run Rscript scripts/03_dmp/run_paired_dmp_limma.R
    run Rscript scripts/03_dmp/summarize_dmp.R
    run Rscript scripts/04_dmr/run_dmr_full.R
    run Rscript scripts/04_dmr/run_dmr_paired.R
    run Rscript scripts/04_dmr/summarize_dmr_overlap.R
    run Rscript scripts/05_enrichment/run_missmethyl_enrichment.R
    [[ -s "$ROOT/annotation/kegg/hsa_pathway_links.tsv" ]] || run bash scripts/05_enrichment/cache_kegg_mapping.sh
    run Rscript scripts/05_enrichment/run_kegg_local.R
    ;;
  models)
    run Rscript scripts/06_model_training/run_baseline_models.R
    run python scripts/07_transformer/build_transformer_dataset.py
    run python scripts/07_transformer/train_methylation_transformer.py
    ;;
  evaluate)
    run python scripts/08_evaluation/evaluate_predictions.py
    run python scripts/08_evaluation/validate_project_evidence.py
    run python scripts/08_evaluation/make_portfolio_figures.py
    ;;
  all)
    for s in prepare qc biology models evaluate; do run bash "$0" "$s"; done
    ;;
  *) echo "Usage: $0 {prepare|qc|biology|models|evaluate|all}" >&2; exit 2 ;;
esac
