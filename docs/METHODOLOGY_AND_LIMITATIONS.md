# Methodology, quality controls, and limitations

This document records the methodological checks used to assess the reliability of the TCGA-BRCA HM450 analysis. It separates verified evidence from interpretation and preserves known technical limitations.

## Data provenance and cohort construction

- Source: open-access TCGA-BRCA Level-3 SeSAMe methylation beta estimates from the NCI Genomic Data Commons.
- All 890 downloaded files passed MD5 verification.
- Technical replicate handling preceded quality control; the lower-missingness aliquot was retained for duplicated biological samples.
- The final autosomal matrix contains 869 samples (773 primary tumors and 96 solid-tissue normals) and 390,644 CpGs.
- The paired analysis contains 90 patient-matched tumor-normal pairs.

Level-3 data do not include raw IDAT intensities. Detection P-values, dye-bias correction, and raw-signal quality control therefore cannot be recomputed in this project.

## Methylation quality control and inference

- Samples and probes with more than 20% missing beta values were excluded.
- Random loci, common CpG/SBE SNP probes, and non-autosomal probes were excluded from the main analysis.
- Differential methylation inference used M-values; beta values were used to report interpretable effect sizes.
- Both full-cohort and patient-paired DMP analyses were performed. Agreement between them is more informative than the unpaired P-values alone because the full cohort is imbalanced and may contain residual confounding.
- DMRs labelled “high confidence” satisfy explicit overlap, direction, CpG-count, and effect-size rules. They are not claimed to be one-to-one consensus regions.
- GO and KEGG enrichment used all QC-passed CpGs as the background and `missMethyl` probe-number bias correction. CpG-to-gene proximity does not establish regulation, causality, pathway activation, or pathway repression.

The workflow does not currently apply a separate cross-reactive-probe blacklist or explicitly adjust for tumor purity, cell composition, or acquisition-center batch effects.

## Leakage controls

- Data were split at the patient level before feature selection. Train, validation, and test sets contain 604, 131, and 134 samples, with no patient overlap.
- The 1,000 modeling CpGs were selected with limma using the training split only.
- Imputation medians and standardization parameters were estimated from training data only.
- Validation data selected classification thresholds and Transformer stopping epochs.
- The held-out TCGA test split was used only for final internal evaluation.
- GSE52865 was evaluated with frozen CpGs, preprocessing statistics, model parameters, and TCGA validation thresholds. External labels were not used for feature selection, retraining, or threshold tuning.

LASSO cross-validation occurs after feature selection on the complete training split. Its internal CV score is therefore not a fully nested estimate, although the held-out validation, test, and external evaluations remain isolated from feature selection.

## Model reliability checks

### Conventional models

The original XGBoost artifact had `best_iteration=0` and no retained evaluation log. It was replaced by five-fold patient-grouped cross-validation within the training set. Log-loss selected `max_depth=2`, `eta=0.02`, and 382 rounds. The corrected frozen test ROC-AUC is 0.9958.

LASSO and random forest achieve test ROC-AUC values of 0.9926 and 0.9947. These very high values should be interpreted in the context of a strong tissue tumor-normal signal rather than as clinical diagnostic performance.

### Transformer

The corrected Transformer constructs 10-CpG patches independently within each chromosome. It contains 111 patches, including 22 padded chromosome-end patches, and zero cross-chromosome patches. Five random seeds were evaluated; the probability ensemble has a test ROC-AUC of 0.9873.

The supervised input panel is still based on 1,000 training-selected CpGs. This makes the model suitable for a controlled classifier comparison but limits conclusions about genome-scale representation learning.

### Masked methylation pretraining

Masked-patch reconstruction used only the 604-sample TCGA training split. Three positional modes, random versus pretrained initialization, four label fractions, and three seeds produced 72 fine-tuning runs. Under genomic positional encoding, mean test ROC-AUC gains from pretraining were small (+0.0009 to +0.0032); other positional settings included neutral or negative changes.

The experiment demonstrates a reproducible pretraining pipeline, not a general advantage of pretraining. The unlabeled corpus is small and the token set remains supervised-selected.

## External validation

GSE52865 contains 40 breast tumors and 17 normal breast tissues. All 1,000 frozen CpGs were present and the extracted matrix had no missing values. LASSO, random forest, and corrected XGBoost each achieved an external ROC-AUC of 1.0000. The genomic Transformer achieved 0.9779 from random initialization and 0.9941 after masked pretraining.

Frozen-threshold sensitivity and specificity reveal calibration shift despite high ROC-AUC. For example, the pretrained genomic Transformer has sensitivity 0.625 and specificity 1.000 at its TCGA validation threshold. Ranking performance and deployment calibration are therefore reported separately.

The external cohort is small and differs from TCGA in cohort composition and preprocessing. One external dataset is not sufficient to establish clinical or cross-platform generalizability.

## Statistical uncertainty and interpretation boundaries

- The internal test set contains only 16 normal samples; one normal error changes specificity by 6.25 percentage points.
- A single patient-level split cannot characterize split-to-split variation. Repeated nested evaluation remains desirable.
- Bootstrap intervals quantify sampling uncertainty within the available test set but cannot replace independent cohorts.
- Tissue methylation classification is not evidence for a blood-based screening assay.
- Results are research findings and are not a validated clinical diagnostic system.

## Reproducible evidence

Machine-readable summaries are stored under:

- `results/evaluation/`
- `results/model_metrics/`
- `results/models/xgboost_grouped_cv/`
- `results/transformer_chromosome_safe/`
- `results/masked_pretraining/`
- `results/external_validation/GSE52865/`

Execution order, software environments, and expected inputs are documented in `docs/reproducibility.md` and `docs/data_dictionary.md`.
