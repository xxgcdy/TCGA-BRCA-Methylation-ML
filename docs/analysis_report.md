# TCGA-BRCA HM450 DNA methylation analysis summary

## Dataset and QC

- Source: GDC TCGA-BRCA SeSAMe level-3 beta values, Illumina HumanMethylation450.
- Downloaded files: 890; MD5 passed: 890; failed/missing: 0.
- Technical replicates: two cases had duplicate aliquots; the lower-missingness aliquot was retained.
- Biological cohort before QC: 888 samples, 790 cases, 791 primary tumors, 97 solid normals, 91 paired cases.
- QC thresholds: sample missingness <=20%; probe missingness <=20%; remove random loci and common CpG/SBE SNP probes; main analysis restricted to autosomes.
- Final analysis matrix: 869 samples (773 tumors, 96 normals) x 390,644 CpGs.
- PCA outliers: 0 using within-sample-type robust z >6 across PCs 1-10.

## Leakage-safe patient split

- Train: 604 samples (538 tumor, 66 normal).
- Validation: 131 samples (117 tumor, 14 normal).
- Test: 134 samples (118 tumor, 16 normal).
- No patient appears in more than one split.

## Differential methylation

### CpG-level DMP

- Full cohort: 275,082 CpGs at FDR <0.05; 28,444 also had |delta beta| >=0.20.
- Paired (90 patient pairs): 244,021 CpGs at FDR <0.05; 28,827 also had |delta beta| >=0.20.
- Paired large-effect sites: 16,775 hypermethylated and 12,052 hypomethylated in tumor.
- Common large-effect sites: 25,975.
- Full versus paired delta-beta correlation over all probes: Pearson r=0.9924.

### Region-level DMR

- DMRcate full cohort: 35,070 DMRs.
- DMRcate paired analysis: 32,541 DMRs.
- Full regions with any paired overlap: 31,218.
- High-confidence full regions: 2,986, requiring overlap, consistent direction, >=3 CpGs in both analyses, |mean difference| >=0.10 and |max difference| >=0.20.
- Leading regions include LBXCOR1, DMRTA2, PAX6, SP9 and UNCX.

## Functional enrichment

- missMethyl `gometh` used 25,975 common large-effect CpGs against all 390,644 tested CpGs, correcting HM450 probe-number bias.
- GO terms at FDR <0.05: 548.
- KEGG used a locally cached official KEGG REST mapping with `missMethyl::gsameth`; 14 of 372 tested pathways passed FDR < 0.05.

## Prediction models

Feature selection used only the training split (limma top 1,000 CpGs). The validation split selected thresholds and early stopping. The test split was used only for final evaluation.

| Model | Test ROC-AUC | Average precision | Accuracy | Sensitivity | Specificity |
|---|---:|---:|---:|---:|---:|
| LASSO | 0.9926 | 0.9990 | 0.9701 | 0.9746 | 0.9375 |
| Random forest | 0.9947 | 0.9993 | 0.9701 | 0.9746 | 0.9375 |
| XGBoost (patient-grouped CV) | 0.9958 | 0.9994 | 0.9776 | 0.9831 | 0.9375 |
| Chromosome-safe Transformer (5-seed ensemble) | 0.9873 | 0.9983 | 0.9627 | 0.9576 | 1.0000 |

The corrected XGBoost uses five-fold patient-grouped cross-validation within the training split. The compact Transformer remains below the strongest conventional baselines on the internal test set. TCGA tumor/normal separation is strong and the number of normal test samples is limited.

## Transformer artifact

- Input: 869 samples x 1,000 training-selected CpGs ordered by chromosome and hg19 coordinate.
- Missing values before imputation: 2,092 (0.2407%); imputed using training-only per-token medians, with the original missing mask retained.
- Corrected tokenization: 111 chromosome-local patches of up to 10 CpGs; 22 chromosome-end patches are padded and no patch crosses a chromosome boundary.
- Architecture: dimension 64, 4 heads and 2 encoder layers; five random seeds were evaluated on CPU.
- Masked methylation pretraining and positional ablations were run locally on an RTX 5090 D. Gains were small and configuration-dependent.

## Independent external validation

Frozen models were evaluated on GSE52865 (40 breast tumors and 17 normal breast tissues). All 1,000 training-selected CpGs were present, with no missing values in the extracted matrix. LASSO, random forest and corrected XGBoost each achieved ROC-AUC 1.0000. The genomic Transformer achieved ROC-AUC 0.9779 from random initialization and 0.9941 after masked pretraining. Frozen-threshold sensitivity and specificity showed calibration shift, so ranking and threshold-dependent performance are reported separately.

## Main output locations

- QC: `results/qc/`
- DMP: `results/dmp/`
- DMR: `results/dmr/`
- Enrichment: `results/enrichment/`
- Baseline models: `results/models/`
- Transformer: `results/transformer/`
- Figures: `results/figures/`
- Transformer-ready data: `processed/transformer/`

## Caveats

- This is a tumor-versus-adjacent-normal proof-of-concept, not a clinical diagnostic validation.
- TCGA-BRCA contains substantial class imbalance and possible center/batch effects.
- GSE52865 is small and differs from TCGA in cohort composition and preprocessing; additional external cohorts and raw-IDAT-level batch-aware analyses are needed before claiming broad generalization.
- KEGG results depend on the cached mapping snapshot and should be regenerated when pathway definitions are updated.
