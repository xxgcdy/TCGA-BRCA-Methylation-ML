# Self-supervised methylation pretraining design

## Research question and implementation status

The experiment asks whether masked methylation pretraining improves classification when only a fraction of labeled training samples is available, and whether the answer depends on positional encoding.

It was implemented with Python 3.10 and PyTorch 2.11.0+cu128 on an NVIDIA RTX 5090 D. CUDA availability and a GPU matrix multiplication were verified before training. Pretraining used only the TCGA training split; model checkpoints and sample-level predictions remain excluded from Git.

## Chromosome-safe tokenization and positional encoding

1. CpGs are sorted by hg19 coordinate within each chromosome.
2. Ten consecutive CpGs form a patch; chromosome-end patches are padded.
3. The implemented representation has 111 patches, including 22 partial patches, and zero cross-chromosome patches.
4. Three positional modes are compared: none, learned sequence position plus chromosome embedding, and chromosome embedding plus hg19 coordinate Fourier features.
5. The current input remains the training-selected 1,000-CpG panel; this supervised vocabulary is an explicit limitation.

## Self-supervised objective

The implemented task uses masked methylation modeling:

- Randomly mask 30% of patches and replace their embeddings with a learned mask token.
- Reconstruct continuous beta values using mean-squared error.
- Compute loss only over masked, non-padding CpG positions.
- Use only training patients for pretraining.

## Evaluation matrix

| 变量 | 设置 |
|---|---|
| Initialization | random / masked pretrained |
| Labeled training fraction | 10% / 25% / 50% / 100% |
| Positional representation | none / sequence / genomic coordinate |
| Random seeds | 3 |
| Evaluation | frozen TCGA test split + frozen GSE52865 cohort |

The comparison contains 72 fine-tuning runs. Validation ROC-AUC selects stopping epochs, and validation-derived thresholds remain frozen for test evaluation.

## Observed results

Under genomic positional encoding, pretraining changed mean internal test ROC-AUC by +0.0009, +0.0032, +0.0019, and +0.0028 at 10%, 25%, 50%, and 100% labeled data. Other positional settings included neutral or negative changes.

On GSE52865, the genomic model achieved ROC-AUC 0.9779 from random initialization and 0.9941 after pretraining. At the frozen TCGA threshold, the pretrained model had sensitivity 0.625 and specificity 1.000, demonstrating cross-cohort calibration shift.

These results do not establish a general benefit of pretraining. The unlabeled corpus is small and the token set is supervised-selected. A stronger follow-up would use a larger multi-cohort unlabeled corpus, a label-independent CpG vocabulary, more random seeds, and prespecified calibration analyses.

## Reproducible artifacts

- `scripts/09_pretraining/run_masked_pretraining_experiments.py`
- `scripts/09_pretraining/summarize_pretraining_results.py`
- `results/masked_pretraining/`
- `results/external_validation/GSE52865/`
