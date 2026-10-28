#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(data.table))

project <- Sys.getenv("METHYLATION_ROOT", unset = normalizePath(".", mustWork = TRUE))
processed <- file.path(project, "processed")
outdir <- file.path(project, "results", "dmp")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Reading QC HDF5 matrix")
beta_pxs <- h5read(file.path(processed, "beta_qc_autosomal.h5"), "beta")
samples <- fread(file.path(processed, "sample_order_qc_autosomal.tsv"), data.table = FALSE)
probes <- fread(file.path(processed, "probe_order_qc_autosomal.tsv.gz"), data.table = FALSE)
stopifnot(nrow(beta_pxs) == nrow(probes), ncol(beta_pxs) == nrow(samples))

group <- factor(samples$sample_type, levels = c("Solid Tissue Normal", "Primary Tumor"))
mean_normal <- rowMeans(beta_pxs[, group == "Solid Tissue Normal", drop = FALSE], na.rm = TRUE)
mean_tumor <- rowMeans(beta_pxs[, group == "Primary Tumor", drop = FALSE], na.rm = TRUE)
delta_beta <- mean_tumor - mean_normal

eps <- 1e-6
message("Transforming beta values to M-values")
m_full <- log2(pmin(pmax(beta_pxs, eps), 1 - eps) / (1 - pmin(pmax(beta_pxs, eps), 1 - eps)))
rownames(m_full) <- probes$probe_id
design <- model.matrix(~0 + group)
colnames(design) <- c("Normal", "Tumor")
contrast <- makeContrasts(TumorVsNormal = Tumor - Normal, levels = design)
message("Fitting full-cohort limma model")
fit <- lmFit(m_full, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit, robust = TRUE, trend = TRUE)
full <- topTable(fit, number = Inf, sort.by = "P", adjust.method = "BH")
full_idx <- match(rownames(full), probes$probe_id)
stopifnot(!anyNA(full_idx))
full$matrix_column <- probes$matrix_column[full_idx]
full$probe_id <- probes$probe_id[full_idx]
full$mean_beta_normal <- mean_normal[full_idx]
full$mean_beta_tumor <- mean_tumor[full_idx]
full$delta_beta <- delta_beta[full_idx]
full$significant_effect <- full$adj.P.Val < 0.05 & abs(full$delta_beta) >= 0.20
full$direction <- ifelse(full$significant_effect, ifelse(full$delta_beta > 0, "hyper", "hypo"), "not_significant_effect")
full <- cbind(full, probes[full_idx, setdiff(names(probes), c("matrix_column", "probe_id")), drop = FALSE])
fwrite(full, file.path(outdir, "dmp_full_cohort.tsv.gz"), sep = "\t", compress = "gzip", na = "NA")

summary <- list(
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  probes_tested = nrow(probes),
  full_cohort = list(
    samples = nrow(samples),
    tumor = sum(group == "Primary Tumor"),
    normal = sum(group == "Solid Tissue Normal"),
    fdr_lt_0_05 = sum(full$adj.P.Val < 0.05, na.rm = TRUE),
    fdr_lt_0_05_abs_delta_beta_ge_0_20 = sum(full$significant_effect, na.rm = TRUE),
    hyper = sum(full$direction == "hyper", na.rm = TRUE),
    hypo = sum(full$direction == "hypo", na.rm = TRUE)
  )
)
jsonlite::write_json(summary, file.path(outdir, "dmp_full_summary.json"), pretty = TRUE, auto_unbox = TRUE)
print(summary)

