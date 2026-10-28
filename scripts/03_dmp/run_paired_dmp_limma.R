suppressPackageStartupMessages({
  library(rhdf5)
  library(data.table)
  library(limma)
  library(jsonlite)
})

root <- Sys.getenv("METHYLATION_ROOT", unset = normalizePath(".", mustWork = TRUE))
h5_path <- file.path(root, "processed/beta_qc_autosomal.h5")
sample_path <- file.path(root, "processed/sample_order_qc_autosomal.tsv")
probe_path <- file.path(root, "processed/probe_order_qc_autosomal.tsv.gz")
out_dir <- file.path(root, "results/dmp")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- fread(sample_path)
probes <- fread(probe_path)
beta <- h5read(h5_path, "beta")

# rhdf5 exposes this dataset as probes x samples.
stopifnot(nrow(beta) == nrow(probes), ncol(beta) == nrow(samples))

paired_cases <- samples[, .(
  n_tumor = sum(sample_type == "Primary Tumor"),
  n_normal = sum(sample_type == "Solid Tissue Normal")
), by = case_id][n_tumor > 0 & n_normal > 0, case_id]

pair_rows <- rbindlist(lapply(sort(paired_cases), function(case) {
  z <- samples[case_id == case]
  tumor <- z[sample_type == "Primary Tumor"][order(sample_id)][1]
  normal <- z[sample_type == "Solid Tissue Normal"][order(sample_id)][1]
  data.table(case_id = case,
             tumor_sample_id = tumor$sample_id,
             normal_sample_id = normal$sample_id,
             tumor_col = match(tumor$sample_id, samples$sample_id),
             normal_col = match(normal$sample_id, samples$sample_id))
}))

fwrite(pair_rows[, !c("tumor_col", "normal_col")],
       file.path(out_dir, "paired_samples_used.tsv"), sep = "\t")

eps <- 1e-6
tumor_beta <- beta[, pair_rows$tumor_col, drop = FALSE]
normal_beta <- beta[, pair_rows$normal_col, drop = FALSE]
delta_beta <- rowMeans(tumor_beta - normal_beta, na.rm = TRUE)

tumor_m <- log2((pmin(pmax(tumor_beta, eps), 1 - eps)) /
                (1 - pmin(pmax(tumor_beta, eps), 1 - eps)))
normal_m <- log2((pmin(pmax(normal_beta, eps), 1 - eps)) /
                 (1 - pmin(pmax(normal_beta, eps), 1 - eps)))
m_diff <- tumor_m - normal_m
rownames(m_diff) <- probes$probe_id

fit <- lmFit(m_diff, matrix(1, nrow = ncol(m_diff), ncol = 1))
fit <- eBayes(fit, robust = TRUE, trend = TRUE)
tt <- as.data.table(topTable(fit, coef = 1, number = Inf, sort.by = "P"), keep.rownames = "probe_id")
tt[, probe_idx := match(probe_id, probes$probe_id)]
tt[, mean_delta_beta := delta_beta[probe_idx]]
ann_cols <- intersect(c("probe_id", "chr", "pos", "gene", "relation_to_island"), names(probes))
tt <- merge(tt, probes[, ..ann_cols], by = "probe_id", all.x = TRUE, sort = FALSE)
tt[, significant_fdr_delta20 := adj.P.Val < 0.05 & abs(mean_delta_beta) >= 0.20]
setorder(tt, P.Value)
fwrite(tt, file.path(out_dir, "dmp_paired_difference.tsv.gz"), sep = "\t")

summary <- list(
  method = "paired tumor-minus-normal M-value differences, limma one-sample moderated test",
  n_pairs = nrow(pair_rows),
  n_probes = nrow(tt),
  significant_fdr_0_05 = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
  significant_fdr_0_05_delta_beta_0_20 = sum(tt$significant_fdr_delta20, na.rm = TRUE),
  hypermethylated_delta20 = sum(tt$significant_fdr_delta20 & tt$mean_delta_beta > 0, na.rm = TRUE),
  hypomethylated_delta20 = sum(tt$significant_fdr_delta20 & tt$mean_delta_beta < 0, na.rm = TRUE)
)
write_json(summary, file.path(out_dir, "dmp_paired_summary.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), "\n")

