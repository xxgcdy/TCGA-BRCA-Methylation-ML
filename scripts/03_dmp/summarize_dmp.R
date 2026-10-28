suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})
root <- Sys.getenv("METHYLATION_ROOT", unset = normalizePath(".", mustWork = TRUE))
out <- file.path(root, "results/dmp")
full_path <- file.path(out, "dmp_full_cohort.tsv.gz")
paired_path <- file.path(out, "dmp_paired_difference.tsv.gz")

full <- fread(full_path)
paired <- fread(paired_path)

delta_col_full <- if ("mean_delta_beta" %in% names(full)) "mean_delta_beta" else if ("delta_beta" %in% names(full)) "delta_beta" else stop("No delta beta column in full result")

full_sig <- full$adj.P.Val < 0.05 & abs(full[[delta_col_full]]) >= 0.20
paired_sig <- paired$adj.P.Val < 0.05 & abs(paired$mean_delta_beta) >= 0.20

summary <- list(
  full_cohort = list(
    n_probes = nrow(full),
    significant_fdr_0_05 = sum(full$adj.P.Val < 0.05, na.rm = TRUE),
    significant_fdr_0_05_delta_beta_0_20 = sum(full_sig, na.rm = TRUE),
    hypermethylated_delta20 = sum(full_sig & full[[delta_col_full]] > 0, na.rm = TRUE),
    hypomethylated_delta20 = sum(full_sig & full[[delta_col_full]] < 0, na.rm = TRUE)
  ),
  paired = list(
    n_probes = nrow(paired),
    significant_fdr_0_05 = sum(paired$adj.P.Val < 0.05, na.rm = TRUE),
    significant_fdr_0_05_delta_beta_0_20 = sum(paired_sig, na.rm = TRUE),
    hypermethylated_delta20 = sum(paired_sig & paired$mean_delta_beta > 0, na.rm = TRUE),
    hypomethylated_delta20 = sum(paired_sig & paired$mean_delta_beta < 0, na.rm = TRUE)
  ),
  overlap = list(
    significant_delta20_both = length(intersect(full$probe_id[full_sig], paired$probe_id[paired_sig])),
    pearson_delta_beta_all_probes = cor(full[[delta_col_full]][match(paired$probe_id, full$probe_id)], paired$mean_delta_beta, use = "complete.obs")
  )
)

write_json(summary, file.path(out, "dmp_combined_summary.json"), pretty = TRUE, auto_unbox = TRUE)
fwrite(full[order(P.Value)][1:20], file.path(out, "top20_dmp_full_cohort.tsv"), sep = "\t")
fwrite(paired[order(P.Value)][1:20], file.path(out, "top20_dmp_paired.tsv"), sep = "\t")
cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), "\n")

