suppressPackageStartupMessages({
  library(data.table)
  library(missMethyl)
  library(jsonlite)
})
root <- Sys.getenv("METHYLATION_ROOT", unset = normalizePath(".", mustWork = TRUE))
out <- file.path(root, "results/enrichment")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
full <- fread(file.path(root, "results/dmp/dmp_full_cohort.tsv.gz"),
              select = c("probe_id", "adj.P.Val", "delta_beta"))
paired <- fread(file.path(root, "results/dmp/dmp_paired_difference.tsv.gz"),
                select = c("probe_id", "adj.P.Val", "mean_delta_beta"))
x <- merge(full, paired, by = "probe_id")
x[, common_large_effect := adj.P.Val.x < 0.05 & abs(delta_beta) >= 0.20 &
                           adj.P.Val.y < 0.05 & abs(mean_delta_beta) >= 0.20 &
                           sign(delta_beta) == sign(mean_delta_beta)]
sig <- x[common_large_effect == TRUE, probe_id]
all <- x$probe_id
fwrite(data.table(probe_id = sig), file.path(out, "common_large_effect_cpgs.tsv"), sep = "\t")

go <- gometh(sig.cpg = sig, all.cpg = all, collection = "GO",
             array.type = "450K", plot.bias = FALSE, prior.prob = TRUE)
go_dt <- as.data.table(go, keep.rownames = "term_id")
fwrite(go_dt, file.path(out, "GO_missMethyl.tsv.gz"), sep = "\t")
fwrite(head(go_dt[order(P.DE)], 50), file.path(out, "GO_top50.tsv"), sep = "\t")

kegg_error <- NULL
kegg_dt <- tryCatch({
  kegg <- gometh(sig.cpg = sig, all.cpg = all, collection = "KEGG",
                 array.type = "450K", plot.bias = FALSE, prior.prob = TRUE)
  z <- as.data.table(kegg, keep.rownames = "term_id")
  fwrite(z, file.path(out, "KEGG_missMethyl.tsv.gz"), sep = "\t")
  fwrite(head(z[order(P.DE)], 50), file.path(out, "KEGG_top50.tsv"), sep = "\t")
  z
}, error = function(e) {
  kegg_error <<- conditionMessage(e)
  NULL
})
summary <- list(n_background_cpgs = length(all), n_common_large_effect_cpgs = length(sig),
                go_terms_fdr_0_05 = sum(p.adjust(go_dt$P.DE, "BH") < 0.05, na.rm = TRUE),
                kegg_terms_fdr_0_05 = if (is.null(kegg_dt)) NA else sum(p.adjust(kegg_dt$P.DE, "BH") < 0.05, na.rm = TRUE),
                kegg_error = kegg_error,
                method = "missMethyl gometh; GO/KEGG; 450K probe-number bias correction")
write_json(summary, file.path(out, "enrichment_summary.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), "\n")

