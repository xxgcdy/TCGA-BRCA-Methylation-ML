suppressPackageStartupMessages({
  library(data.table)
  library(missMethyl)
  library(jsonlite)
})

root <- Sys.getenv("METHYLATION_ROOT", unset = normalizePath(".", mustWork = TRUE))
out <- file.path(root, "results/enrichment")
fig <- file.path(root, "results/figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)

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

links <- fread(file.path(root, "annotation/kegg/hsa_pathway_links.tsv"),
               header = FALSE, col.names = c("gene_id", "pathway_id"))
paths <- fread(file.path(root, "annotation/kegg/hsa_pathways.tsv"),
               header = FALSE, col.names = c("pathway_id", "pathway_name"))
links[, gene_id := sub("^hsa:", "", gene_id)]
links[, pathway_id := sub("^path:", "", pathway_id)]
collection <- split(links$gene_id, links$pathway_id)

kegg <- gsameth(sig.cpg = sig, all.cpg = all, collection = collection,
                array.type = "450K", plot.bias = FALSE, prior.prob = TRUE)
kegg_dt <- as.data.table(kegg, keep.rownames = "pathway_id")
kegg_dt <- merge(kegg_dt, paths, by = "pathway_id", all.x = TRUE, sort = FALSE)
setorder(kegg_dt, P.DE)
fwrite(kegg_dt, file.path(out, "KEGG_missMethyl.tsv.gz"), sep = "\t")
fwrite(head(kegg_dt, 50), file.path(out, "KEGG_top50.tsv"), sep = "\t")

summary_path <- file.path(out, "enrichment_summary.json")
summary <- if (file.exists(summary_path)) fromJSON(summary_path, simplifyVector = FALSE) else list()
summary$n_background_cpgs <- length(all)
summary$n_common_large_effect_cpgs <- length(sig)
summary$kegg_terms_tested <- nrow(kegg_dt)
summary$kegg_terms_fdr_0_05 <- sum(kegg_dt$FDR < 0.05, na.rm = TRUE)
summary$kegg_error <- NULL
summary$kegg_method <- "missMethyl gsameth; locally cached KEGG REST mapping; 450K probe-number bias correction"
summary$kegg_mapping_gene_pathway_links <- nrow(links)
summary$kegg_mapping_pathways <- nrow(paths)
write_json(summary, summary_path, pretty = TRUE, auto_unbox = TRUE, null = "null")

top <- head(kegg_dt, 20)
top[, label := sub(" - Homo sapiens.*$", "", pathway_name)]
top[, label := ifelse(nchar(label) > 55, paste0(substr(label, 1, 52), "..."), label)]
png(file.path(fig, "KEGG_top20.png"), 1800, 1500, res = 180)
par(mar = c(5, 19, 2, 2))
barplot(rev(-log10(top$FDR)), names.arg = rev(top$label), horiz = TRUE, las = 1,
        col = "#009E73", border = NA, xlab = expression(-log[10](FDR)),
        cex.names = 0.75, main = "Top KEGG enrichment (missMethyl bias-corrected)")
dev.off()

cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE, null = "null"), "\n")
