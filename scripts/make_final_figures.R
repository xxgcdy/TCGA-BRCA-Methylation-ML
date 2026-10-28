suppressPackageStartupMessages({library(data.table); library(pROC)})
root <- Sys.getenv("METHYLATION_ROOT", unset = normalizePath(".", mustWork = TRUE))
fig <- file.path(root, "results/figures")
dir.create(fig, recursive = TRUE, showWarnings = FALSE)

# Test ROC curves for all models.
b <- fread(file.path(root, "results/models/baseline_predictions.tsv.gz"))[split == "test"]
t <- fread(file.path(root, "results/transformer/transformer_predictions.tsv.gz"))[split == "test"]
setnames(t, "model", "model")
pred <- rbind(b[, .(model, label, probability)], t[, .(model, label, probability)])
cols <- c(LASSO = "#0072B2", RandomForest = "#009E73", XGBoost = "#D55E00", CpGPatchTransformer = "#CC79A7")
png(file.path(fig, "test_roc_all_models.png"), 1600, 1300, res = 180)
par(mar = c(5, 5, 2, 2))
first <- TRUE; labs <- character()
for (m in names(cols)) {
  z <- pred[model == m]; r <- roc(z$label, z$probability, quiet = TRUE, direction = "<")
  if (first) { plot(r, col = cols[m], lwd = 3, legacy.axes = TRUE, main = "Held-out test ROC"); first <- FALSE }
  else lines(r, col = cols[m], lwd = 3)
  labs <- c(labs, sprintf("%s (AUC %.3f)", m, as.numeric(auc(r))))
}
abline(0, 1, lty = 2, col = "grey60")
legend("bottomright", legend = labs, col = cols, lwd = 3, bty = "n", cex = 0.85)
dev.off()

# Paired DMP volcano.
d <- fread(file.path(root, "results/dmp/dmp_paired_difference.tsv.gz"),
           select = c("mean_delta_beta", "adj.P.Val", "significant_fdr_delta20"))
y <- pmin(-log10(pmax(d$adj.P.Val, 1e-300)), 60)
cc <- ifelse(d$significant_fdr_delta20 & d$mean_delta_beta > 0, "#D55E00",
             ifelse(d$significant_fdr_delta20, "#0072B2", "grey75"))
png(file.path(fig, "paired_dmp_volcano.png"), 1600, 1300, res = 180)
plot(d$mean_delta_beta, y, pch = 16, cex = 0.25, col = adjustcolor(cc, 0.5),
     xlab = expression(Delta*beta~"(tumor - normal)"), ylab = expression(-log[10](FDR)),
     main = "Paired differential methylation (90 patient pairs)")
abline(v = c(-0.2, 0.2), lty = 2); abline(h = -log10(0.05), lty = 2)
legend("topright", c("Hyper |Δβ|≥0.2", "Hypo |Δβ|≥0.2", "Other"),
       col = c("#D55E00", "#0072B2", "grey75"), pch = 16, bty = "n")
dev.off()

# Top GO terms.
go <- fread(file.path(root, "results/enrichment/GO_top50.tsv"))[1:20]
go[, label := ifelse(nchar(TERM) > 55, paste0(substr(TERM, 1, 52), "..."), TERM)]
png(file.path(fig, "GO_top20.png"), 1800, 1500, res = 180)
par(mar = c(5, 19, 2, 2))
barplot(rev(-log10(go$FDR)), names.arg = rev(go$label), horiz = TRUE, las = 1,
        col = "#56B4E9", border = NA, xlab = expression(-log[10](FDR)), cex.names = 0.75,
        main = "Top GO enrichment (missMethyl bias-corrected)")
dev.off()
cat("Figures written to", fig, "\n")

