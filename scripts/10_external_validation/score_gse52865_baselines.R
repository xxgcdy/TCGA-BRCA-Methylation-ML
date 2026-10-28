suppressPackageStartupMessages({
  library(rhdf5); library(data.table); library(matrixStats); library(glmnet)
  library(ranger); library(xgboost); library(pROC); library(jsonlite)
})

root <- Sys.getenv("METHYLATION_ROOT", unset = normalizePath(".", mustWork = TRUE))
ext_dir <- file.path(root, "external/GSE52865/processed")
out <- file.path(root, "results/external_validation/GSE52865")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

ext <- h5read(file.path(ext_dir, "gse52865_top1000.h5"), "beta")
y <- as.integer(h5read(file.path(ext_dir, "gse52865_top1000.h5"), "label"))
sample_ids <- as.character(h5read(file.path(ext_dir, "gse52865_top1000.h5"), "sample_id"))
features <- fread(file.path(root, "results/models/training_only_top1000_features.tsv.gz"))
beta <- h5read(file.path(root, "processed/beta_qc_autosomal.h5"), "beta")
split_meta <- fread(file.path(root, "metadata/sample_split.tsv"))
samples <- fread(file.path(root, "processed/sample_order_qc_autosomal.tsv"))
meta <- merge(samples[, .(sample_id)], split_meta[, .(sample_id, split)], by = "sample_id", sort = FALSE)
meta <- meta[match(samples$sample_id, sample_id)]
train_cols <- which(meta$split == "train")
train_median <- rowMedians(beta[as.integer(features$probe_idx), train_cols, drop = FALSE], na.rm = TRUE)

# Python wrote samples x features; rhdf5 exposes the transposed orientation.
if (nrow(ext) == 1000 && ncol(ext) == length(y)) ext <- t(ext)
stopifnot(nrow(ext) == length(y), ncol(ext) == 1000)
colnames(ext) <- features$probe_id
miss <- which(is.na(ext), arr.ind = TRUE)
if (nrow(miss)) ext[miss] <- train_median[miss[, 2]]

metric <- function(model, p, threshold) {
  pred <- as.integer(p >= threshold)
  data.table(model, n = length(y), n_tumor = sum(y == 1), n_normal = sum(y == 0), threshold,
    roc_auc = as.numeric(auc(roc(y, p, quiet = TRUE, direction = "<"))),
    accuracy = mean(pred == y), sensitivity = sum(pred == 1 & y == 1) / sum(y == 1),
    specificity = sum(pred == 0 & y == 0) / sum(y == 0))
}

baseline_metrics <- fread(file.path(root, "results/models/baseline_metrics.tsv"))
threshold_for <- function(name) baseline_metrics[model == name & split == "test", threshold][1]

lasso <- readRDS(file.path(root, "results/models/lasso_glmnet.rds"))
rf <- readRDS(file.path(root, "results/models/ranger_random_forest.rds"))
reliable_xgb <- xgb.load(file.path(root, "results/models/xgboost_grouped_cv/xgboost_grouped_cv.ubj"))
xgb_metrics <- fread(file.path(root, "results/models/xgboost_grouped_cv/metrics.tsv"))

predictions <- list(
  LASSO = as.numeric(predict(lasso, ext, s = "lambda.1se", type = "response")),
  RandomForest = predict(rf, data.frame(ext, check.names = FALSE), num.threads = 8)$predictions[, "1"],
  XGBoostGroupedCV = predict(reliable_xgb, ext)
)
thresholds <- c(LASSO = threshold_for("LASSO"), RandomForest = threshold_for("RandomForest"),
                XGBoostGroupedCV = xgb_metrics[split == "test", threshold][1])
metrics <- rbindlist(lapply(names(predictions), function(name) metric(name, predictions[[name]], thresholds[[name]])))
pred_dt <- rbindlist(lapply(names(predictions), function(name) data.table(
  model = name, sample_id = sample_ids, label = y, probability = predictions[[name]], threshold = thresholds[[name]])))
fwrite(metrics, file.path(out, "baseline_external_metrics.tsv"), sep = "\t")
fwrite(pred_dt, file.path(out, "baseline_external_predictions.tsv.gz"), sep = "\t")
write_json(list(cohort = "GSE52865", frozen_features = 1000, threshold_source = "TCGA validation",
                metrics = metrics), file.path(out, "baseline_external_summary.json"), pretty = TRUE, auto_unbox = TRUE)
print(metrics)
