suppressPackageStartupMessages({
  library(rhdf5)
  library(data.table)
  library(matrixStats)
  library(xgboost)
  library(pROC)
  library(jsonlite)
})

set.seed(20261028)
root <- Sys.getenv("METHYLATION_ROOT", unset = normalizePath(".", mustWork = TRUE))
out <- file.path(root, "results/models/xgboost_grouped_cv")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

samples <- fread(file.path(root, "processed/sample_order_qc_autosomal.tsv"))
probes <- fread(file.path(root, "processed/probe_order_qc_autosomal.tsv.gz"))
split_meta <- fread(file.path(root, "metadata/sample_split.tsv"))
features <- fread(file.path(root, "results/models/training_only_top1000_features.tsv.gz"))
beta <- h5read(file.path(root, "processed/beta_qc_autosomal.h5"), "beta")

meta <- merge(samples[, .(sample_id, case_id, sample_type)],
              split_meta[, .(sample_id, split, label)], by = "sample_id", sort = FALSE)
meta <- meta[match(samples$sample_id, sample_id)]
stopifnot(!anyNA(meta$split), nrow(beta) == nrow(probes), ncol(beta) == nrow(samples))

tr <- which(meta$split == "train")
va <- which(meta$split == "validation")
te <- which(meta$split == "test")
sel <- as.integer(features$probe_idx)

train_beta <- beta[sel, tr, drop = FALSE]
train_median <- rowMedians(train_beta, na.rm = TRUE)
make_x <- function(rows) {
  z <- beta[sel, rows, drop = FALSE]
  miss <- which(is.na(z), arr.ind = TRUE)
  if (nrow(miss)) z[miss] <- train_median[miss[, 1]]
  t(z)
}
xtr <- make_x(tr); xva <- make_x(va); xte <- make_x(te)
ytr <- meta$label[tr]; yva <- meta$label[va]; yte <- meta$label[te]

# Five deterministic patient-grouped folds; paired samples remain together.
cases <- unique(meta$case_id[tr])
case_labels <- vapply(cases, function(z) max(ytr[meta$case_id[tr] == z]), numeric(1))
fold_map <- integer(length(cases))
for (lab in sort(unique(case_labels))) {
  idx <- which(case_labels == lab)
  fold_map[idx] <- sample(rep(1:5, length.out = length(idx)))
}
names(fold_map) <- cases
fold_id <- unname(fold_map[meta$case_id[tr]])
stopifnot(all(vapply(cases, function(z) length(unique(fold_id[meta$case_id[tr] == z])) == 1, logical(1))))
folds <- lapply(1:5, function(k) which(fold_id == k))

weights <- ifelse(ytr == 1, length(ytr) / (2 * sum(ytr == 1)),
                  length(ytr) / (2 * sum(ytr == 0)))
dtr <- xgb.DMatrix(xtr, label = ytr, weight = weights)

grid <- CJ(max_depth = c(2L, 3L, 4L), eta = c(0.02, 0.05))
grid[, `:=`(min_child_weight = 3, subsample = 0.8, colsample_bytree = 0.6)]
cv_summaries <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  params <- list(objective = "binary:logistic", eval_metric = c("auc", "logloss"),
                 max_depth = grid$max_depth[i], eta = grid$eta[i],
                 min_child_weight = grid$min_child_weight[i],
                 subsample = grid$subsample[i], colsample_bytree = grid$colsample_bytree[i],
                 nthread = 8, seed = 20261028)
  cv <- xgb.cv(params = params, data = dtr, folds = folds, nrounds = 600,
               early_stopping_rounds = 50, maximize = FALSE, verbose = 0)
  log <- as.data.table(cv$evaluation_log)
  best <- which.min(log$test_logloss_mean)
  log[, grid_id := i]
  fwrite(log, file.path(out, sprintf("cv_grid_%02d.tsv", i)), sep = "\t")
  cv_summaries[[i]] <- data.table(grid_id = i, best_round = best,
      cv_logloss = log$test_logloss_mean[best], cv_logloss_sd = log$test_logloss_std[best],
      cv_auc = log$test_auc_mean[best], cv_auc_sd = log$test_auc_std[best])
}

cv_summary <- cbind(grid, rbindlist(cv_summaries)[, -"grid_id"])
setorder(cv_summary, cv_logloss, -cv_auc)
fwrite(cv_summary, file.path(out, "grouped_cv_grid_summary.tsv"), sep = "\t")
best <- cv_summary[1]

params <- list(objective = "binary:logistic", eval_metric = c("auc", "logloss"),
               max_depth = best$max_depth, eta = best$eta,
               min_child_weight = best$min_child_weight, subsample = best$subsample,
               colsample_bytree = best$colsample_bytree, nthread = 8, seed = 20261028)
model <- xgb.train(params = params, data = dtr, nrounds = best$best_round, verbose = 0)
xgb.save(model, file.path(out, "xgboost_grouped_cv.ubj"))

pva <- predict(model, xva); pte <- predict(model, xte)
roc_va <- roc(yva, pva, quiet = TRUE, direction = "<")
threshold <- as.numeric(coords(roc_va, "best", best.method = "youden", ret = "threshold", transpose = FALSE)[1])
average_precision <- function(y, p) {
  o <- order(p, decreasing = TRUE); yy <- y[o]
  sum(cumsum(yy) / seq_along(yy) * yy) / sum(yy)
}
metric_row <- function(split, y, p) {
  pred <- as.integer(p >= threshold)
  data.table(model = "XGBoostGroupedCV", split, n = length(y), n_tumor = sum(y == 1),
    n_normal = sum(y == 0), threshold,
    roc_auc = as.numeric(auc(roc(y, p, quiet = TRUE, direction = "<"))),
    average_precision = average_precision(y, p), accuracy = mean(pred == y),
    sensitivity = sum(pred == 1 & y == 1) / sum(y == 1),
    specificity = sum(pred == 0 & y == 0) / sum(y == 0))
}
metrics <- rbind(metric_row("validation", yva, pva), metric_row("test", yte, pte))
predictions <- rbind(
  data.table(model = "XGBoostGroupedCV", split = "validation", sample_id = meta$sample_id[va],
             case_id = meta$case_id[va], label = yva, probability = pva),
  data.table(model = "XGBoostGroupedCV", split = "test", sample_id = meta$sample_id[te],
             case_id = meta$case_id[te], label = yte, probability = pte))
fwrite(metrics, file.path(out, "metrics.tsv"), sep = "\t")
fwrite(predictions, file.path(out, "predictions.tsv.gz"), sep = "\t")
fwrite(data.table(case_id = names(fold_map), fold = as.integer(fold_map)),
       file.path(out, "patient_grouped_folds.tsv"), sep = "\t")
summary <- list(seed = 20261028, feature_set = "training-only limma top 1000",
  tuning = "5-fold patient-grouped CV on training set; logloss selects parameters and rounds",
  selected = as.list(best), validation_threshold = threshold, metrics = metrics)
write_json(summary, file.path(out, "summary.json"), pretty = TRUE, auto_unbox = TRUE)
print(cv_summary)
print(metrics)
