suppressPackageStartupMessages({
  library(rhdf5)
  library(data.table)
  library(matrixStats)
  library(limma)
  library(glmnet)
  library(ranger)
  library(xgboost)
  library(pROC)
  library(jsonlite)
})
set.seed(20261027)
root <- Sys.getenv("METHYLATION_ROOT", unset = normalizePath(".", mustWork = TRUE))
out <- file.path(root, "results/models")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

samples <- fread(file.path(root, "processed/sample_order_qc_autosomal.tsv"))
probes <- fread(file.path(root, "processed/probe_order_qc_autosomal.tsv.gz"))
split_meta <- fread(file.path(root, "metadata/sample_split.tsv"))
beta <- h5read(file.path(root, "processed/beta_qc_autosomal.h5"), "beta")
stopifnot(nrow(beta) == nrow(probes), ncol(beta) == nrow(samples))
meta <- merge(samples[, .(sample_id, case_id, sample_type)],
              split_meta[, .(sample_id, split, label)], by = "sample_id", sort = FALSE)
meta <- meta[match(samples$sample_id, sample_id)]
stopifnot(!anyNA(meta$split), all(meta$label %in% c(0, 1)))

tr <- which(meta$split == "train")
va <- which(meta$split == "validation")
te <- which(meta$split == "test")

# Training-only feature selection by limma on M values.
btr <- beta[, tr, drop = FALSE]
med <- rowMedians(btr, na.rm = TRUE)
miss <- which(is.na(btr), arr.ind = TRUE)
if (nrow(miss)) btr[miss] <- med[miss[, 1]]
eps <- 1e-6
mtr <- log2(pmin(pmax(btr, eps), 1 - eps) / (1 - pmin(pmax(btr, eps), 1 - eps)))
design <- model.matrix(~ factor(meta$label[tr], levels = c(0, 1)))
fit <- eBayes(lmFit(mtr, design), robust = TRUE, trend = TRUE)
fs <- as.data.table(topTable(fit, coef = 2, number = Inf, sort.by = "P"), keep.rownames = "probe_idx")
fs[, probe_idx := as.integer(probe_idx)]
fs[, probe_id := probes$probe_id[probe_idx]]
train_delta_all <- rowMeans(btr[, meta$label[tr] == 1, drop = FALSE]) -
                   rowMeans(btr[, meta$label[tr] == 0, drop = FALSE])
fs[, train_delta_beta := train_delta_all[probe_idx]]
selected <- fs[order(P.Value)][1:1000]
fwrite(selected, file.path(out, "training_only_top1000_features.tsv.gz"), sep = "\t")
sel <- selected$probe_idx

make_x <- function(cols) {
  z <- beta[sel, cols, drop = FALSE]
  for (i in seq_len(nrow(z))) z[i, is.na(z[i, ])] <- med[sel[i]]
  t(z)
}
xtr <- make_x(tr); xva <- make_x(va); xte <- make_x(te)
colnames(xtr) <- colnames(xva) <- colnames(xte) <- selected$probe_id
ytr <- meta$label[tr]; yva <- meta$label[va]; yte <- meta$label[te]

auc_value <- function(y, p) as.numeric(auc(roc(y, p, quiet = TRUE, direction = "<")))
average_precision <- function(y, p) {
  o <- order(p, decreasing = TRUE); yy <- y[o]
  sum(cumsum(yy) / seq_along(yy) * yy) / sum(yy)
}
best_threshold <- function(y, p) {
  r <- roc(y, p, quiet = TRUE, direction = "<")
  as.numeric(coords(r, "best", best.method = "youden", ret = "threshold", transpose = FALSE)[1])
}
metrics <- function(model, split, y, p, threshold) {
  pred <- as.integer(p >= threshold)
  data.table(model, split, n = length(y), n_tumor = sum(y == 1), n_normal = sum(y == 0),
             threshold, roc_auc = auc_value(y, p), average_precision = average_precision(y, p),
             accuracy = mean(pred == y), sensitivity = sum(pred == 1 & y == 1) / sum(y == 1),
             specificity = sum(pred == 0 & y == 0) / sum(y == 0))
}
pred_tables <- list(); metric_tables <- list()
add_results <- function(name, pva, pte) {
  th <- best_threshold(yva, pva)
  metric_tables[[name]] <<- rbind(metrics(name, "validation", yva, pva, th),
                                  metrics(name, "test", yte, pte, th))
  pred_tables[[name]] <<- rbind(
    data.table(model = name, split = "validation", sample_id = meta$sample_id[va], case_id = meta$case_id[va], label = yva, probability = pva),
    data.table(model = name, split = "test", sample_id = meta$sample_id[te], case_id = meta$case_id[te], label = yte, probability = pte))
}

# Patient-grouped folds keep tumor/normal from a case in the same CV fold.
cases <- unique(meta$case_id[tr])
case_fold <- setNames(sample(rep(1:5, length.out = length(cases))), cases)
foldid <- unname(case_fold[meta$case_id[tr]])
weights <- ifelse(ytr == 1, length(ytr) / (2 * sum(ytr == 1)), length(ytr) / (2 * sum(ytr == 0)))
cvfit <- cv.glmnet(xtr, ytr, family = "binomial", alpha = 1, foldid = foldid,
                   weights = weights, type.measure = "auc", parallel = FALSE)
saveRDS(cvfit, file.path(out, "lasso_glmnet.rds"))
add_results("LASSO", as.numeric(predict(cvfit, xva, s = "lambda.1se", type = "response")),
            as.numeric(predict(cvfit, xte, s = "lambda.1se", type = "response")))

rf_dat <- data.frame(label = factor(ytr, levels = c(0, 1)), xtr, check.names = FALSE)
rf <- ranger(label ~ ., data = rf_dat, probability = TRUE, num.trees = 500,
             mtry = floor(sqrt(ncol(xtr))), min.node.size = 5,
             case.weights = weights, num.threads = 8, seed = 20261027,
             importance = "permutation")
saveRDS(rf, file.path(out, "ranger_random_forest.rds"))
add_results("RandomForest", predict(rf, data.frame(xva, check.names = FALSE), num.threads = 8)$predictions[, "1"],
            predict(rf, data.frame(xte, check.names = FALSE), num.threads = 8)$predictions[, "1"])

dtr <- xgb.DMatrix(xtr, label = ytr, weight = weights)
dva <- xgb.DMatrix(xva, label = yva)
params <- list(objective = "binary:logistic", eval_metric = "auc", eta = 0.03,
               max_depth = 4, min_child_weight = 3, subsample = 0.8,
               colsample_bytree = 0.6, nthread = 8, seed = 20261027)
xgb <- xgb.train(params = params, data = dtr, nrounds = 500,
                 watchlist = list(train = dtr, validation = dva),
                 early_stopping_rounds = 30, verbose = 0)
xgb.save(xgb, file.path(out, "xgboost.ubj"))
add_results("XGBoost", predict(xgb, xva), predict(xgb, xte))

metric_dt <- rbindlist(metric_tables)
pred_dt <- rbindlist(pred_tables)
fwrite(metric_dt, file.path(out, "baseline_metrics.tsv"), sep = "\t")
fwrite(pred_dt, file.path(out, "baseline_predictions.tsv.gz"), sep = "\t")
summary <- list(n_train = length(tr), n_validation = length(va), n_test = length(te),
                n_features = length(sel), feature_selection = "training-only limma top 1000",
                patient_grouped_cv = TRUE, metrics = metric_dt)
write_json(summary, file.path(out, "baseline_summary.json"), pretty = TRUE, auto_unbox = TRUE)
print(metric_dt)

