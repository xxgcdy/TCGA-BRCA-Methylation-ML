suppressPackageStartupMessages(library(xgboost))

root <- Sys.getenv("METHYLATION_ROOT", unset = getwd())
model_path <- file.path(root, "results", "models", "xgboost.ubj")
if (!file.exists(model_path)) stop("Missing model: ", model_path)

model <- xgb.load(model_path)
cat("Model attributes\n")
print(xgb.attributes(model))
cat("Parameters retained by R object\n")
print(model$params)
cat("Evaluation log\n")
print(model$evaluation_log)
