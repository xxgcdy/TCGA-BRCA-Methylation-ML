suppressPackageStartupMessages({
  library(rhdf5)
  library(data.table)
  library(DMRcate)
  library(jsonlite)
  library(matrixStats)
})

root <- Sys.getenv("METHYLATION_ROOT", unset = normalizePath(".", mustWork = TRUE))
out_dir <- file.path(root, "results/dmr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
samples <- fread(file.path(root, "processed/sample_order_qc_autosomal.tsv"))
probes <- fread(file.path(root, "processed/probe_order_qc_autosomal.tsv.gz"))
beta <- h5read(file.path(root, "processed/beta_qc_autosomal.h5"), "beta")
stopifnot(nrow(beta) == nrow(probes), ncol(beta) == nrow(samples))

# Median imputation per probe is used only for residual missing beta values (<1%).
probe_med <- rowMedians(beta, na.rm = TRUE)
idx <- which(is.na(beta), arr.ind = TRUE)
if (nrow(idx) > 0) beta[idx] <- probe_med[idx[, 1]]
eps <- 1e-6
m <- log2(pmin(pmax(beta, eps), 1 - eps) / (1 - pmin(pmax(beta, eps), 1 - eps)))
rownames(m) <- probes$probe_id

group <- factor(samples$sample_type,
                levels = c("Solid Tissue Normal", "Primary Tumor"),
                labels = c("Normal", "Tumor"))
design <- model.matrix(~ group)
colnames(design)[2] <- "TumorVsNormal"

ann <- cpg.annotate(datatype = "array", object = m, what = "M",
                    arraytype = "450K", analysis.type = "differential",
                    design = design, contrasts = FALSE, coef = "TumorVsNormal",
                    fdr = 0.05)
dmr <- dmrcate(ann, lambda = 1000, C = 2, pcutoff = "fdr", min.cpgs = 2)
saveRDS(dmr, file.path(out_dir, "dmrcate_full_cohort.rds"), compress = FALSE)

# Local equivalent of extractRanges(), intentionally omitting ExperimentHub network calls.
coords <- DMRcate:::extractCoords(dmr@coord)
coords <- cbind(coords, dmr@no.cpgs, dmr@min_smoothed_fdr,
                dmr@Stouffer, dmr@HMFDR, dmr@Fisher,
                dmr@maxdiff, dmr@meandiff)
setDT(coords)
# DMRcate currently contributes seven statistics after coordinates; normalize robustly.
stat_names <- c("no_cpgs", "min_smoothed_fdr", "Stouffer", "HMFDR",
                "Fisher", "maxdiff", "meandiff")
setnames(coords, (ncol(coords)-6):ncol(coords), stat_names)
coords[, region_id := sprintf("DMR%06d", .I)]

# Annotate regions with genes represented by HM450 probes inside each interval.
regions <- coords[, .(region_id, chr = as.character(chrom),
                      start = as.integer(chromStart), end = as.integer(chromEnd))]
setkey(regions, chr, start, end)
ann_probes <- probes[!is.na(chr) & !is.na(pos),
                     .(chr = as.character(chr), start = as.integer(pos), end = as.integer(pos),
                       gene = if ("UCSC_RefGene_Name" %in% names(probes)) UCSC_RefGene_Name else "")]
setkey(ann_probes, chr, start, end)
hits <- foverlaps(ann_probes, regions, type = "within", nomatch = 0L)
gene_map <- hits[, .(overlapping_genes = paste(sort(unique(unlist(strsplit(gene[gene != "" & !is.na(gene)], ";", fixed = TRUE)))), collapse = ";")),
                 by = region_id]
tbl <- merge(coords, gene_map, by = "region_id", all.x = TRUE, sort = FALSE)
fwrite(tbl, file.path(out_dir, "dmr_full_cohort.tsv.gz"), sep = "\t")

summary <- list(
  method = "DMRcate full cohort Tumor vs Normal",
  n_samples = nrow(samples),
  n_tumor = sum(group == "Tumor"),
  n_normal = sum(group == "Normal"),
  n_input_probes = nrow(m),
  n_dmrs = nrow(tbl),
  lambda = 1000,
  C = 2,
  min_cpgs = 2
)
write_json(summary, file.path(out_dir, "dmr_full_summary.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), "\n")

