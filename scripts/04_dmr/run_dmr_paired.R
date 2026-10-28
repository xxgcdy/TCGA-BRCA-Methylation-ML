suppressPackageStartupMessages({
  library(rhdf5)
  library(data.table)
  library(DMRcate)
  library(jsonlite)
  library(matrixStats)
})

root <- Sys.getenv("METHYLATION_ROOT", unset = normalizePath(".", mustWork = TRUE))
out_dir <- file.path(root, "results/dmr")
samples <- fread(file.path(root, "processed/sample_order_qc_autosomal.tsv"))
probes <- fread(file.path(root, "processed/probe_order_qc_autosomal.tsv.gz"))
pairs <- fread(file.path(root, "results/dmp/paired_samples_used.tsv"))
beta <- h5read(file.path(root, "processed/beta_qc_autosomal.h5"), "beta")
stopifnot(nrow(beta) == nrow(probes), ncol(beta) == nrow(samples))

tumor_cols <- match(pairs$tumor_sample_id, samples$sample_id)
normal_cols <- match(pairs$normal_sample_id, samples$sample_id)
stopifnot(!anyNA(tumor_cols), !anyNA(normal_cols))
probe_med <- rowMedians(beta, na.rm = TRUE)
idx <- which(is.na(beta), arr.ind = TRUE)
if (nrow(idx) > 0) beta[idx] <- probe_med[idx[, 1]]
eps <- 1e-6
to_m <- function(x) log2(pmin(pmax(x, eps), 1 - eps) / (1 - pmin(pmax(x, eps), 1 - eps)))
m_diff <- to_m(beta[, tumor_cols, drop = FALSE]) - to_m(beta[, normal_cols, drop = FALSE])
rownames(m_diff) <- probes$probe_id
design <- matrix(1, nrow = ncol(m_diff), ncol = 1,
                 dimnames = list(NULL, "(Intercept)"))

ann <- cpg.annotate(datatype = "array", object = m_diff, what = "M",
                    arraytype = "450K", analysis.type = "differential",
                    design = design, contrasts = FALSE, coef = "(Intercept)",
                    fdr = 0.05)
dmr <- dmrcate(ann, lambda = 1000, C = 2, pcutoff = "fdr", min.cpgs = 2)
saveRDS(dmr, file.path(out_dir, "dmrcate_paired.rds"), compress = FALSE)

coords <- DMRcate:::extractCoords(dmr@coord)
coords <- cbind(coords, dmr@no.cpgs, dmr@min_smoothed_fdr,
                dmr@Stouffer, dmr@HMFDR, dmr@Fisher,
                dmr@maxdiff, dmr@meandiff)
setDT(coords)
setnames(coords, (ncol(coords)-6):ncol(coords),
         c("no_cpgs", "min_smoothed_fdr", "Stouffer", "HMFDR",
           "Fisher", "maxdiff", "meandiff"))
coords[, region_id := sprintf("PDMR%06d", .I)]
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
fwrite(tbl, file.path(out_dir, "dmr_paired.tsv.gz"), sep = "\t")

summary <- list(method = "DMRcate paired tumor-minus-normal M-value differences",
                n_pairs = nrow(pairs), n_input_probes = nrow(m_diff), n_dmrs = nrow(tbl),
                lambda = 1000, C = 2, min_cpgs = 2)
write_json(summary, file.path(out_dir, "dmr_paired_summary.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), "\n")

