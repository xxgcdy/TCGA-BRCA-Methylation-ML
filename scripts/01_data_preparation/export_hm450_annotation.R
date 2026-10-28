#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))

project <- Sys.getenv("METHYLATION_ROOT", unset = normalizePath(".", mustWork = TRUE))
outdir <- file.path(project, "annotation")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

anno <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
anno$probe_id <- rownames(anno)
anno <- anno[, c("probe_id", setdiff(colnames(anno), "probe_id"))]

outfile <- gzfile(file.path(outdir, "hm450_annotation.tsv.gz"), "wt")
write.table(anno, outfile, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
close(outfile)

summary <- data.frame(
  metric = c(
    "annotation_rows",
    "autosomal",
    "chrX",
    "chrY",
    "missing_chr",
    "probe_snp_maf_ge_0.05",
    "cpg_snp_maf_ge_0.05",
    "sbe_snp_maf_ge_0.05"
  ),
  value = c(
    nrow(anno),
    sum(anno$chr %in% paste0("chr", 1:22), na.rm = TRUE),
    sum(anno$chr == "chrX", na.rm = TRUE),
    sum(anno$chr == "chrY", na.rm = TRUE),
    sum(is.na(anno$chr) | anno$chr == ""),
    sum(!is.na(anno$Probe_maf) & anno$Probe_maf >= 0.05),
    sum(!is.na(anno$CpG_maf) & anno$CpG_maf >= 0.05),
    sum(!is.na(anno$SBE_maf) & anno$SBE_maf >= 0.05)
  )
)
write.table(summary, file.path(outdir, "hm450_annotation_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary, row.names = FALSE)

