suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})
root <- Sys.getenv("METHYLATION_ROOT", unset = normalizePath(".", mustWork = TRUE))
out <- file.path(root, "results/dmr")
full <- fread(file.path(out, "dmr_full_cohort.tsv.gz"))
paired <- fread(file.path(out, "dmr_paired.tsv.gz"))

f <- full[, .(full_region_id = region_id, chr = chrom,
              start = as.integer(chromStart), end = as.integer(chromEnd),
              full_no_cpgs = no_cpgs, full_fdr = min_smoothed_fdr,
              full_maxdiff = maxdiff, full_meandiff = meandiff,
              full_genes = overlapping_genes)]
p <- paired[, .(paired_region_id = region_id, chr = chrom,
                start = as.integer(chromStart), end = as.integer(chromEnd),
                paired_no_cpgs = no_cpgs, paired_fdr = min_smoothed_fdr,
                paired_maxdiff = maxdiff, paired_meandiff = meandiff,
                paired_genes = overlapping_genes)]
setkey(f, chr, start, end)
setkey(p, chr, start, end)
ov <- foverlaps(f, p, type = "any", nomatch = 0L)
ov[, overlap_bp := pmax(0L, pmin(end, i.end) - pmax(start, i.start) + 1L)]
ov[, direction_consistent := sign(full_meandiff) == sign(paired_meandiff)]
ov[, high_confidence := direction_consistent & full_no_cpgs >= 3 & paired_no_cpgs >= 3 &
                         abs(full_meandiff) >= 0.10 & abs(paired_meandiff) >= 0.10 &
                         abs(full_maxdiff) >= 0.20 & abs(paired_maxdiff) >= 0.20]
ov[, `:=`(abs_full_meandiff = abs(full_meandiff),
          abs_paired_meandiff = abs(paired_meandiff))]
setorder(ov, -high_confidence, -overlap_bp, -abs_full_meandiff, -abs_paired_meandiff)
fwrite(ov, file.path(out, "dmr_full_paired_overlaps.tsv.gz"), sep = "\t")
fwrite(ov[high_confidence == TRUE], file.path(out, "dmr_high_confidence.tsv.gz"), sep = "\t")

summary <- list(
  full_dmrs = nrow(full), paired_dmrs = nrow(paired),
  full_regions_with_any_paired_overlap = uniqueN(ov$full_region_id),
  paired_regions_with_any_full_overlap = uniqueN(ov$paired_region_id),
  overlap_pairs = nrow(ov),
  direction_consistent_overlap_pairs = sum(ov$direction_consistent),
  high_confidence_overlap_pairs = sum(ov$high_confidence),
  high_confidence_unique_full_regions = uniqueN(ov[high_confidence == TRUE]$full_region_id),
  high_confidence_unique_paired_regions = uniqueN(ov[high_confidence == TRUE]$paired_region_id),
  high_confidence_definition = "overlap, same direction, both >=3 CpGs, abs mean diff >=0.10, abs max diff >=0.20"
)
write_json(summary, file.path(out, "dmr_overlap_summary.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), "\n")

