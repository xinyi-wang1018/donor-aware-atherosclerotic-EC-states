#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 8 * 1024^3)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

in_rds  <- "/root/AS_external_results/GSE159677/stage1_raw_qc/GSE159677_raw_unfiltered_with_metadata.rds"
out_dir <- "/root/AS_external_results/GSE159677/stage2_qc_filtered"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "00_stage2_qc_log.txt")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("===== GSE159677 STAGE 2: SAMPLE-WISE QC FILTERING =====\n")
cat("Start:", format(Sys.time()), "\n")
cat("Input:", in_rds, "\n")
cat("Output:", out_dir, "\n\n")

if (!file.exists(in_rds)) {
  stop("Input RDS not found: ", in_rds)
}

obj <- readRDS(in_rds)

required_meta <- c(
  "sample_id", "patient_id", "region",
  "nFeature_RNA", "nCount_RNA", "percent.mt"
)
missing_meta <- setdiff(required_meta, colnames(obj[[]]))
if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
}

# Conservative first-pass QC:
# - keep cells with >=200 detected genes and >=500 UMIs
# - remove cells with >20% mitochondrial reads
# - remove the top 1% nFeature/nCount outliers within each biological sample
# The high-count upper tail is not treated as final doublet removal; doublets are handled later.
base_min_features <- 200
base_min_counts   <- 500
max_percent_mt    <- 20
upper_quantile    <- 0.99

meta <- obj[[]]
meta$cell <- rownames(meta)

sample_ids <- unique(as.character(meta$sample_id))
sample_ids <- sample_ids[order(sample_ids)]

threshold_rows <- list()
keep_cells <- character(0)

for (sid in sample_ids) {
  df <- meta[meta$sample_id == sid, , drop = FALSE]

  feature_upper <- as.numeric(
    quantile(df$nFeature_RNA, probs = upper_quantile, na.rm = TRUE)
  )
  count_upper <- as.numeric(
    quantile(df$nCount_RNA, probs = upper_quantile, na.rm = TRUE)
  )

  keep <- (
    df$nFeature_RNA >= base_min_features &
    df$nCount_RNA >= base_min_counts &
    df$percent.mt <= max_percent_mt &
    df$nFeature_RNA <= feature_upper &
    df$nCount_RNA <= count_upper
  )

  keep_cells <- c(keep_cells, df$cell[keep])

  threshold_rows[[sid]] <- data.frame(
    sample_id = sid,
    patient_id = as.character(df$patient_id[1]),
    region = as.character(df$region[1]),
    n_before = nrow(df),
    min_nFeature_RNA = base_min_features,
    max_nFeature_RNA = feature_upper,
    min_nCount_RNA = base_min_counts,
    max_nCount_RNA = count_upper,
    max_percent_mt = max_percent_mt,
    n_after_qc = sum(keep),
    retain_rate = sum(keep) / nrow(df),
    stringsAsFactors = FALSE
  )
}

thresholds <- do.call(rbind, threshold_rows)
rownames(thresholds) <- NULL

write.csv(
  thresholds,
  file.path(out_dir, "01_sample_specific_qc_thresholds_and_retention.csv"),
  row.names = FALSE
)

cat("Sample-specific QC summary:\n")
print(thresholds, row.names = FALSE)

obj_qc <- subset(obj, cells = keep_cells)

# Save per-cell filter status for full auditability.
meta$pass_qc <- meta$cell %in% keep_cells
write.csv(
  meta[, c(
    "cell", "sample_id", "patient_id", "region",
    "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo",
    "pass_qc"
  )],
  gzfile(file.path(out_dir, "02_cell_level_qc_status.csv.gz")),
  row.names = FALSE
)

saveRDS(
  obj_qc,
  file.path(out_dir, "GSE159677_qc_filtered_before_doublet.rds"),
  compress = FALSE
)

# Create split objects so later doublet detection is run per biological sample.
obj_list <- SplitObject(obj_qc, split.by = "sample_id")
saveRDS(
  obj_list,
  file.path(out_dir, "GSE159677_qc_filtered_split_by_sample.rds"),
  compress = FALSE
)

# Before/after counts.
after_meta <- obj_qc[[]]
after_counts <- as.data.frame(table(after_meta$sample_id))
colnames(after_counts) <- c("sample_id", "n_after_qc")
after_counts$sample_id <- as.character(after_counts$sample_id)

before_counts <- as.data.frame(table(meta$sample_id))
colnames(before_counts) <- c("sample_id", "n_before")
before_counts$sample_id <- as.character(before_counts$sample_id)

count_summary <- merge(before_counts, after_counts, by = "sample_id", all = TRUE)
count_summary$n_after_qc[is.na(count_summary$n_after_qc)] <- 0
count_summary$retain_rate <- count_summary$n_after_qc / count_summary$n_before
count_summary <- count_summary[order(count_summary$sample_id), ]

write.csv(
  count_summary,
  file.path(out_dir, "03_cell_counts_before_after_qc.csv"),
  row.names = FALSE
)

# QC visualizations after filtering.
plot_violin <- function(data, variable, y_label, filename, log_scale = FALSE) {
  p <- ggplot(data, aes(x = sample_id, y = .data[[variable]])) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.shape = NA) +
    labs(
      x = NULL,
      y = y_label,
      title = paste("GSE159677", y_label, "after sample-wise QC")
    ) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (log_scale) {
    p <- p + scale_y_log10()
  }

  ggsave(
    file.path(out_dir, filename),
    p,
    width = 8,
    height = 5.5,
    dpi = 300
  )
}

plot_violin(after_meta, "nFeature_RNA", "Detected genes", "04_nFeature_RNA_after_qc.png")
plot_violin(after_meta, "nCount_RNA", "UMI counts", "05_nCount_RNA_after_qc_log10.png", log_scale = TRUE)
plot_violin(after_meta, "percent.mt", "Mitochondrial percentage", "06_percent_mt_after_qc.png")

# Before/after retention bar plot.
retention_plot_df <- count_summary
p_ret <- ggplot(
  retention_plot_df,
  aes(x = sample_id, y = retain_rate)
) +
  geom_col() +
  geom_text(
    aes(label = sprintf("%.1f%%", 100 * retain_rate)),
    vjust = -0.3,
    size = 3.5
  ) +
  scale_y_continuous(
    limits = c(0, 1.08),
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  labs(
    x = NULL,
    y = "QC retention rate",
    title = "GSE159677 cell retention after sample-wise QC"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(out_dir, "07_qc_retention_rate_by_sample.png"),
  p_ret,
  width = 8,
  height = 5.5,
  dpi = 300
)

cat("\nTotal cells before QC:", ncol(obj), "\n")
cat("Total cells after QC:", ncol(obj_qc), "\n")
cat("Overall retain rate:", round(ncol(obj_qc) / ncol(obj), 4), "\n")

cat("\nSaved outputs:\n")
for (f in list.files(out_dir, full.names = TRUE)) {
  cat(f, "\n")
}

cat("\nIMPORTANT: this object is before doublet removal.\n")
cat("===== STAGE 2 COMPLETE =====\n")
cat("End:", format(Sys.time()), "\n")
