#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 8 * 1024^3)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

data_dir <- "/root/AS_external_raw/GSE159677"
out_dir  <- "/root/AS_external_results/GSE159677/stage1_raw_qc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "00_stage1_log.txt")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("===== GSE159677 STAGE 1: RAW IMPORT AND QC AUDIT =====\n")
cat("Start:", format(Sys.time()), "\n")
cat("Input:", data_dir, "\n")
cat("Output:", out_dir, "\n\n")

required <- file.path(
  data_dir,
  c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
)
if (!all(file.exists(required))) {
  stop(
    "Missing required files:\n",
    paste(required[!file.exists(required)], collapse = "\n")
  )
}

# Published suffix mapping for the aggregated GSE159677 matrix:
# -1 Patient 1 PA; -2 Patient 1 AC
# -3 Patient 2 PA; -4 Patient 2 AC
# -5 Patient 3 PA; -6 Patient 3 AC
sample_map <- data.frame(
  suffix = as.character(1:6),
  patient_id = rep(c("P1", "P2", "P3"), each = 2),
  region = rep(c("PA", "AC"), 3),
  sample_id = c("P1_PA", "P1_AC", "P2_PA", "P2_AC", "P3_PA", "P3_AC"),
  stringsAsFactors = FALSE
)

write.csv(
  sample_map,
  file.path(out_dir, "01_suffix_to_sample_mapping.csv"),
  row.names = FALSE
)

cat("Reading 10x matrix...\n")
counts <- Read10X(
  data.dir = data_dir,
  gene.column = 2,
  unique.features = TRUE,
  strip.suffix = FALSE
)

if (is.list(counts)) {
  if ("Gene Expression" %in% names(counts)) {
    counts <- counts[["Gene Expression"]]
  } else if (length(counts) == 1) {
    counts <- counts[[1]]
  } else {
    stop(
      "Read10X returned multiple feature types, but Gene Expression was not found. Types: ",
      paste(names(counts), collapse = ", ")
    )
  }
}

cat("Matrix dimensions:", nrow(counts), "genes x", ncol(counts), "cells\n")

barcodes <- colnames(counts)
suffix <- sub("^.*-([1-6])$", "\\1", barcodes)

invalid <- !(suffix %in% sample_map$suffix)
if (any(invalid)) {
  bad <- unique(barcodes[invalid])
  stop(
    "Some barcodes do not end in -1 to -6. Examples:\n",
    paste(head(bad, 20), collapse = "\n")
  )
}

cell_meta <- data.frame(
  cell_barcode = barcodes,
  suffix = suffix,
  stringsAsFactors = FALSE
)
cell_meta <- merge(
  cell_meta,
  sample_map,
  by = "suffix",
  all.x = TRUE,
  sort = FALSE
)
cell_meta <- cell_meta[match(barcodes, cell_meta$cell_barcode), , drop = FALSE]
rownames(cell_meta) <- cell_meta$cell_barcode

if (anyNA(cell_meta$sample_id)) {
  stop("Sample metadata mapping produced NA values.")
}

sample_counts <- as.data.frame(table(cell_meta$sample_id))
colnames(sample_counts) <- c("sample_id", "n_cells")
sample_counts$sample_id <- as.character(sample_counts$sample_id)
sample_counts <- merge(sample_map, sample_counts, by = "sample_id", all.x = TRUE)
sample_counts <- sample_counts[match(sample_map$sample_id, sample_counts$sample_id), ]

write.csv(
  sample_counts,
  file.path(out_dir, "02_raw_cell_counts_by_sample.csv"),
  row.names = FALSE
)
write.csv(
  cell_meta,
  gzfile(file.path(out_dir, "03_cell_level_metadata.csv.gz")),
  row.names = FALSE
)

cat("\nRaw cell counts by sample:\n")
print(sample_counts, row.names = FALSE)

cat("\nCreating Seurat object without filtering...\n")
obj <- CreateSeuratObject(
  counts = counts,
  project = "GSE159677",
  min.cells = 0,
  min.features = 0
)

obj$dataset_id <- "GSE159677"
obj$suffix <- cell_meta[colnames(obj), "suffix"]
obj$patient_id <- cell_meta[colnames(obj), "patient_id"]
obj$region <- cell_meta[colnames(obj), "region"]
obj$sample_id <- cell_meta[colnames(obj), "sample_id"]
obj$orig.ident <- obj$sample_id
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Additional basic metrics useful for audit.
obj[["percent.ribo"]] <- PercentageFeatureSet(
  obj,
  pattern = "^RP[SL][0-9]"
)

meta <- obj[[]]
meta$cell_barcode <- rownames(meta)

quantile_summary <- function(df, variable) {
  values <- df[[variable]]
  qs <- quantile(
    values,
    probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),
    na.rm = TRUE,
    names = FALSE
  )
  data.frame(
    metric = variable,
    min = qs[1],
    q01 = qs[2],
    q05 = qs[3],
    q25 = qs[4],
    median = qs[5],
    q75 = qs[6],
    q95 = qs[7],
    q99 = qs[8],
    max = qs[9],
    stringsAsFactors = FALSE
  )
}

qc_by_sample <- do.call(
  rbind,
  lapply(
    split(meta, meta$sample_id),
    function(df) {
      rbind(
        cbind(sample_id = df$sample_id[1], quantile_summary(df, "nFeature_RNA")),
        cbind(sample_id = df$sample_id[1], quantile_summary(df, "nCount_RNA")),
        cbind(sample_id = df$sample_id[1], quantile_summary(df, "percent.mt")),
        cbind(sample_id = df$sample_id[1], quantile_summary(df, "percent.ribo"))
      )
    }
  )
)
rownames(qc_by_sample) <- NULL

write.csv(
  qc_by_sample,
  file.path(out_dir, "04_qc_quantiles_by_sample.csv"),
  row.names = FALSE
)

saveRDS(
  obj,
  file.path(out_dir, "GSE159677_raw_unfiltered_with_metadata.rds"),
  compress = FALSE
)

cat("\nSaved raw unfiltered Seurat object.\n")

plot_violin <- function(variable, y_label, filename, log_scale = FALSE) {
  p <- ggplot(meta, aes(x = sample_id, y = .data[[variable]])) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.shape = NA) +
    labs(x = NULL, y = y_label, title = paste("GSE159677", y_label, "before QC")) +
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

plot_violin("nFeature_RNA", "Detected genes", "05_nFeature_RNA_by_sample.png")
plot_violin("nCount_RNA", "UMI counts", "06_nCount_RNA_by_sample_log10.png", log_scale = TRUE)
plot_violin("percent.mt", "Mitochondrial percentage", "07_percent_mt_by_sample.png")
plot_violin("percent.ribo", "Ribosomal percentage", "08_percent_ribo_by_sample.png")

cat("\nOutputs:\n")
for (f in list.files(out_dir, full.names = TRUE)) {
  cat(f, "\n")
}

cat("\n===== STAGE 1 COMPLETE =====\n")
cat("End:", format(Sys.time()), "\n")
