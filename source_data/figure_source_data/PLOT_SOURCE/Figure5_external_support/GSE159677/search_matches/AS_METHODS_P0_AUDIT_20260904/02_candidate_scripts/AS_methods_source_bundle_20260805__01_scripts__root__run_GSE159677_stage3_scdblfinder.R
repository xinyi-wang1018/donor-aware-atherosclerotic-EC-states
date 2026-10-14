#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 12 * 1024^3)

required_pkgs <- c(
  "Seurat",
  "SingleCellExperiment",
  "scDblFinder",
  "BiocParallel",
  "ggplot2"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them first; this script will not modify the R environment automatically."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
  library(ggplot2)
})

in_rds <- paste0(
  "/root/AS_external_results/GSE159677/stage2_qc_filtered/",
  "GSE159677_qc_filtered_before_doublet.rds"
)
out_dir <- "/root/AS_external_results/GSE159677/stage3_doublet"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "per_sample"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "00_stage3_doublet_log.txt")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("===== GSE159677 STAGE 3: SAMPLE-WISE scDblFinder =====\n")
cat("Start:", format(Sys.time()), "\n")
cat("Input:", in_rds, "\n")
cat("Output:", out_dir, "\n")
cat("scDblFinder version:", as.character(packageVersion("scDblFinder")), "\n\n")

if (!file.exists(in_rds)) {
  stop("Input RDS not found: ", in_rds)
}

obj <- readRDS(in_rds)

required_meta <- c("sample_id", "patient_id", "region")
missing_meta <- setdiff(required_meta, colnames(obj[[]]))
if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
}

get_counts <- function(seurat_obj) {
  assay_name <- "RNA"

  if (inherits(seurat_obj[[assay_name]], "Assay5")) {
    layers <- SeuratObject::Layers(seurat_obj[[assay_name]])
    if (!("counts" %in% layers)) {
      stop("RNA assay has no counts layer.")
    }
    return(
      SeuratObject::LayerData(
        object = seurat_obj,
        assay = assay_name,
        layer = "counts"
      )
    )
  }

  Seurat::GetAssayData(
    object = seurat_obj,
    assay = assay_name,
    slot = "counts"
  )
}

sample_ids <- sort(unique(as.character(obj$sample_id)))
cat("Samples:", paste(sample_ids, collapse = ", "), "\n")
cat("Total QC-passed cells:", ncol(obj), "\n\n")

sample_objects <- Seurat::SplitObject(obj, split.by = "sample_id")

annotated_list <- list()
singlet_list <- list()
summary_rows <- list()
cell_rows <- list()

for (i in seq_along(sample_ids)) {
  sid <- sample_ids[i]
  cat("\n----------------------------------------\n")
  cat("Sample", i, "of", length(sample_ids), ":", sid, "\n")

  sobj <- sample_objects[[sid]]
  if (is.null(sobj)) {
    stop("Sample object not found after SplitObject: ", sid)
  }

  counts <- get_counts(sobj)
  cat("Dimensions:", nrow(counts), "genes x", ncol(counts), "cells\n")

  if (ncol(counts) < 100) {
    stop("Too few cells for reliable doublet detection in ", sid)
  }

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts)
  )

  set.seed(20260729 + i)

  sce <- scDblFinder::scDblFinder(
    sce,
    clusters = TRUE,
    samples = NULL,
    dbr = NULL,
    dbr.per1k = 0.008,
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
    verbose = TRUE
  )

  score <- as.numeric(SummarizedExperiment::colData(sce)$scDblFinder.score)
  cls <- as.character(
    SummarizedExperiment::colData(sce)$scDblFinder.class
  )

  if (length(score) != ncol(sobj) || length(cls) != ncol(sobj)) {
    stop("scDblFinder output length mismatch in ", sid)
  }

  names(score) <- colnames(sce)
  names(cls) <- colnames(sce)

  sobj$scDblFinder_score <- score[colnames(sobj)]
  sobj$scDblFinder_class <- cls[colnames(sobj)]

  n_doublet <- sum(tolower(sobj$scDblFinder_class) == "doublet")
  n_singlet <- sum(tolower(sobj$scDblFinder_class) == "singlet")
  n_total <- ncol(sobj)

  if ((n_doublet + n_singlet) != n_total) {
    stop(
      "Unexpected scDblFinder classes in ", sid, ": ",
      paste(unique(sobj$scDblFinder_class), collapse = ", ")
    )
  }

  doublet_rate <- n_doublet / n_total
  expected_rate_default <- min(1, 0.008 * n_total / 1000)

  summary_rows[[sid]] <- data.frame(
    sample_id = sid,
    patient_id = as.character(unique(sobj$patient_id))[1],
    region = as.character(unique(sobj$region))[1],
    n_before_doublet = n_total,
    n_predicted_singlet = n_singlet,
    n_predicted_doublet = n_doublet,
    predicted_doublet_rate = doublet_rate,
    default_expected_rate_from_cell_count = expected_rate_default,
    median_doublet_score = median(sobj$scDblFinder_score, na.rm = TRUE),
    q95_doublet_score = as.numeric(
      quantile(sobj$scDblFinder_score, 0.95, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  cell_rows[[sid]] <- data.frame(
    cell = colnames(sobj),
    sample_id = sid,
    patient_id = as.character(sobj$patient_id),
    region = as.character(sobj$region),
    nFeature_RNA = sobj$nFeature_RNA,
    nCount_RNA = sobj$nCount_RNA,
    percent.mt = sobj$percent.mt,
    scDblFinder_score = sobj$scDblFinder_score,
    scDblFinder_class = sobj$scDblFinder_class,
    stringsAsFactors = FALSE
  )

  annotated_list[[sid]] <- sobj
  singlet_list[[sid]] <- subset(
    sobj,
    cells = colnames(sobj)[
      tolower(sobj$scDblFinder_class) == "singlet"
    ]
  )

  saveRDS(
    sobj,
    file.path(out_dir, "per_sample", paste0(sid, "_doublet_annotated.rds")),
    compress = FALSE
  )
  saveRDS(
    singlet_list[[sid]],
    file.path(out_dir, "per_sample", paste0(sid, "_singlets.rds")),
    compress = FALSE
  )

  cat(
    "Predicted doublets:", n_doublet,
    sprintf("(%.2f%%)", 100 * doublet_rate), "\n"
  )
  cat("Predicted singlets:", n_singlet, "\n")
}

summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL
summary_df <- summary_df[order(summary_df$sample_id), ]

cell_df <- do.call(rbind, cell_rows)
rownames(cell_df) <- NULL

write.csv(
  summary_df,
  file.path(out_dir, "01_doublet_summary_by_sample.csv"),
  row.names = FALSE
)

write.csv(
  cell_df,
  gzfile(file.path(out_dir, "02_cell_level_doublet_calls.csv.gz")),
  row.names = FALSE
)

cat("\nDoublet summary:\n")
print(summary_df, row.names = FALSE)

cat("\nMerging doublet-annotated sample objects...\n")
annotated_merged <- merge(
  x = annotated_list[[1]],
  y = annotated_list[-1],
  project = "GSE159677_doublet_annotated",
  merge.data = FALSE
)

cat("Merging singlet sample objects...\n")
singlets_merged <- merge(
  x = singlet_list[[1]],
  y = singlet_list[-1],
  project = "GSE159677_singlets",
  merge.data = FALSE
)

saveRDS(
  annotated_merged,
  file.path(out_dir, "GSE159677_qc_doublet_annotated_all_cells.rds"),
  compress = FALSE
)

saveRDS(
  singlets_merged,
  file.path(out_dir, "GSE159677_qc_doublet_clean_singlets.rds"),
  compress = FALSE
)

saveRDS(
  singlet_list,
  file.path(out_dir, "GSE159677_singlets_split_by_sample.rds"),
  compress = FALSE
)

# Plots
p_rate <- ggplot(
  summary_df,
  aes(x = sample_id, y = predicted_doublet_rate)
) +
  geom_col() +
  geom_point(
    aes(y = default_expected_rate_from_cell_count),
    size = 2.5
  ) +
  geom_text(
    aes(
      label = sprintf("%.1f%%", 100 * predicted_doublet_rate)
    ),
    vjust = -0.35,
    size = 3.5
  ) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    x = NULL,
    y = "Predicted doublet rate",
    title = "GSE159677 scDblFinder doublet rate by sample",
    subtitle = "Bars: predicted rate; points: default expected rate from recovered cell count"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(out_dir, "03_doublet_rate_by_sample.png"),
  p_rate,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

p_score <- ggplot(
  cell_df,
  aes(
    x = sample_id,
    y = scDblFinder_score,
    fill = scDblFinder_class
  )
) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.10, outlier.shape = NA) +
  labs(
    x = NULL,
    y = "scDblFinder score",
    title = "GSE159677 doublet score distributions"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(out_dir, "04_doublet_score_by_sample_and_class.png"),
  p_score,
  width = 9,
  height = 5.5,
  dpi = 300
)

p_count <- ggplot(
  cell_df,
  aes(
    x = nCount_RNA,
    y = scDblFinder_score,
    shape = scDblFinder_class
  )
) +
  geom_point(alpha = 0.25, size = 0.6) +
  scale_x_log10() +
  facet_wrap(~ sample_id, scales = "free_x") +
  labs(
    x = "UMI counts (log10 scale)",
    y = "scDblFinder score",
    title = "Doublet score versus RNA content"
  ) +
  theme_classic(base_size = 10)

ggsave(
  file.path(out_dir, "05_doublet_score_vs_nCount.png"),
  p_count,
  width = 11,
  height = 7,
  dpi = 300
)

cat("\nTotal cells before doublet removal:", ncol(annotated_merged), "\n")
cat("Total singlets after doublet removal:", ncol(singlets_merged), "\n")
cat(
  "Overall singlet retain rate:",
  round(ncol(singlets_merged) / ncol(annotated_merged), 4),
  "\n"
)

cat("\nSaved outputs:\n")
for (f in list.files(out_dir, full.names = TRUE, recursive = TRUE)) {
  cat(f, "\n")
}

cat("\n===== STAGE 3 COMPLETE =====\n")
cat("End:", format(Sys.time()), "\n")
