suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(ggplot2)
  library(Matrix)
})

set.seed(20260623)

# =========================================================
# 1. Input / output
# =========================================================
input_rds <- paste0(
  "/root/AS_RNAonly_clean_objects/",
  "all_53_qc_doublet_clean_geneid_CITEHTO_RNAonly_objs.rds"
)

out_dir <- "/root/AS_final_RNAonly_harmony"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

final_rds <- file.path(
  out_dir,
  "AS_53samples_final_RNAonly_harmony_umap_clustered.rds"
)

# =========================================================
# 2. Exact checks for remaining technical features
# =========================================================
known_cite_tags <- c(
  "CD3", "CD4.1", "CD8",
  "CD14.1", "CD16", "CD19.1",
  "CD24.1", "CD27.1", "CD33.1",
  "CD38.1", "CD56", "CD62L",
  "CD64", "CD66b", "CD123",
  "CD141", "CD206", "CD11b",
  "CD11c", "HLA-DR"
)

is_technical_feature <- function(genes) {
  grepl("^HTO[-_]?[0-9]+$", genes, ignore.case = TRUE) |
    grepl("(-ADT|_ADT)$", genes, ignore.case = TRUE) |
    genes %in% known_cite_tags
}

# =========================================================
# 3. Load RNA-only object list
# =========================================================
obj_list <- readRDS(input_rds)

if (!is.list(obj_list) || length(obj_list) != 53) {
  stop("Expected exactly 53 Seurat objects.")
}

cat("Loaded RNA-only object list:", length(obj_list), "objects.\n")

validation_list <- vector("list", length(obj_list))

# =========================================================
# 4. Validate and prepare each sample
# =========================================================
for (i in seq_along(obj_list)) {

  obj <- obj_list[[i]]

  if (!inherits(obj, "Seurat")) {
    stop("Object ", i, " is not a Seurat object.")
  }

  if (!"RNA" %in% Assays(obj)) {
    stop("Object ", i, " has no RNA assay.")
  }

  if (!"dataset_id" %in% colnames(obj@meta.data)) {
    stop("Object ", i, " has no dataset_id metadata.")
  }

  DefaultAssay(obj) <- "RNA"

  # 保留原条形码；后续统一加样本前缀，避免 merge 时细胞名重复
  obj$original_cell_barcode <- Cells(obj)
  obj$final_reintegration_sample_index <- i

  obj <- RenameCells(
    object = obj,
    add.cell.id = sprintf("S%02d", i)
  )

  layers <- Layers(obj[["RNA"]])
  count_layers <- grep("^counts", layers, value = TRUE)

  if (length(count_layers) == 0) {
    stop("Object ", i, " has no RNA counts layer.")
  }

  if (length(count_layers) > 1 || !"counts" %in% layers) {
    obj <- JoinLayers(obj, assay = "RNA")
  }

  genes <- rownames(obj[["RNA"]])

  n_ensg <- sum(grepl("^ENSG[0-9]+(\\.[0-9]+)?$", genes))
  n_duplicates <- sum(duplicated(genes))
  n_technical <- sum(is_technical_feature(genes))

  if (n_ensg > 0) {
    stop("Object ", i, " still contains ENSG feature names.")
  }

  if (n_duplicates > 0) {
    stop("Object ", i, " contains duplicated gene names.")
  }

  if (n_technical > 0) {
    stop(
      "Object ", i,
      " still contains ", n_technical,
      " CITE/HTO/ADT technical features."
    )
  }

  sample_id <- if ("sample_id" %in% colnames(obj@meta.data)) {
    paste(unique(as.character(obj$sample_id)), collapse = ";")
  } else {
    paste0("object_", i)
  }

  dataset_id <- paste(
    unique(as.character(obj$dataset_id)),
    collapse = ";"
  )

  validation_list[[i]] <- data.frame(
    object_index = i,
    sample_id = sample_id,
    dataset_id = dataset_id,
    n_cells = ncol(obj),
    n_RNA_features = nrow(obj),
    n_ENSG_features = n_ensg,
    n_technical_features = n_technical,
    duplicated_gene_names = n_duplicates,
    stringsAsFactors = FALSE
  )

  obj_list[[i]] <- obj

  cat(
    sprintf(
      "[%02d/53] %s | %s | %d cells | %d RNA genes\n",
      i, sample_id, dataset_id, ncol(obj), nrow(obj)
    )
  )
}

validation_df <- do.call(rbind, validation_list)

write.csv(
  validation_df,
  file.path(out_dir, "final_input_validation.csv"),
  row.names = FALSE
)

cat("\nAll 53 samples passed RNA-only validation.\n")

# =========================================================
# 5. Merge
# =========================================================
cat("\nMerging 53 RNA-only samples...\n")

combined <- merge(
  x = obj_list[[1]],
  y = obj_list[2:length(obj_list)],
  merge.data = FALSE
)

rm(obj_list)
gc(verbose = FALSE)

DefaultAssay(combined) <- "RNA"

layers <- Layers(combined[["RNA"]])

if (length(layers) > 1) {
  cat("Joining RNA layers after merge...\n")
  combined <- JoinLayers(combined, assay = "RNA")
}

final_genes <- rownames(combined[["RNA"]])

if (sum(is_technical_feature(final_genes)) > 0) {
  stop("Technical CITE/HTO features remain after merging.")
}

if (sum(grepl("^ENSG[0-9]+(\\.[0-9]+)?$", final_genes)) > 0) {
  stop("ENSG features remain after merging.")
}

if (sum(duplicated(final_genes)) > 0) {
  stop("Duplicated RNA features remain after merging.")
}

cat("Merged cells:", ncol(combined), "\n")
cat("Merged RNA genes:", nrow(combined), "\n")

# =========================================================
# 6. RNA-only preprocessing
# =========================================================
combined[["percent.mt"]] <- PercentageFeatureSet(
  combined,
  pattern = "^MT-"
)

cat("\nNormalizeData...\n")

combined <- NormalizeData(
  combined,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

cat("FindVariableFeatures...\n")

combined <- FindVariableFeatures(
  combined,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

cat("ScaleData...\n")

combined <- ScaleData(
  combined,
  features = VariableFeatures(combined),
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

cat("RunPCA...\n")

combined <- RunPCA(
  combined,
  features = VariableFeatures(combined),
  npcs = 50,
  verbose = FALSE
)

# =========================================================
# 7. Harmony
# 当前服务器 harmony 版本必须用 reduction.use
# =========================================================
cat("\nRunning Harmony using dataset_id...\n")

combined <- harmony::RunHarmony(
  object = combined,
  group.by.vars = "dataset_id",
  reduction.use = "pca",
  dims.use = 1:40,
  reduction.save = "harmony",
  plot_convergence = FALSE,
  verbose = TRUE
)

# =========================================================
# 8. UMAP and clustering
# =========================================================
cat("\nRunning UMAP...\n")

combined <- RunUMAP(
  combined,
  reduction = "harmony",
  dims = 1:40,
  reduction.name = "umap",
  reduction.key = "UMAP_",
  n.neighbors = 30,
  min.dist = 0.30,
  seed.use = 20260623,
  verbose = TRUE
)

cat("Finding neighbors...\n")

combined <- FindNeighbors(
  combined,
  reduction = "harmony",
  dims = 1:40,
  verbose = TRUE
)

cat("Finding clusters...\n")

combined <- FindClusters(
  combined,
  resolution = 0.50,
  algorithm = 1,
  random.seed = 20260623,
  verbose = TRUE
)

# =========================================================
# 9. Diagnostics and plots
# =========================================================
cluster_sizes <- as.data.frame(table(combined$seurat_clusters))
colnames(cluster_sizes) <- c("cluster", "n_cells")

write.csv(
  cluster_sizes,
  file.path(out_dir, "cluster_cell_numbers.csv"),
  row.names = FALSE
)

dataset_counts <- table(
  combined$seurat_clusters,
  combined$dataset_id
)

dataset_props <- prop.table(dataset_counts, margin = 1)

write.csv(
  as.data.frame.matrix(dataset_counts),
  file.path(out_dir, "cluster_by_dataset_counts.csv")
)

write.csv(
  as.data.frame.matrix(dataset_props),
  file.path(out_dir, "cluster_by_dataset_proportions.csv")
)

p_cluster <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 4,
  pt.size = 0.03,
  raster = TRUE
) +
  ggtitle("AS single-cell atlas: final RNA-only Harmony clusters") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(out_dir, "UMAP_clusters_final_RNAonly.png"),
  p_cluster,
  width = 14,
  height = 10,
  dpi = 300,
  bg = "white"
)

p_dataset <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "dataset_id",
  pt.size = 0.03,
  raster = TRUE
) +
  ggtitle("AS single-cell atlas: final RNA-only Harmony by dataset") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(out_dir, "UMAP_dataset_final_RNAonly.png"),
  p_dataset,
  width = 14,
  height = 10,
  dpi = 300,
  bg = "white"
)

# =========================================================
# 10. Save final object
# =========================================================
saveRDS(combined, final_rds)

report_lines <- c(
  "========== FINAL RNA-ONLY REINTEGRATION SUMMARY ==========",
  paste("Input RNA-only list:", input_rds),
  paste("Final object:", final_rds),
  paste("Total cells:", ncol(combined)),
  paste("Total RNA genes:", nrow(combined)),
  paste("Number of clusters:", length(unique(combined$seurat_clusters))),
  "Gene-ID repair: completed before integration.",
  "CITE/ADT/HTO cleanup: exact audited blacklist removed before integration.",
  "Normalization: LogNormalize; 3000 HVGs; percent.mt regressed.",
  "Harmony batch variable: dataset_id.",
  "Harmony dimensions: 1-40.",
  "Clustering resolution: 0.50.",
  "=========================================================="
)

writeLines(
  report_lines,
  file.path(out_dir, "run_summary.txt")
)

cat(paste(report_lines, collapse = "\n"), "\n")
cat("\nStep 8B finished successfully.\n")
