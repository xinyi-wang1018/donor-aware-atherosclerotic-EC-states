suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(ggplot2)
  library(Matrix)
})

set.seed(20260623)

# =========================================================
# 1. Paths
# =========================================================
input_rds <- "/root/AS_geneid_fixed_objects/all_53_qc_doublet_clean_geneid_fixed_objs.rds"

out_dir <- "/root/AS_53samples_geneidfixed_RNA_harmony_cleanADT"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

final_rds <- file.path(
  out_dir,
  "AS_53samples_geneidfixed_cleanADT_RNA_harmony_umap_clustered.rds"
)

# =========================================================
# 2. Helper: remove ADT-like features accidentally stored in RNA
# =========================================================
clean_rna_features <- function(obj, object_index) {

  DefaultAssay(obj) <- "RNA"

  layer_names <- Layers(obj[["RNA"]])
  count_layers <- grep("^counts", layer_names, value = TRUE)

  if (length(count_layers) == 0) {
    stop("Object ", object_index, " has no RNA counts layer.")
  }

  if (length(count_layers) > 1 || !"counts" %in% layer_names) {
    obj <- JoinLayers(obj, assay = "RNA")
  }

  counts_before <- LayerData(
    object = obj,
    assay = "RNA",
    layer = "counts"
  )

  genes_before <- rownames(counts_before)

  # 只移除明确以 -ADT 或 _ADT 结尾的抗体标签 feature
  adt_idx <- grepl("(-ADT|_ADT)$", genes_before, ignore.case = TRUE)

  n_adt_features <- sum(adt_idx)
  total_counts_before <- sum(counts_before)

  adt_counts_removed <- if (n_adt_features > 0) {
    sum(counts_before[adt_idx, , drop = FALSE])
  } else {
    0
  }

  if (n_adt_features > 0) {

    if ("nCount_RNA" %in% colnames(obj@meta.data)) {
      obj$nCount_RNA_before_ADT_removal <- obj$nCount_RNA
    }

    if ("nFeature_RNA" %in% colnames(obj@meta.data)) {
      obj$nFeature_RNA_before_ADT_removal <- obj$nFeature_RNA
    }

    keep_genes <- genes_before[!adt_idx]

    obj <- subset(
      x = obj,
      features = keep_genes
    )

    counts_after <- LayerData(
      object = obj,
      assay = "RNA",
      layer = "counts"
    )

    obj$nCount_RNA <- as.numeric(Matrix::colSums(counts_after))
    obj$nFeature_RNA <- as.numeric(Matrix::colSums(counts_after > 0))
    obj$ADT_feature_cleanup <- "removed_ADT_features_from_RNA"

    total_counts_after <- sum(counts_after)

  } else {

    obj$ADT_feature_cleanup <- "no_ADT_features_detected"
    total_counts_after <- total_counts_before
  }

  report <- data.frame(
    object_index = object_index,
    n_cells = ncol(obj),
    n_features_before = length(genes_before),
    n_ADT_features_removed = n_adt_features,
    ADT_counts_removed = adt_counts_removed,
    total_counts_before = total_counts_before,
    total_counts_after = total_counts_after,
    count_retention_fraction = total_counts_after / max(total_counts_before, 1),
    stringsAsFactors = FALSE
  )

  list(object = obj, report = report)
}

# =========================================================
# 3. Load repaired 53-sample list
# =========================================================
obj_list <- readRDS(input_rds)

if (!is.list(obj_list) || length(obj_list) != 53) {
  stop("Expected exactly 53 Seurat objects.")
}

cat("Loaded repaired object list: ", length(obj_list), " objects.\n", sep = "")

validation_list <- vector("list", length(obj_list))
adt_cleanup_list <- vector("list", length(obj_list))

# =========================================================
# 4. Validate, remove ADT features, and prepare each sample
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

  # 删除额外 assay；本轮仅做 RNA 图谱
  extra_assays <- setdiff(Assays(obj), "RNA")

  if (length(extra_assays) > 0) {
    for (assay_name in extra_assays) {
      obj[[assay_name]] <- NULL
    }
  }

  cleaned <- clean_rna_features(obj, i)
  obj <- cleaned$object
  adt_cleanup_list[[i]] <- cleaned$report

  # 保留原 barcode，并让合并后的细胞名唯一
  obj$original_cell_barcode <- Cells(obj)
  obj$reintegrated_sample_index <- i

  obj <- RenameCells(
    object = obj,
    add.cell.id = sprintf("S%02d", i)
  )

  DefaultAssay(obj) <- "RNA"

  genes <- rownames(obj[["RNA"]])

  n_ensg <- sum(grepl("^ENSG[0-9]+(\\.[0-9]+)?$", genes))
  n_duplicates <- sum(duplicated(genes))
  n_adt_after <- sum(grepl("(-ADT|_ADT)$", genes, ignore.case = TRUE))

  if (n_ensg > 0) {
    stop(
      "Object ", i, " still has ", n_ensg,
      " ENSG feature names."
    )
  }

  if (n_duplicates > 0) {
    stop(
      "Object ", i, " has ", n_duplicates,
      " duplicated gene names."
    )
  }

  if (n_adt_after > 0) {
    stop(
      "Object ", i, " still has ", n_adt_after,
      " ADT-like features in RNA."
    )
  }

  sample_id <- paste(
    unique(as.character(obj$sample_id)),
    collapse = ";"
  )

  dataset_id <- paste(
    unique(as.character(obj$dataset_id)),
    collapse = ";"
  )

  validation_list[[i]] <- data.frame(
    object_index = i,
    sample_id = sample_id,
    dataset_id = dataset_id,
    n_cells = ncol(obj),
    n_features = nrow(obj),
    n_ENSG_features = n_ensg,
    n_ADT_features_after_cleanup = n_adt_after,
    duplicated_gene_names = n_duplicates,
    stringsAsFactors = FALSE
  )

  obj_list[[i]] <- obj

  cat(
    sprintf(
      "[%02d/53] %s | %s | %d cells | %d genes\n",
      i, sample_id, dataset_id, ncol(obj), nrow(obj)
    )
  )
}

validation_df <- do.call(rbind, validation_list)
adt_cleanup_df <- do.call(rbind, adt_cleanup_list)

write.csv(
  validation_df,
  file.path(out_dir, "input_samples_validation.csv"),
  row.names = FALSE
)

write.csv(
  adt_cleanup_df,
  file.path(out_dir, "ADT_feature_cleanup_by_sample.csv"),
  row.names = FALSE
)

cat("\nAll samples passed ENSG / ADT / duplicate-gene checks.\n")

# =========================================================
# 5. Merge all samples
# =========================================================
cat("\nMerging 53 samples...\n")

combined <- merge(
  x = obj_list[[1]],
  y = obj_list[2:length(obj_list)],
  merge.data = FALSE
)

rm(obj_list)
gc(verbose = FALSE)

DefaultAssay(combined) <- "RNA"

rna_layers <- Layers(combined[["RNA"]])

if (length(rna_layers) > 1) {
  cat("Joining RNA layers...\n")
  combined <- JoinLayers(combined, assay = "RNA")
}

rna_layers <- Layers(combined[["RNA"]])

if (!"counts" %in% rna_layers) {
  stop(
    "Merged object does not contain standard counts layer. Layers: ",
    paste(rna_layers, collapse = ", ")
  )
}

merged_genes <- rownames(combined[["RNA"]])

if (sum(grepl("(-ADT|_ADT)$", merged_genes, ignore.case = TRUE)) > 0) {
  stop("ADT-like features remain after merging.")
}

if (sum(grepl("^ENSG[0-9]+(\\.[0-9]+)?$", merged_genes)) > 0) {
  stop("ENSG feature names remain after merging.")
}

cat("Merged cells: ", ncol(combined), "\n", sep = "")
cat("Merged genes: ", nrow(combined), "\n", sep = "")

# =========================================================
# 6. RNA preprocessing and PCA
# =========================================================
cat("\nCalculating mitochondrial percentage...\n")

combined[["percent.mt"]] <- PercentageFeatureSet(
  combined,
  pattern = "^MT-"
)

cat("NormalizeData...\n")

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

saveRDS(
  combined,
  file.path(out_dir, "AS_53samples_geneidfixed_cleanADT_RNA_pca.rds")
)

# =========================================================
# 7. Harmony
# 注意：此服务器 Harmony 版本使用 reduction.use，不使用 reduction
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
# 9. Save object and diagnostics
# =========================================================
cluster_sizes <- as.data.frame(table(combined$seurat_clusters))
colnames(cluster_sizes) <- c("cluster", "n_cells")

write.csv(
  cluster_sizes,
  file.path(out_dir, "cluster_cell_numbers.csv"),
  row.names = FALSE
)

composition_count <- table(
  combined$seurat_clusters,
  combined$dataset_id
)

composition_prop <- prop.table(composition_count, margin = 1)

write.csv(
  as.data.frame.matrix(composition_count),
  file.path(out_dir, "cluster_by_dataset_cell_counts.csv")
)

write.csv(
  as.data.frame.matrix(composition_prop),
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
  ggtitle("AS atlas after gene-ID and ADT cleanup: clusters") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(out_dir, "UMAP_clusters_cleanADT.png"),
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
  ggtitle("AS atlas after gene-ID and ADT cleanup: datasets") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(out_dir, "UMAP_dataset_cleanADT.png"),
  p_dataset,
  width = 14,
  height = 10,
  dpi = 300,
  bg = "white"
)

saveRDS(combined, final_rds)

writeLines(
  c(
    "========== REINTEGRATION SUMMARY ==========",
    paste("Input repaired object list:", input_rds),
    paste("Final object:", final_rds),
    paste("Total cells:", ncol(combined)),
    paste("Total genes:", nrow(combined)),
    paste("Number of clusters:", length(unique(combined$seurat_clusters))),
    "RNA preprocessing: LogNormalize + 3000 HVGs + percent.mt regression",
    "Harmony batch variable: dataset_id",
    "Harmony interface: reduction.use = pca",
    "RNA cleanup: features ending in -ADT or _ADT removed before PCA",
    "Clustering resolution: 0.50",
    "==========================================="
  ),
  file.path(out_dir, "run_summary.txt")
)

cat("\n====================================================\n")
cat("Step 6B finished successfully.\n")
cat("Final object:", final_rds, "\n")
cat("Output directory:", out_dir, "\n")
cat("====================================================\n")
