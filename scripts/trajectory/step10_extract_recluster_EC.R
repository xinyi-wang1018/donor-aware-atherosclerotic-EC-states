suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(Matrix)
  library(ggplot2)
})

set.seed(20260623)

# =========================================================
# 1. Path
# =========================================================
input_rds <- paste0(
  "/root/AS_final_RNAonly_harmony/",
  "AS_53samples_final_RNAonly_harmony_umap_clustered.rds"
)

out_dir <- "/root/AS_EC_reclustering"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ec_final_rds <- file.path(
  out_dir,
  "AS_EC_Harmony_umap_clustered.rds"
)

# =========================================================
# 2. Load final whole atlas
# =========================================================
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"

if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
  stop("seurat_clusters not found in final whole-atlas object.")
}

if (!"dataset_id" %in% colnames(obj@meta.data)) {
  stop("dataset_id not found in final whole-atlas object.")
}

if (!"sample_id" %in% colnames(obj@meta.data)) {
  stop("sample_id not found in final whole-atlas object.")
}

# =========================================================
# 3. Final first-level annotation table
#    This saves a formal annotation record without overwriting input RDS.
# =========================================================
annotation_map <- data.frame(
  cluster = as.character(0:29),
  major_celltype = c(
    "IL7R_CD4_T",
    "Cytotoxic_T_GZMK",
    "VSMC_fibromyocyte",
    "C1QC_macrophage",
    "NK_cell",
    "SPP1_macrophage",
    "B_cell",
    "Inflammatory_monocyte",
    "Activated_CD4_T_dataset_enriched",
    "Cytotoxic_T_NK_like",
    "Activated_cytotoxic_lymphocyte_dataset_enriched",
    "Endothelial_EC1_PLVAP_SPARCL1",
    "Contractile_VSMC_pericyte_like",
    "SPP1_lipid_associated_macrophage",
    "Neutrophil_like",
    "cDC2_like",
    "CD1C_immune_mixed_dataset_enriched",
    "Plasma_cell",
    "Mast_cell",
    "Endothelial_EC2_ITLN1_EDN1",
    "pDC_like",
    "Cycling_cells",
    "Interferon_stimulated_T",
    "Inflammatory_monocyte_dataset_enriched",
    "Macrophage_monocyte",
    "B_cell_dataset_enriched",
    "cDC1_like",
    "MHCII_B_cell_like_small_cluster",
    "Naive_B_cell_like_small_cluster",
    "Plasma_cell_small_cluster"
  ),
  confidence = c(
    "high", "high", "high", "high", "high",
    "high", "high", "high", "moderate", "high",
    "moderate", "high", "high", "high", "high",
    "high", "low", "high", "high", "high",
    "high", "moderate", "high", "moderate", "high",
    "moderate", "high", "moderate", "moderate", "high"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  annotation_map,
  file.path(out_dir, "final_whole_atlas_cluster_annotation.csv"),
  row.names = FALSE
)

obj$whole_atlas_cluster <- as.character(obj$seurat_clusters)

obj$whole_atlas_major_celltype <- annotation_map$major_celltype[
  match(
    obj$whole_atlas_cluster,
    annotation_map$cluster
  )
]

# =========================================================
# 4. Extract EC
# =========================================================
ec_source_clusters <- c("11", "19")

ec_cells <- Cells(obj)[
  obj$whole_atlas_cluster %in% ec_source_clusters
]

if (length(ec_cells) < 500) {
  stop("Too few EC candidate cells found. Expected clusters 11 and 19.")
}

ec <- subset(obj, cells = ec_cells)

cat("EC extraction completed.\n")
cat("EC cells:", ncol(ec), "\n")
cat(
  "Source clusters:",
  paste(sort(unique(ec$whole_atlas_cluster)), collapse = ", "),
  "\n"
)

# =========================================================
# 5. Save EC source composition before reclustering
# =========================================================
source_cluster_counts <- table(
  ec$whole_atlas_cluster,
  ec$dataset_id
)

source_cluster_props <- prop.table(
  source_cluster_counts,
  margin = 1
)

write.csv(
  as.data.frame.matrix(source_cluster_counts),
  file.path(out_dir, "EC_source_cluster_by_dataset_counts.csv")
)

write.csv(
  as.data.frame.matrix(source_cluster_props),
  file.path(out_dir, "EC_source_cluster_by_dataset_proportions.csv")
)

source_sample_counts <- table(
  ec$whole_atlas_cluster,
  ec$sample_id
)

write.csv(
  as.data.frame.matrix(source_sample_counts),
  file.path(out_dir, "EC_source_cluster_by_sample_counts.csv")
)

# =========================================================
# 6. Prepare RNA assay
# =========================================================
DefaultAssay(ec) <- "RNA"

rna_layers <- Layers(ec[["RNA"]])
count_layers <- grep("^counts", rna_layers, value = TRUE)

if (length(count_layers) == 0) {
  stop("No RNA counts layer found in extracted EC object.")
}

if (length(count_layers) > 1 || !"counts" %in% rna_layers) {
  ec <- JoinLayers(ec, assay = "RNA")
}

if (!"counts" %in% Layers(ec[["RNA"]])) {
  stop("RNA counts layer unavailable after JoinLayers.")
}

# =========================================================
# 7. EC-specific RNA preprocessing
# =========================================================
ec[["percent.mt"]] <- PercentageFeatureSet(
  ec,
  pattern = "^MT-"
)

cat("NormalizeData for EC...\n")

ec <- NormalizeData(
  ec,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

cat("FindVariableFeatures for EC...\n")

ec <- FindVariableFeatures(
  ec,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

cat("ScaleData for EC...\n")

ec <- ScaleData(
  ec,
  features = VariableFeatures(ec),
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

cat("RunPCA for EC...\n")

ec <- RunPCA(
  ec,
  features = VariableFeatures(ec),
  npcs = 40,
  verbose = FALSE
)

# =========================================================
# 8. Harmony correction within EC subset
# =========================================================
cat("RunHarmony for EC using dataset_id...\n")

ec <- harmony::RunHarmony(
  object = ec,
  group.by.vars = "dataset_id",
  reduction.use = "pca",
  dims.use = 1:30,
  reduction.save = "harmony",
  plot_convergence = FALSE,
  verbose = TRUE
)

# =========================================================
# 9. EC UMAP and reclustering
# =========================================================
cat("RunUMAP for EC...\n")

ec <- RunUMAP(
  ec,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap",
  reduction.key = "ECUMAP_",
  n.neighbors = 25,
  min.dist = 0.25,
  seed.use = 20260623,
  verbose = TRUE
)

cat("FindNeighbors for EC...\n")

ec <- FindNeighbors(
  ec,
  reduction = "harmony",
  dims = 1:30,
  verbose = TRUE
)

cat("FindClusters for EC...\n")

ec <- FindClusters(
  ec,
  resolution = 0.40,
  algorithm = 1,
  random.seed = 20260623,
  verbose = TRUE
)

ec$EC_cluster <- as.character(ec$seurat_clusters)
Idents(ec) <- "EC_cluster"

# =========================================================
# 10. EC composition by dataset and source cluster
# =========================================================
ec_cluster_dataset_counts <- table(
  ec$EC_cluster,
  ec$dataset_id
)

ec_cluster_dataset_props <- prop.table(
  ec_cluster_dataset_counts,
  margin = 1
)

ec_cluster_source_counts <- table(
  ec$EC_cluster,
  ec$whole_atlas_cluster
)

ec_cluster_sample_counts <- table(
  ec$EC_cluster,
  ec$sample_id
)

write.csv(
  as.data.frame.matrix(ec_cluster_dataset_counts),
  file.path(out_dir, "EC_cluster_by_dataset_counts.csv")
)

write.csv(
  as.data.frame.matrix(ec_cluster_dataset_props),
  file.path(out_dir, "EC_cluster_by_dataset_proportions.csv")
)

write.csv(
  as.data.frame.matrix(ec_cluster_source_counts),
  file.path(out_dir, "EC_cluster_by_whole_atlas_source_counts.csv")
)

write.csv(
  as.data.frame.matrix(ec_cluster_sample_counts),
  file.path(out_dir, "EC_cluster_by_sample_counts.csv")
)

ec_cluster_size <- as.data.frame(table(ec$EC_cluster))
colnames(ec_cluster_size) <- c("EC_cluster", "n_cells")

write.csv(
  ec_cluster_size,
  file.path(out_dir, "EC_cluster_cell_numbers.csv"),
  row.names = FALSE
)

# =========================================================
# 11. Lightweight EC marker enrichment
# =========================================================
cat("Calculating EC cluster-average RNA expression...\n")

rna_data <- LayerData(
  object = ec,
  assay = "RNA",
  layer = "data"
)

ec_cluster_levels <- as.character(
  sort(unique(as.numeric(as.character(ec$EC_cluster))))
)

ec_factor <- factor(
  as.character(ec$EC_cluster),
  levels = ec_cluster_levels
)

cluster_n <- as.numeric(table(ec_factor))
n_clusters <- length(ec_cluster_levels)

group_matrix <- Matrix::sparseMatrix(
  i = seq_len(ncol(rna_data)),
  j = as.integer(ec_factor),
  x = 1,
  dims = c(ncol(rna_data), n_clusters)
)

avg_expr <- as.matrix(rna_data %*% group_matrix)
avg_expr <- sweep(avg_expr, 2, cluster_n, FUN = "/")

rownames(avg_expr) <- rownames(rna_data)
colnames(avg_expr) <- ec_cluster_levels

sum_expr <- sweep(avg_expr, 2, cluster_n, FUN = "*")
total_cells <- sum(cluster_n)

exclude_pattern <- "^(MT-|RPL|RPS|MALAT1$|HBA|HBB)"

marker_list <- vector("list", n_clusters)

for (j in seq_len(n_clusters)) {

  other_cells <- total_cells - cluster_n[j]

  other_mean <- (
    rowSums(sum_expr) - sum_expr[, j]
  ) / max(other_cells, 1)

  enrichment_score <- avg_expr[, j] - other_mean

  valid <- (
    avg_expr[, j] >= 0.10 &
    !grepl(exclude_pattern, rownames(avg_expr), ignore.case = TRUE)
  )

  idx <- which(valid)

  idx <- idx[
    order(enrichment_score[idx], decreasing = TRUE)
  ]

  idx <- head(idx, 30)

  marker_list[[j]] <- data.frame(
    EC_cluster = ec_cluster_levels[j],
    gene = rownames(avg_expr)[idx],
    avg_expression = avg_expr[idx, j],
    mean_expression_other_clusters = other_mean[idx],
    enrichment_score = enrichment_score[idx],
    stringsAsFactors = FALSE
  )
}

all_ec_markers <- do.call(rbind, marker_list)

top20_ec_markers <- do.call(
  rbind,
  lapply(
    split(all_ec_markers, all_ec_markers$EC_cluster),
    head,
    20
  )
)

write.csv(
  all_ec_markers,
  file.path(out_dir, "all_EC_cluster_enriched_genes.csv"),
  row.names = FALSE
)

write.csv(
  top20_ec_markers,
  file.path(out_dir, "top20_markers_each_EC_cluster.csv"),
  row.names = FALSE
)

# =========================================================
# 12. EC-state marker DotPlot
# =========================================================
ec_state_markers <- c(
  # Pan-EC
  "PECAM1", "VWF", "KDR", "CDH5", "CLDN5", "RAMP2", "EMCN",

  # Angiogenic / capillary-like
  "PLVAP", "RGCC", "CA4", "SPARCL1", "EGFL7", "KDR",

  # Arterial / quiescent
  "EFNB2", "SOX17", "GJA5", "HEY1", "NOTCH4", "EDN1",

  # Activated / inflammatory
  "ICAM1", "VCAM1", "SELE", "SELP", "CXCL8", "CXCL12",

  # ACKR1 / venous-like
  "ACKR1",

  # EndMT-like
  "VIM", "TAGLN", "ACTA2", "COL1A1", "COL1A2",
  "COL3A1", "SNAI1", "SNAI2", "TWIST1", "ZEB1",

  # IFN / stress
  "ISG15", "IFI27", "IFITM3"
)

ec_state_markers <- unique(
  intersect(ec_state_markers, rownames(rna_data))
)

if (length(ec_state_markers) > 0) {

  marker_expr <- avg_expr[ec_state_markers, , drop = FALSE]

  marker_detected <- Matrix::Matrix(
    rna_data[ec_state_markers, , drop = FALSE] > 0,
    sparse = TRUE
  )

  pct_expr <- as.matrix(marker_detected %*% group_matrix)
  pct_expr <- sweep(pct_expr, 2, cluster_n, FUN = "/")

  scaled_expr <- t(scale(t(marker_expr)))
  scaled_expr[is.na(scaled_expr)] <- 0

  dot_df <- data.frame(
    gene = rep(rownames(marker_expr), times = ncol(marker_expr)),
    EC_cluster = rep(colnames(marker_expr), each = nrow(marker_expr)),
    scaled_avg_expression = as.vector(scaled_expr),
    pct_expression = as.vector(pct_expr),
    stringsAsFactors = FALSE
  )

  dot_df$EC_cluster <- factor(
    dot_df$EC_cluster,
    levels = rev(ec_cluster_levels)
  )

  p_dot <- ggplot(
    dot_df,
    aes(
      x = gene,
      y = EC_cluster,
      size = pct_expression,
      colour = scaled_avg_expression
    )
  ) +
    geom_point(alpha = 0.90) +
    scale_size(range = c(0.2, 7)) +
    labs(
      title = "EC reclustering: state-marker overview",
      subtitle = "Colour: scaled average RNA expression; size: expressing-cell proportion",
      x = "EC state markers",
      y = "EC cluster"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 55, hjust = 1)
    )

  ggsave(
    file.path(out_dir, "DotPlot_EC_state_markers.png"),
    p_dot,
    width = 22,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}

# =========================================================
# 13. EC UMAP figures
# =========================================================
p_ec_cluster <- DimPlot(
  ec,
  reduction = "umap",
  group.by = "EC_cluster",
  label = TRUE,
  label.size = 5,
  pt.size = 0.15
) +
  ggtitle("AS endothelial cells: Harmony UMAP reclustering") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(out_dir, "UMAP_EC_clusters.png"),
  p_ec_cluster,
  width = 12,
  height = 9,
  dpi = 300,
  bg = "white"
)

p_ec_source <- DimPlot(
  ec,
  reduction = "umap",
  group.by = "whole_atlas_cluster",
  pt.size = 0.15
) +
  ggtitle("EC subset coloured by original whole-atlas cluster") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(out_dir, "UMAP_EC_by_original_cluster.png"),
  p_ec_source,
  width = 12,
  height = 9,
  dpi = 300,
  bg = "white"
)

p_ec_dataset <- DimPlot(
  ec,
  reduction = "umap",
  group.by = "dataset_id",
  pt.size = 0.15
) +
  ggtitle("EC subset coloured by dataset") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(out_dir, "UMAP_EC_by_dataset.png"),
  p_ec_dataset,
  width = 12,
  height = 9,
  dpi = 300,
  bg = "white"
)

# =========================================================
# 14. Save EC object and report
# =========================================================
saveRDS(ec, ec_final_rds)

report_lines <- c(
  "========== EC EXTRACTION AND RECLUSTERING SUMMARY ==========",
  paste("Whole-atlas input:", input_rds),
  paste("EC output object:", ec_final_rds),
  paste("Whole-atlas EC source clusters:", paste(ec_source_clusters, collapse = ", ")),
  paste("EC cells extracted:", ncol(ec)),
  paste("EC reclusters:", length(unique(ec$EC_cluster))),
  "EC preprocessing: LogNormalize + 2000 HVGs + percent.mt regression.",
  "EC Harmony batch variable: dataset_id.",
  "EC Harmony dimensions: 1-30.",
  "EC clustering resolution: 0.40.",
  "============================================================="
)

writeLines(
  report_lines,
  file.path(out_dir, "run_summary.txt")
)

cat(paste(report_lines, collapse = "\n"), "\n")
cat("\nStep 10 finished successfully.\n")
