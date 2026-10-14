#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 16 * 1024^3)

required_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "harmony",
  "ggplot2",
  "Matrix"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nThis script does not install packages automatically."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(ggplot2)
  library(Matrix)
})

set.seed(20260729)

in_rds <- paste0(
  "/root/AS_external_results/GSE159677/stage3_doublet/",
  "GSE159677_qc_doublet_clean_singlets.rds"
)

out_dir <- "/root/AS_external_results/GSE159677/stage4_fullcell_harmony"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "objects"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "00_stage4_log.txt")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("===== GSE159677 STAGE 4: FULL-CELL HARMONY AND PRELIMINARY ANNOTATION =====\n")
cat("Start:", format(Sys.time()), "\n")
cat("Input:", in_rds, "\n")
cat("Output:", out_dir, "\n")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("Harmony:", as.character(packageVersion("harmony")), "\n\n")

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

DefaultAssay(obj) <- "RNA"

# Seurat v5 merged objects can contain multiple counts layers.
if (inherits(obj[["RNA"]], "Assay5")) {
  current_layers <- SeuratObject::Layers(obj[["RNA"]])
  cat("RNA layers before JoinLayers:", paste(current_layers, collapse = ", "), "\n")

  count_layers <- grep("^counts", current_layers, value = TRUE)
  if (length(count_layers) > 1) {
    cat("Joining multiple RNA layers...\n")
    obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
  }

  cat(
    "RNA layers after JoinLayers:",
    paste(SeuratObject::Layers(obj[["RNA"]]), collapse = ", "),
    "\n"
  )
}

cat("\nCells:", ncol(obj), "\n")
cat("Genes:", nrow(obj), "\n")
cat("Samples:", paste(sort(unique(obj$sample_id)), collapse = ", "), "\n\n")

# Preserve all biological labels. Region is not regressed out.
cat("Normalizing RNA counts...\n")
obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

cat("Selecting variable genes...\n")
obj <- FindVariableFeatures(
  obj,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 3000,
  verbose = TRUE
)

hvg <- VariableFeatures(obj)
writeLines(
  hvg,
  file.path(out_dir, "tables", "01_variable_features_3000.txt")
)

cat("Scaling variable genes and regressing mitochondrial percentage...\n")
obj <- ScaleData(
  obj,
  assay = "RNA",
  features = hvg,
  vars.to.regress = "percent.mt",
  verbose = TRUE
)

cat("Running PCA...\n")
obj <- RunPCA(
  obj,
  assay = "RNA",
  features = hvg,
  npcs = 50,
  verbose = TRUE
)

p_elbow <- ElbowPlot(obj, ndims = 50) +
  ggtitle("GSE159677 PCA elbow plot")

ggsave(
  file.path(out_dir, "plots", "01_PCA_elbow.png"),
  p_elbow,
  width = 7,
  height = 5,
  dpi = 300
)

cat("Running Harmony by biological sample...\n")
obj <- harmony::RunHarmony(
  object = obj,
  group.by.vars = "sample_id",
  reduction.use = "pca",
  dims.use = 1:40,
  reduction.save = "harmony",
  project.dim = FALSE,
  max.iter.harmony = 30,
  verbose = TRUE
)

if (!("harmony" %in% Reductions(obj))) {
  stop("Harmony reduction was not created.")
}

cat("Running UMAP and graph clustering...\n")
dims_use <- 1:30

obj <- FindNeighbors(
  obj,
  reduction = "harmony",
  dims = dims_use,
  k.param = 20,
  verbose = TRUE
)

obj <- FindClusters(
  obj,
  resolution = c(0.2, 0.4, 0.6, 0.8),
  algorithm = 1,
  random.seed = 20260729,
  verbose = TRUE
)

final_resolution_col <- "RNA_snn_res.0.6"
if (!(final_resolution_col %in% colnames(obj[[]]))) {
  candidate_cols <- grep(
    "snn_res",
    colnames(obj[[]]),
    value = TRUE
  )
  stop(
    "Expected clustering column not found: ",
    final_resolution_col,
    "\nAvailable resolution columns: ",
    paste(candidate_cols, collapse = ", ")
  )
}

Idents(obj) <- final_resolution_col
obj$seurat_clusters <- as.character(Idents(obj))

obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = dims_use,
  n.neighbors = 30,
  min.dist = 0.30,
  spread = 1,
  seed.use = 20260729,
  reduction.name = "umap_harmony",
  reduction.key = "hUMAP_",
  verbose = TRUE
)

# Standard plots
save_dimplot <- function(group_by, filename, label = FALSE, split_by = NULL) {
  p <- DimPlot(
    obj,
    reduction = "umap_harmony",
    group.by = group_by,
    split.by = split_by,
    label = label,
    repel = TRUE,
    raster = TRUE
  ) +
    ggtitle(paste("GSE159677:", group_by))
  ggsave(
    file.path(out_dir, "plots", filename),
    p,
    width = if (is.null(split_by)) 8.5 else 13,
    height = 6.5,
    dpi = 300
  )
}

save_dimplot("seurat_clusters", "02_UMAP_clusters_res0.6.png", label = TRUE)
save_dimplot("sample_id", "03_UMAP_sample_id.png")
save_dimplot("patient_id", "04_UMAP_patient_id.png")
save_dimplot("region", "05_UMAP_region_AC_PA.png")
save_dimplot("seurat_clusters", "06_UMAP_clusters_split_by_region.png", label = TRUE, split_by = "region")

# Major-lineage marker panels
marker_sets <- list(
  Endothelial = c("PECAM1", "VWF", "KDR", "ESAM", "EMCN", "CLDN5", "RAMP2", "ENG", "CDH5"),
  Smooth_muscle = c("ACTA2", "TAGLN", "MYH11", "CNN1", "LMOD1", "MYL9"),
  Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "C7"),
  Pericyte = c("RGS5", "CSPG4", "PDGFRB", "MCAM", "NOTCH3", "DES"),
  Myeloid = c("LST1", "TYROBP", "FCER1G", "CTSS", "AIF1", "CD68", "CSF1R"),
  T_cell = c("CD3D", "CD3E", "TRAC", "IL7R", "LTB", "MALAT1"),
  NK_cell = c("NKG7", "GNLY", "KLRD1", "PRF1", "GZMB"),
  B_cell = c("MS4A1", "CD79A", "CD79B", "CD37", "CD74", "HLA-DRA"),
  Mast_cell = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2")
)

present_marker_sets <- lapply(
  marker_sets,
  function(x) intersect(x, rownames(obj))
)

present_marker_sets <- present_marker_sets[
  lengths(present_marker_sets) >= 2
]

if (length(present_marker_sets) < 5) {
  stop(
    "Too few major-lineage marker sets were detected. Present sets: ",
    paste(names(present_marker_sets), collapse = ", ")
  )
}

all_dot_genes <- unique(unlist(present_marker_sets, use.names = FALSE))

p_dot <- DotPlot(
  obj,
  features = all_dot_genes,
  group.by = "seurat_clusters",
  assay = "RNA",
  dot.scale = 5
) +
  RotatedAxis() +
  ggtitle("GSE159677 major cell-type markers")

ggsave(
  file.path(out_dir, "plots", "07_DotPlot_major_celltype_markers.png"),
  p_dot,
  width = 18,
  height = 7.5,
  dpi = 300,
  limitsize = FALSE
)

feature_genes <- intersect(
  c(
    "PECAM1", "VWF", "CLDN5",
    "ACTA2", "MYH11",
    "COL1A1", "DCN",
    "LST1", "CD3D", "NKG7", "MS4A1", "TPSAB1"
  ),
  rownames(obj)
)

p_feat <- FeaturePlot(
  obj,
  features = feature_genes,
  reduction = "umap_harmony",
  ncol = 3,
  raster = TRUE,
  order = TRUE
)

ggsave(
  file.path(out_dir, "plots", "08_FeaturePlot_major_markers.png"),
  p_feat,
  width = 12,
  height = 4 * ceiling(length(feature_genes) / 3),
  dpi = 300,
  limitsize = FALSE
)

# Preliminary program scores based on mean normalized log-expression.
cat("Calculating preliminary major-lineage scores...\n")
data_mat <- SeuratObject::LayerData(
  obj,
  assay = "RNA",
  layer = "data"
)

score_cols <- character(0)
for (nm in names(present_marker_sets)) {
  genes <- present_marker_sets[[nm]]
  score_name <- paste0(nm, "_Score")
  obj[[score_name]] <- Matrix::colMeans(
    data_mat[genes, , drop = FALSE]
  )
  score_cols <- c(score_cols, score_name)
}

meta <- obj[[]]
cluster_score_means <- aggregate(
  meta[, score_cols, drop = FALSE],
  by = list(cluster = as.character(meta$seurat_clusters)),
  FUN = mean
)

score_matrix <- as.matrix(
  cluster_score_means[, score_cols, drop = FALSE]
)

top_index <- apply(score_matrix, 1, which.max)
top_score <- score_matrix[cbind(seq_len(nrow(score_matrix)), top_index)]

second_score <- apply(
  score_matrix,
  1,
  function(x) sort(x, decreasing = TRUE)[min(2, length(x))]
)

cluster_score_means$preliminary_major_celltype <- sub(
  "_Score$",
  "",
  score_cols[top_index]
)
cluster_score_means$top_score <- top_score
cluster_score_means$top_minus_second <- top_score - second_score

write.csv(
  cluster_score_means,
  file.path(out_dir, "tables", "02_cluster_preliminary_major_annotation.csv"),
  row.names = FALSE
)

cluster_to_label <- setNames(
  cluster_score_means$preliminary_major_celltype,
  cluster_score_means$cluster
)

obj$preliminary_major_celltype <- unname(
  cluster_to_label[as.character(obj$seurat_clusters)]
)

p_prelim <- DimPlot(
  obj,
  reduction = "umap_harmony",
  group.by = "preliminary_major_celltype",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle("GSE159677 preliminary major cell-type annotation")

ggsave(
  file.path(out_dir, "plots", "09_UMAP_preliminary_major_celltype.png"),
  p_prelim,
  width = 9,
  height = 6.5,
  dpi = 300
)

# Composition summaries
cluster_sample <- as.data.frame(
  table(
    cluster = obj$seurat_clusters,
    sample_id = obj$sample_id
  )
)
write.csv(
  cluster_sample,
  file.path(out_dir, "tables", "03_cluster_by_sample_counts.csv"),
  row.names = FALSE
)

cluster_region <- as.data.frame(
  table(
    cluster = obj$seurat_clusters,
    region = obj$region
  )
)
write.csv(
  cluster_region,
  file.path(out_dir, "tables", "04_cluster_by_region_counts.csv"),
  row.names = FALSE
)

major_sample <- as.data.frame(
  table(
    preliminary_major_celltype = obj$preliminary_major_celltype,
    sample_id = obj$sample_id
  )
)
write.csv(
  major_sample,
  file.path(out_dir, "tables", "05_major_celltype_by_sample_counts.csv"),
  row.names = FALSE
)

major_region <- as.data.frame(
  table(
    preliminary_major_celltype = obj$preliminary_major_celltype,
    region = obj$region
  )
)
write.csv(
  major_region,
  file.path(out_dir, "tables", "06_major_celltype_by_region_counts.csv"),
  row.names = FALSE
)

# Stacked sample composition plot
major_sample$sample_total <- ave(
  major_sample$Freq,
  major_sample$sample_id,
  FUN = sum
)
major_sample$fraction <- major_sample$Freq / major_sample$sample_total

p_comp <- ggplot(
  major_sample,
  aes(
    x = sample_id,
    y = fraction,
    fill = preliminary_major_celltype
  )
) +
  geom_col() +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(
    x = NULL,
    y = "Cell fraction",
    fill = "Preliminary type",
    title = "GSE159677 preliminary cell-type composition by sample"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(out_dir, "plots", "10_major_celltype_fraction_by_sample.png"),
  p_comp,
  width = 10,
  height = 6,
  dpi = 300
)

# Cluster markers. Use capped cells per cluster for speed and reproducibility.
cat("Finding positive markers for each cluster...\n")
Idents(obj) <- "seurat_clusters"

markers <- FindAllMarkers(
  obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25,
  max.cells.per.ident = 3000,
  random.seed = 20260729,
  verbose = TRUE
)

write.csv(
  markers,
  file.path(out_dir, "tables", "07_all_positive_cluster_markers.csv"),
  row.names = FALSE
)

if (nrow(markers) > 0) {
  fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(markers))[1]
  if (is.na(fc_col)) {
    stop("Could not identify the fold-change column in marker results.")
  }

  split_markers <- split(markers, markers$cluster)
  top10 <- do.call(
    rbind,
    lapply(
      split_markers,
      function(df) {
        ord <- order(df$p_val_adj, -df[[fc_col]], na.last = TRUE)
        head(df[ord, , drop = FALSE], 10)
      }
    )
  )
  rownames(top10) <- NULL

  write.csv(
    top10,
    file.path(out_dir, "tables", "08_top10_markers_each_cluster.csv"),
    row.names = FALSE
  )
}

# Save formal full-cell object.
saveRDS(
  obj,
  file.path(out_dir, "objects", "GSE159677_fullcell_harmony_clustered_preliminary_annotated.rds"),
  compress = FALSE
)

# Candidate endothelial clusters only; this is NOT yet the final high-confidence EC set.
candidate_ec_clusters <- cluster_score_means$cluster[
  cluster_score_means$preliminary_major_celltype == "Endothelial"
]

writeLines(
  candidate_ec_clusters,
  file.path(out_dir, "tables", "09_candidate_endothelial_clusters.txt")
)

if (length(candidate_ec_clusters) > 0) {
  candidate_ec_cells <- colnames(obj)[
    as.character(obj$seurat_clusters) %in% candidate_ec_clusters
  ]

  candidate_ec_obj <- subset(obj, cells = candidate_ec_cells)

  saveRDS(
    candidate_ec_obj,
    file.path(out_dir, "objects", "GSE159677_candidate_EC_before_manual_review.rds"),
    compress = FALSE
  )

  writeLines(
    candidate_ec_cells,
    file.path(out_dir, "tables", "10_candidate_EC_cell_barcodes.txt")
  )

  cat(
    "\nCandidate endothelial clusters:",
    paste(candidate_ec_clusters, collapse = ", "),
    "\n"
  )
  cat("Candidate endothelial cells:", length(candidate_ec_cells), "\n")
} else {
  cat("\nWARNING: no cluster received Endothelial as its top preliminary program.\n")
}

cat("\nFinal cluster count:", length(unique(obj$seurat_clusters)), "\n")
cat("Saved formal full-cell object.\n")
cat("IMPORTANT: preliminary labels must be manually reviewed before final EC extraction.\n")
cat("\n===== STAGE 4 COMPLETE =====\n")
cat("End:", format(Sys.time()), "\n")
