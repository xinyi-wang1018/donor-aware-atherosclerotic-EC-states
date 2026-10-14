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

input_candidates <- c(
  "/root/AS_external_results/GSE159677/stage4_fullcell_harmony/objects/GSE159677_fullcell_harmony_clustered_preliminary_annotated.rds",
  "/root/AS_external_results/GSE159677/stage4_fullcell_harmony/objects/GSE159677_checkpoint_after_UMAP_clustering.rds"
)

existing_inputs <- input_candidates[file.exists(input_candidates)]
if (length(existing_inputs) == 0) {
  stop(
    "No usable Stage 4 object was found.\nChecked:\n",
    paste(input_candidates, collapse = "\n")
  )
}
in_rds <- existing_inputs[1]

out_dir <- "/root/AS_external_results/GSE159677/stage5_manual_annotation_EC"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "objects"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "00_stage5_log.txt")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("===== GSE159677 STAGE 5: MANUAL MAJOR ANNOTATION AND EC REANALYSIS =====\n")
cat("Start:", format(Sys.time()), "\n")
cat("Input:", in_rds, "\n")
cat("Output:", out_dir, "\n")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("Harmony:", as.character(packageVersion("harmony")), "\n\n")

obj <- readRDS(in_rds)
DefaultAssay(obj) <- "RNA"

if (!("umap_harmony" %in% Reductions(obj))) {
  stop("The input object does not contain the umap_harmony reduction.")
}

if (!("seurat_clusters" %in% colnames(obj[[]]))) {
  candidate_cols <- grep("snn_res", colnames(obj[[]]), value = TRUE)
  if ("RNA_snn_res.0.6" %in% candidate_cols) {
    obj$seurat_clusters <- as.character(obj$RNA_snn_res.0.6)
  } else {
    stop(
      "No seurat_clusters metadata column was found. Available clustering columns: ",
      paste(candidate_cols, collapse = ", ")
    )
  }
}

obj$seurat_clusters <- as.character(obj$seurat_clusters)

# Manual review based on the uploaded DotPlot and cluster top markers.
manual_map <- c(
  "0"  = "T_NK_cytotoxic",
  "1"  = "Fibroblast",
  "2"  = "Endothelial_ACKR1_venous",
  "3"  = "T_CD4_IL7R",
  "4"  = "T_CCR7_naive_memory",
  "5"  = "Macrophage_FOLR2",
  "6"  = "SMC_contractile",
  "7"  = "Macrophage_SPP1",
  "8"  = "Immune_stressed_highMT",
  "9"  = "Monocyte_FCN1",
  "10" = "B_cell",
  "11" = "NK_cytotoxic",
  "12" = "T_NK_gammadelta",
  "13" = "Endothelial_ITLN1_arterial",
  "14" = "SMC_modulated",
  "15" = "Pericyte",
  "16" = "Dendritic_CD1C",
  "17" = "Mast_cell",
  "18" = "Endothelial_SOX17_arterial",
  "19" = "Plasma_cell",
  "20" = "Cycling_immune",
  "21" = "Macrophage_C1QC",
  "22" = "Dendritic_plasmacytoid"
)

observed_clusters <- sort(unique(obj$seurat_clusters))
missing_in_map <- setdiff(observed_clusters, names(manual_map))
if (length(missing_in_map) > 0) {
  stop(
    "The manual map is missing cluster(s): ",
    paste(missing_in_map, collapse = ", ")
  )
}

obj$major_celltype_manual <- unname(manual_map[obj$seurat_clusters])

major_broad_map <- c(
  T_NK_cytotoxic = "T_NK",
  Fibroblast = "Fibroblast",
  Endothelial_ACKR1_venous = "Endothelial",
  T_CD4_IL7R = "T_NK",
  T_CCR7_naive_memory = "T_NK",
  Macrophage_FOLR2 = "Myeloid",
  SMC_contractile = "Smooth_muscle",
  Macrophage_SPP1 = "Myeloid",
  Immune_stressed_highMT = "Immune_stressed",
  Monocyte_FCN1 = "Myeloid",
  B_cell = "B_cell",
  NK_cytotoxic = "T_NK",
  T_NK_gammadelta = "T_NK",
  Endothelial_ITLN1_arterial = "Endothelial",
  SMC_modulated = "Smooth_muscle",
  Pericyte = "Pericyte",
  Dendritic_CD1C = "Myeloid",
  Mast_cell = "Mast_cell",
  Endothelial_SOX17_arterial = "Endothelial",
  Plasma_cell = "Plasma_cell",
  Cycling_immune = "Cycling_immune",
  Macrophage_C1QC = "Myeloid",
  Dendritic_plasmacytoid = "Myeloid"
)

obj$major_celltype_broad_manual <- unname(
  major_broad_map[obj$major_celltype_manual]
)

write.csv(
  data.frame(
    cluster = names(manual_map),
    manual_annotation = unname(manual_map),
    broad_annotation = unname(
      major_broad_map[unname(manual_map)]
    ),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "tables", "01_manual_cluster_annotation.csv"),
  row.names = FALSE
)

p_manual <- DimPlot(
  obj,
  reduction = "umap_harmony",
  group.by = "major_celltype_manual",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggtitle("GSE159677 manually reviewed major cell types")

ggsave(
  file.path(out_dir, "plots", "01_UMAP_manual_major_celltype.png"),
  p_manual,
  width = 12,
  height = 7.5,
  dpi = 300
)

p_broad <- DimPlot(
  obj,
  reduction = "umap_harmony",
  group.by = "major_celltype_broad_manual",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggtitle("GSE159677 broad manually reviewed cell types")

ggsave(
  file.path(out_dir, "plots", "02_UMAP_manual_broad_celltype.png"),
  p_broad,
  width = 10,
  height = 7,
  dpi = 300
)

saveRDS(
  obj,
  file.path(out_dir, "objects", "GSE159677_fullcell_manual_annotation.rds"),
  compress = FALSE
)

# High-confidence endothelial clusters identified by pan-endothelial markers.
high_confidence_ec_clusters <- c("2", "13", "18")
ec_cells <- colnames(obj)[
  obj$seurat_clusters %in% high_confidence_ec_clusters
]

if (length(ec_cells) < 200) {
  stop(
    "Too few endothelial cells were extracted: ",
    length(ec_cells)
  )
}

ec <- subset(obj, cells = ec_cells)
DefaultAssay(ec) <- "RNA"

ec$parent_cluster <- as.character(ec$seurat_clusters)
ec$parent_EC_state <- unname(
  c(
    "2" = "ACKR1_venous_EC",
    "13" = "ITLN1_arterial_EC",
    "18" = "SOX17_arterial_EC"
  )[ec$parent_cluster]
)

cat("\nHigh-confidence EC clusters:", paste(high_confidence_ec_clusters, collapse = ", "), "\n")
cat("High-confidence EC cells:", ncol(ec), "\n")

ec_composition <- as.data.frame(
  table(
    parent_EC_state = ec$parent_EC_state,
    sample_id = ec$sample_id,
    region = ec$region
  )
)
write.csv(
  ec_composition,
  file.path(out_dir, "tables", "02_EC_parent_state_by_sample_region_counts.csv"),
  row.names = FALSE
)

# Join Seurat v5 layers after subsetting/merging.
if (inherits(ec[["RNA"]], "Assay5")) {
  current_layers <- SeuratObject::Layers(ec[["RNA"]])
  count_layers <- grep("^counts", current_layers, value = TRUE)
  if (length(count_layers) > 1) {
    ec <- SeuratObject::JoinLayers(ec, assay = "RNA")
  }
}

cat("\nReprocessing endothelial cells...\n")
ec <- NormalizeData(
  ec,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

ec <- FindVariableFeatures(
  ec,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = min(2500, nrow(ec)),
  verbose = TRUE
)

hvg <- VariableFeatures(ec)

ec <- ScaleData(
  ec,
  assay = "RNA",
  features = hvg,
  vars.to.regress = "percent.mt",
  verbose = TRUE
)

ec <- RunPCA(
  ec,
  assay = "RNA",
  features = hvg,
  npcs = 40,
  verbose = TRUE
)

saveRDS(
  ec,
  file.path(out_dir, "objects", "GSE159677_EC_checkpoint_after_PCA.rds"),
  compress = FALSE
)

ec <- harmony::RunHarmony(
  object = ec,
  group.by.vars = "sample_id",
  reduction.use = "pca",
  dims.use = 1:30,
  reduction.save = "harmony_ec",
  project.dim = FALSE,
  max_iter = 30,
  early_stop = TRUE,
  verbose = TRUE
)

if (!("harmony_ec" %in% Reductions(ec))) {
  stop("The endothelial Harmony reduction was not created.")
}

saveRDS(
  ec,
  file.path(out_dir, "objects", "GSE159677_EC_checkpoint_after_Harmony.rds"),
  compress = FALSE
)

dims_use <- 1:min(25, ncol(Embeddings(ec, "harmony_ec")))

ec <- FindNeighbors(
  ec,
  reduction = "harmony_ec",
  dims = dims_use,
  k.param = min(20, max(5, floor(sqrt(ncol(ec))))),
  verbose = TRUE
)

ec <- FindClusters(
  ec,
  resolution = c(0.2, 0.4, 0.6, 0.8),
  algorithm = 1,
  random.seed = 20260729,
  verbose = TRUE
)

final_col <- "RNA_snn_res.0.4"
if (!(final_col %in% colnames(ec[[]]))) {
  stop("Expected endothelial clustering column was not created: ", final_col)
}

Idents(ec) <- final_col
ec$EC_cluster <- as.character(Idents(ec))

ec <- RunUMAP(
  ec,
  reduction = "harmony_ec",
  dims = dims_use,
  n.neighbors = min(30, max(10, floor(sqrt(ncol(ec))))),
  min.dist = 0.25,
  spread = 1,
  seed.use = 20260729,
  reduction.name = "umap_ec",
  reduction.key = "ECUMAP_",
  verbose = TRUE
)

saveRDS(
  ec,
  file.path(out_dir, "objects", "GSE159677_EC_checkpoint_after_UMAP_clustering.rds"),
  compress = FALSE
)

save_ec_dimplot <- function(group_by, filename, label = FALSE, split_by = NULL) {
  p <- DimPlot(
    ec,
    reduction = "umap_ec",
    group.by = group_by,
    split.by = split_by,
    label = label,
    repel = TRUE,
    raster = FALSE
  ) +
    ggtitle(paste("GSE159677 EC:", group_by))

  ggsave(
    file.path(out_dir, "plots", filename),
    p,
    width = if (is.null(split_by)) 8.5 else 13,
    height = 6.5,
    dpi = 300
  )
}

save_ec_dimplot("EC_cluster", "03_EC_UMAP_clusters_res0.4.png", label = TRUE)
save_ec_dimplot("parent_EC_state", "04_EC_UMAP_parent_states.png", label = TRUE)
save_ec_dimplot("sample_id", "05_EC_UMAP_sample_id.png")
save_ec_dimplot("patient_id", "06_EC_UMAP_patient_id.png")
save_ec_dimplot("region", "07_EC_UMAP_region_AC_PA.png")
save_ec_dimplot("EC_cluster", "08_EC_UMAP_clusters_split_by_region.png", label = TRUE, split_by = "region")

# Endothelial identity and pathology programs.
programs <- list(
  Pan_endothelial = c(
    "PECAM1", "VWF", "KDR", "ESAM", "EMCN", "CLDN5",
    "RAMP2", "ENG", "CDH5", "CLEC14A", "ADGRL4", "EGFL7"
  ),
  Arterial_quiescent = c(
    "SOX17", "SOX18", "EFNB2", "GJA4", "GJA5", "HEY1",
    "DLL4", "BMX", "CXCL12", "F5", "ITLN1", "SELP"
  ),
  Capillary_angiogenic = c(
    "RGCC", "CA4", "EDNRB", "ESM1", "APLN", "PGF",
    "ANGPT2", "PLVAP", "FLT1", "KDR", "ADM"
  ),
  Venous_repair = c(
    "ACKR1", "NR2F2", "PLAT", "CCL14", "VWA1",
    "MMRN1", "RAMP3", "STC1", "MAFB", "POSTN"
  ),
  Inflammatory_activation = c(
    "SELE", "VCAM1", "ICAM1", "CCL2", "CXCL2",
    "IL6", "NFKBIA", "ISG15", "IFIT1", "IFI27"
  ),
  ECM_EndMT = c(
    "FN1", "BGN", "COL8A1", "ELN", "COL1A1", "COL1A2",
    "COL3A1", "SPARC", "POSTN", "VIM", "TAGLN", "ACTA2"
  ),
  Fibroblast_contamination = c(
    "DCN", "LUM", "COL1A1", "COL1A2", "COL3A1",
    "PDGFRA", "SFRP2", "FBLN1"
  ),
  SMC_pericyte_contamination = c(
    "ACTA2", "TAGLN", "MYH11", "CNN1",
    "RGS5", "PDGFRB", "CSPG4", "NOTCH3"
  )
)

programs_present <- lapply(
  programs,
  function(g) intersect(g, rownames(ec))
)
programs_present <- programs_present[lengths(programs_present) >= 3]

writeLines(
  unlist(
    lapply(
      names(programs_present),
      function(nm) {
        paste0(nm, "\t", paste(programs_present[[nm]], collapse = ","))
      }
    )
  ),
  file.path(out_dir, "tables", "03_EC_program_genes_present.txt")
)

data_mat <- SeuratObject::LayerData(
  ec,
  assay = "RNA",
  layer = "data"
)

score_cols <- character(0)
for (nm in names(programs_present)) {
  score_name <- paste0(nm, "_Score")
  ec[[score_name]] <- Matrix::colMeans(
    data_mat[programs_present[[nm]], , drop = FALSE]
  )
  score_cols <- c(score_cols, score_name)
}

ec_meta <- ec[[]]

cluster_score_mean <- aggregate(
  ec_meta[, score_cols, drop = FALSE],
  by = list(EC_cluster = ec_meta$EC_cluster),
  FUN = mean
)
write.csv(
  cluster_score_mean,
  file.path(out_dir, "tables", "04_EC_cluster_program_scores.csv"),
  row.names = FALSE
)

# Row-scaled program heatmap using base R only.
score_matrix <- as.matrix(
  cluster_score_mean[, score_cols, drop = FALSE]
)
rownames(score_matrix) <- cluster_score_mean$EC_cluster
score_matrix_scaled <- t(scale(t(score_matrix)))
score_matrix_scaled[is.na(score_matrix_scaled)] <- 0

png(
  file.path(out_dir, "plots", "09_EC_cluster_program_score_heatmap.png"),
  width = 2200,
  height = 1700,
  res = 220
)
heatmap(
  score_matrix_scaled,
  Rowv = NA,
  Colv = NA,
  scale = "none",
  margins = c(12, 14),
  xlab = "EC cluster",
  ylab = "Program",
  main = "GSE159677 EC program scores (row-scaled)"
)
dev.off()

dot_genes <- unique(unlist(programs_present, use.names = FALSE))
p_ec_dot <- DotPlot(
  ec,
  features = dot_genes,
  group.by = "EC_cluster",
  assay = "RNA",
  dot.scale = 5
) +
  RotatedAxis() +
  ggtitle("GSE159677 endothelial identity and pathology markers")

ggsave(
  file.path(out_dir, "plots", "10_EC_DotPlot_identity_pathology_markers.png"),
  p_ec_dot,
  width = 22,
  height = 7.5,
  dpi = 300,
  limitsize = FALSE
)

feature_genes <- intersect(
  c(
    "PECAM1", "VWF", "CLDN5",
    "SOX17", "GJA4", "CXCL12", "ITLN1",
    "ACKR1", "CCL14",
    "PGF", "ANGPT2", "PLVAP",
    "SELE", "VCAM1",
    "FN1", "BGN", "COL8A1", "POSTN"
  ),
  rownames(ec)
)

p_ec_feat <- FeaturePlot(
  ec,
  features = feature_genes,
  reduction = "umap_ec",
  ncol = 3,
  raster = FALSE,
  order = TRUE
)

ggsave(
  file.path(out_dir, "plots", "11_EC_FeaturePlot_key_markers.png"),
  p_ec_feat,
  width = 12,
  height = 4 * ceiling(length(feature_genes) / 3),
  dpi = 300,
  limitsize = FALSE
)

# EC cluster composition by sample/region.
ec_cluster_sample <- as.data.frame(
  table(
    EC_cluster = ec$EC_cluster,
    sample_id = ec$sample_id
  )
)
write.csv(
  ec_cluster_sample,
  file.path(out_dir, "tables", "05_EC_cluster_by_sample_counts.csv"),
  row.names = FALSE
)

ec_cluster_region <- as.data.frame(
  table(
    EC_cluster = ec$EC_cluster,
    region = ec$region
  )
)
write.csv(
  ec_cluster_region,
  file.path(out_dir, "tables", "06_EC_cluster_by_region_counts.csv"),
  row.names = FALSE
)

ec_cluster_patient <- as.data.frame(
  table(
    EC_cluster = ec$EC_cluster,
    patient_id = ec$patient_id
  )
)
write.csv(
  ec_cluster_patient,
  file.path(out_dir, "tables", "07_EC_cluster_by_patient_counts.csv"),
  row.names = FALSE
)

# Positive EC-cluster markers.
Idents(ec) <- "EC_cluster"
ec_markers <- FindAllMarkers(
  ec,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.15,
  logfc.threshold = 0.20,
  max.cells.per.ident = 2500,
  random.seed = 20260729,
  verbose = TRUE
)

write.csv(
  ec_markers,
  file.path(out_dir, "tables", "08_all_positive_EC_cluster_markers.csv"),
  row.names = FALSE
)

if (nrow(ec_markers) > 0) {
  fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(ec_markers))[1]
  if (is.na(fc_col)) {
    stop("Could not identify the fold-change column in EC marker results.")
  }

  split_markers <- split(ec_markers, ec_markers$cluster)
  top20 <- do.call(
    rbind,
    lapply(
      split_markers,
      function(df) {
        ord <- order(df$p_val_adj, -df[[fc_col]], na.last = TRUE)
        head(df[ord, , drop = FALSE], 20)
      }
    )
  )
  rownames(top20) <- NULL

  write.csv(
    top20,
    file.path(out_dir, "tables", "09_top20_markers_each_EC_cluster.csv"),
    row.names = FALSE
  )
}

# AC/PA fractions at patient level.
ec_count_sample_state <- as.data.frame(
  table(
    patient_id = ec$patient_id,
    region = ec$region,
    EC_cluster = ec$EC_cluster
  )
)

sample_totals <- aggregate(
  Freq ~ patient_id + region,
  data = ec_count_sample_state,
  FUN = sum
)
colnames(sample_totals)[3] <- "sample_EC_total"

ec_fraction <- merge(
  ec_count_sample_state,
  sample_totals,
  by = c("patient_id", "region"),
  all.x = TRUE
)
ec_fraction$fraction <- ec_fraction$Freq / ec_fraction$sample_EC_total

write.csv(
  ec_fraction,
  file.path(out_dir, "tables", "10_EC_cluster_fraction_by_patient_region.csv"),
  row.names = FALSE
)

p_frac <- ggplot(
  ec_fraction,
  aes(
    x = interaction(patient_id, region, sep = "_"),
    y = fraction,
    fill = EC_cluster
  )
) +
  geom_col() +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(
    x = NULL,
    y = "Fraction of endothelial cells",
    fill = "EC cluster",
    title = "GSE159677 EC composition by patient and AC/PA region"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(out_dir, "plots", "12_EC_cluster_fraction_by_patient_region.png"),
  p_frac,
  width = 9,
  height = 6,
  dpi = 300
)

saveRDS(
  ec,
  file.path(out_dir, "objects", "GSE159677_high_confidence_EC_reclustered_scored.rds"),
  compress = FALSE
)

cat("\nFinal high-confidence EC cells:", ncol(ec), "\n")
cat("Final local EC clusters:", length(unique(ec$EC_cluster)), "\n")
cat("Formal EC object saved.\n")
cat("IMPORTANT: EC clusters are not yet assigned final biological names.\n")
cat("===== STAGE 5 COMPLETE =====\n")
cat("End:", format(Sys.time()), "\n")
