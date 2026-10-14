#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260813)

BASE <- "/root/GSE218166"
IN_RDS <- file.path(
  BASE,
  "step2A_strict_EC_identity_audit",
  "objects",
  "GSE218166_step2A_ABC_candidate_EC.rds"
)

OUT <- file.path(BASE, "step2B_EC_recluster_lineage_audit")
TAB <- file.path(OUT, "tables")
PLOT <- file.path(OUT, "plots")
OBJ <- file.path(OUT, "objects")
LOG <- file.path(OUT, "logs")

for (d in c(OUT, TAB, PLOT, OBJ, LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

sink(file.path(LOG, "step2B_console.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

cat("============================================================\n")
cat("GSE218166 Step 2B: EC reclustering + lineage-conflict audit\n")
cat("Start:", format(Sys.time()), "\n")
cat("Input:", IN_RDS, "\n")
cat("============================================================\n\n")

if (!file.exists(IN_RDS)) {
  stop("Input RDS not found: ", IN_RDS)
}

obj <- readRDS(IN_RDS)
DefaultAssay(obj) <- "RNA"

cat("Input candidate cells:", ncol(obj), "\n")
cat("Genes:", nrow(obj), "\n")
cat("Tech:", paste(unique(as.character(obj$tech)), collapse = ", "), "\n\n")

# -------------------------------------------------------------------------
# Important design rule:
# disease_stage is carried as metadata for later diagnostics ONLY.
# It is NOT used in any cell identity gate, cluster selection, PCA, Harmony,
# neighbor graph, clustering, or marker-based EC decision.
# -------------------------------------------------------------------------

required_md <- c(
  "donor_id", "coronary_position", "lesion_process",
  "disease_stage", "tech", "seurat_clusters", "plate"
)
missing_md <- setdiff(required_md, colnames(obj@meta.data))
if (length(missing_md) > 0) {
  stop("Missing metadata: ", paste(missing_md, collapse = ", "))
}

obj$disease_stage <- factor(
  as.character(obj$disease_stage),
  levels = c("Normal", "Early", "Advanced"),
  ordered = TRUE
)

# -------------------------------------------------------------------------
# 1. Audit possible physical segment identifiers.
#    We do NOT automatically assume 'plate' is the segment ID.
# -------------------------------------------------------------------------
md <- obj@meta.data
md$unit_DPP <- paste(
  md$donor_id, md$coronary_position, md$lesion_process, md$tech,
  sep = "__"
)
md$unit_DPPP <- paste(
  md$donor_id, md$coronary_position, md$lesion_process, md$plate, md$tech,
  sep = "__"
)

unit_def_summary <- data.frame(
  definition = c(
    "donor+position+process+tech",
    "donor+position+process+plate+tech"
  ),
  n_units = c(
    length(unique(md$unit_DPP)),
    length(unique(md$unit_DPPP))
  )
)

write.csv(
  unit_def_summary,
  file.path(TAB, "01_candidate_sample_unit_definitions.csv"),
  row.names = FALSE
)

unit_plate_audit <- md %>%
  count(
    unit_DPPP,
    donor_id,
    coronary_position,
    lesion_process,
    plate,
    disease_stage,
    tech,
    name = "n_candidate_EC"
  ) %>%
  arrange(tech, donor_id, coronary_position, lesion_process, plate)

write.csv(
  unit_plate_audit,
  file.path(TAB, "02_plate_level_candidate_unit_audit.csv"),
  row.names = FALSE
)

unit_stage_consistency <- unit_plate_audit %>%
  group_by(unit_DPPP) %>%
  summarise(
    n_stage = n_distinct(disease_stage),
    stages = paste(unique(as.character(disease_stage)), collapse = ";"),
    .groups = "drop"
  )

write.csv(
  unit_stage_consistency,
  file.path(TAB, "03_plate_unit_stage_consistency.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 2. Robust RNA layer access
# -------------------------------------------------------------------------
get_layer <- function(x, what = c("data", "counts")) {
  what <- match.arg(what)
  out <- tryCatch(
    GetAssayData(x, assay = "RNA", layer = what),
    error = function(e1) {
      tryCatch(
        GetAssayData(x, assay = "RNA", slot = what),
        error = function(e2) NULL
      )
    }
  )
  out
}

# Normalize the candidate object independently.
obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

rna_data <- get_layer(obj, "data")
rna_counts <- get_layer(obj, "counts")

if (is.null(rna_data)) stop("RNA normalized data layer unavailable.")
if (is.null(rna_counts)) {
  cat("WARNING: counts unavailable; using normalized data > 0 for detection.\n")
  rna_counts <- rna_data
}

# -------------------------------------------------------------------------
# 3. Marker panels.
#    Core EC markers include markers used by the source paper
#    (TIE1/PECAM1) plus canonical vascular EC markers.
# -------------------------------------------------------------------------
marker_sets <- list(
  EC_core = c(
    "TIE1", "PECAM1", "VWF", "KDR", "CDH5", "CLDN5",
    "ESAM", "ENG", "EMCN", "RAMP2", "EPAS1", "ERG", "FLT1"
  ),
  EC_key = c(
    "TIE1", "PECAM1", "KDR", "CDH5", "CLDN5", "ERG", "ESAM", "EMCN"
  ),

  Immune_myeloid = c(
    "PTPRC", "LST1", "TYROBP", "FCER1G", "AIF1", "CTSS", "LYZ", "C5AR1"
  ),
  T_NK = c(
    "CD3D", "CD3E", "CD2", "TRBC2", "NKG7", "GNLY", "CD69"
  ),
  Epithelial = c(
    "EPCAM", "KRT8", "KRT18", "KRT19", "KRT14"
  ),
  Platelet_RBC = c(
    "PPBP", "PF4", "NRGN", "HBB", "HBA1", "HBA2"
  ),

  Mural_contractile = c(
    "ACTA2", "TAGLN", "MYH11", "CNN1", "RGS5", "CSPG4",
    "PDGFRB", "NOTCH3"
  ),
  Fibroblast = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"
  ),

  Arterial = c(
    "CXCL12", "EFNB2", "GJA4", "GJA5", "SOX17", "HEY1", "SEMA3G", "BMX"
  ),
  Venous_microvascular = c(
    "ACKR1", "PLVAP", "COL15A1", "VWF", "CA4", "RGCC"
  ),
  Lymphatic = c(
    "PROX1", "PDPN", "CCL21", "LYVE1", "FLT4"
  ),

  Author_EC5 = c(
    "SLCO4A1", "PRDM15", "PROM1"
  ),
  EndMT_like = c(
    "FN1", "COL1A1", "COL3A1", "COL8A1", "BGN", "MMP2", "VIM"
  ),

  Project_stable_arterial = c(
    "CXCL12", "EFNB2", "GJA4", "SOX17", "HEY1", "KLF2", "KLF4", "THBD", "EPAS1"
  ),
  Project_early_activation = c(
    "ANGPT2", "CCL14", "PGF", "LOX", "CTHRC1", "THBS1"
  ),
  Project_late_remodeling = c(
    "SULF1", "ITLN1", "MGP", "ELN", "COL8A1", "FN1",
    "SPARC", "COL6A2", "IGFBP5"
  )
)

marker_presence <- do.call(
  rbind,
  lapply(names(marker_sets), function(nm) {
    data.frame(
      marker_set = nm,
      gene = marker_sets[[nm]],
      present = marker_sets[[nm]] %in% rownames(obj),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  marker_presence,
  file.path(TAB, "04_marker_panel_presence.csv"),
  row.names = FALSE
)

mean_expr <- function(mat, genes) {
  gg <- intersect(genes, rownames(mat))
  if (length(gg) == 0) return(rep(NA_real_, ncol(mat)))
  Matrix::colMeans(mat[gg, , drop = FALSE])
}

n_detected <- function(mat, genes) {
  gg <- intersect(genes, rownames(mat))
  if (length(gg) == 0) return(rep(0L, ncol(mat)))
  as.integer(Matrix::colSums(mat[gg, , drop = FALSE] > 0))
}

cell_qc <- data.frame(
  cell = colnames(obj),
  ec_core_mean = mean_expr(rna_data, marker_sets$EC_core),
  ec_core_detected = n_detected(rna_counts, marker_sets$EC_core),
  ec_key_detected = n_detected(rna_counts, marker_sets$EC_key),

  immune_mean = mean_expr(rna_data, marker_sets$Immune_myeloid),
  immune_detected = n_detected(rna_counts, marker_sets$Immune_myeloid),

  tnk_mean = mean_expr(rna_data, marker_sets$T_NK),
  tnk_detected = n_detected(rna_counts, marker_sets$T_NK),

  epithelial_mean = mean_expr(rna_data, marker_sets$Epithelial),
  epithelial_detected = n_detected(rna_counts, marker_sets$Epithelial),

  platelet_rbc_mean = mean_expr(rna_data, marker_sets$Platelet_RBC),
  platelet_rbc_detected = n_detected(rna_counts, marker_sets$Platelet_RBC),

  mural_mean = mean_expr(rna_data, marker_sets$Mural_contractile),
  mural_detected = n_detected(rna_counts, marker_sets$Mural_contractile),

  fibro_mean = mean_expr(rna_data, marker_sets$Fibroblast),
  fibro_detected = n_detected(rna_counts, marker_sets$Fibroblast),

  arterial_mean = mean_expr(rna_data, marker_sets$Arterial),
  venous_microvascular_mean = mean_expr(rna_data, marker_sets$Venous_microvascular),
  lymphatic_mean = mean_expr(rna_data, marker_sets$Lymphatic),
  author_EC5_mean = mean_expr(rna_data, marker_sets$Author_EC5),
  endmt_like_mean = mean_expr(rna_data, marker_sets$EndMT_like),

  project_stable_mean = mean_expr(rna_data, marker_sets$Project_stable_arterial),
  project_early_mean = mean_expr(rna_data, marker_sets$Project_early_activation),
  project_late_mean = mean_expr(rna_data, marker_sets$Project_late_remodeling),

  stringsAsFactors = FALSE
)

meta_keep <- c(
  "donor_id", "coronary_position", "lesion_process",
  "disease_stage", "tech", "plate",
  "seurat_clusters", "EC_review_tier"
)

cell_qc <- cbind(
  cell_qc,
  obj@meta.data[cell_qc$cell, meta_keep, drop = FALSE]
)

# Hard non-EC conflicts require MULTIPLE detected markers.
# Single-marker positivity is intentionally not enough because ambient RNA
# is expected in vascular plaque preparations.
cell_qc <- cell_qc %>%
  mutate(
    immune_conflict = immune_detected >= 3,
    tnk_conflict = tnk_detected >= 3,
    epithelial_conflict = epithelial_detected >= 3,
    platelet_rbc_conflict = platelet_rbc_detected >= 3,

    hard_nonEC_conflict =
      immune_conflict |
      tnk_conflict |
      epithelial_conflict |
      platelet_rbc_conflict,

    canonical_EC_positive =
      ec_core_detected >= 5 &
      ec_key_detected >= 2,

    # Mural/fibro signals are NOT hard exclusions, because genuine
    # EndMT-like or remodeling EC may co-express them.
    mural_high = mural_detected >= 5,
    fibro_high = fibro_detected >= 4,

    strict_cell_identity_support =
      canonical_EC_positive &
      !hard_nonEC_conflict
  )

write.csv(
  cell_qc,
  gzfile(file.path(TAB, "05_cell_level_lineage_audit.csv.gz")),
  row.names = FALSE
)

# Add fields to Seurat object
for (cc in setdiff(colnames(cell_qc), c("cell", meta_keep))) {
  obj[[cc]] <- cell_qc[[cc]][match(colnames(obj), cell_qc$cell)]
}

# -------------------------------------------------------------------------
# 4. Stage-blind EC-only reclustering
# -------------------------------------------------------------------------
obj <- FindVariableFeatures(
  obj,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 2500,
  verbose = FALSE
)

obj <- ScaleData(
  obj,
  assay = "RNA",
  features = VariableFeatures(obj),
  verbose = FALSE
)

npcs_use <- min(30, length(VariableFeatures(obj)) - 1)
obj <- RunPCA(
  obj,
  assay = "RNA",
  features = VariableFeatures(obj),
  npcs = npcs_use,
  reduction.name = "pca.ec",
  verbose = FALSE
)

reduction_use <- "pca.ec"

# Harmony is used only to control the known FL vs 3'DE technology split.
# disease_stage is NEVER passed to Harmony.
if (requireNamespace("harmony", quietly = TRUE)) {
  cat("Running Harmony by tech only...\n")
  obj <- harmony::RunHarmony(
    object = obj,
    group.by.vars = "tech",
    reduction.use = "pca.ec",
    reduction.save = "harmony.ec",
    assay.use = "RNA",
    verbose = FALSE
  )
  reduction_use <- "harmony.ec"
} else {
  cat("WARNING: harmony package not available; using PCA for EC reclustering.\n")
}

dims_use <- 1:min(20, ncol(Embeddings(obj, reduction = reduction_use)))

obj <- FindNeighbors(
  obj,
  reduction = reduction_use,
  dims = dims_use,
  graph.name = c("EC_nn", "EC_snn"),
  verbose = FALSE
)

obj <- FindClusters(
  obj,
  graph.name = "EC_snn",
  resolution = 0.45,
  random.seed = 20260813,
  verbose = FALSE
)

obj$EC_recluster <- as.character(Idents(obj))

obj <- RunUMAP(
  obj,
  reduction = reduction_use,
  dims = dims_use,
  reduction.name = "umap.ec",
  seed.use = 20260813,
  verbose = FALSE
)

# -------------------------------------------------------------------------
# 5. Cluster-level lineage audit
# -------------------------------------------------------------------------
cluster_qc <- obj@meta.data %>%
  group_by(EC_recluster) %>%
  summarise(
    n_cells = n(),
    n_donors = n_distinct(donor_id),
    n_tech = n_distinct(tech),

    median_ec_core = median(ec_core_mean, na.rm = TRUE),
    median_ec_core_detected = median(ec_core_detected, na.rm = TRUE),
    median_ec_key_detected = median(ec_key_detected, na.rm = TRUE),

    frac_canonical_EC_positive = mean(canonical_EC_positive, na.rm = TRUE),
    frac_strict_identity_support = mean(strict_cell_identity_support, na.rm = TRUE),
    frac_hard_nonEC_conflict = mean(hard_nonEC_conflict, na.rm = TRUE),

    frac_mural_high = mean(mural_high, na.rm = TRUE),
    frac_fibro_high = mean(fibro_high, na.rm = TRUE),

    median_arterial = median(arterial_mean, na.rm = TRUE),
    median_venous_microvascular = median(venous_microvascular_mean, na.rm = TRUE),
    median_lymphatic = median(lymphatic_mean, na.rm = TRUE),
    median_author_EC5 = median(author_EC5_mean, na.rm = TRUE),
    median_EndMT_like = median(endmt_like_mean, na.rm = TRUE),

    median_project_stable = median(project_stable_mean, na.rm = TRUE),
    median_project_early = median(project_early_mean, na.rm = TRUE),
    median_project_late = median(project_late_mean, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  mutate(
    auto_review_flag = case_when(
      frac_canonical_EC_positive >= 0.70 &
        frac_hard_nonEC_conflict <= 0.15 ~ "HIGH_CONFIDENCE_EC_REVIEW",

      frac_canonical_EC_positive >= 0.50 &
        frac_hard_nonEC_conflict <= 0.25 ~ "POSSIBLE_EC_REVIEW",

      TRUE ~ "MIXED_OR_NON_EC_REVIEW"
    )
  ) %>%
  arrange(
    factor(
      auto_review_flag,
      levels = c(
        "HIGH_CONFIDENCE_EC_REVIEW",
        "POSSIBLE_EC_REVIEW",
        "MIXED_OR_NON_EC_REVIEW"
      )
    ),
    desc(frac_canonical_EC_positive),
    frac_hard_nonEC_conflict
  )

write.csv(
  cluster_qc,
  file.path(TAB, "06_EC_recluster_lineage_QC_summary.csv"),
  row.names = FALSE
)

flag_map <- setNames(cluster_qc$auto_review_flag, cluster_qc$EC_recluster)
obj$EC_cluster_review_flag <- unname(flag_map[as.character(obj$EC_recluster)])

# -------------------------------------------------------------------------
# 6. Cluster markers
# -------------------------------------------------------------------------
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "EC_recluster"

# Seurat v5 may require joined layers for marker testing.
if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  obj[["RNA"]] <- tryCatch(
    SeuratObject::JoinLayers(obj[["RNA"]]),
    error = function(e) obj[["RNA"]]
  )
}

markers <- tryCatch(
  FindAllMarkers(
    obj,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = 0.20,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    verbose = FALSE
  ),
  error = function(e) {
    cat("WARNING FindAllMarkers failed:", conditionMessage(e), "\n")
    data.frame()
  }
)

if (nrow(markers) > 0) {
  write.csv(
    markers,
    file.path(TAB, "07_EC_recluster_all_positive_markers.csv"),
    row.names = FALSE
  )

  fc_col <- intersect(
    c("avg_log2FC", "avg_logFC"),
    colnames(markers)
  )[1]

  if (!is.na(fc_col) && length(fc_col) == 1) {
    top_markers <- markers %>%
      group_by(cluster) %>%
      arrange(desc(.data[[fc_col]]), .by_group = TRUE) %>%
      slice_head(n = 20) %>%
      ungroup()

    write.csv(
      top_markers,
      file.path(TAB, "08_EC_recluster_top20_markers.csv"),
      row.names = FALSE
    )
  }
}

# -------------------------------------------------------------------------
# 7. Stage distribution is diagnostic only and cannot affect identity freeze
# -------------------------------------------------------------------------
stage_cluster <- obj@meta.data %>%
  count(tech, disease_stage, EC_recluster, name = "n_cells") %>%
  group_by(tech, disease_stage) %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(
  stage_cluster,
  file.path(TAB, "09_stage_distribution_by_EC_recluster_DIAGNOSTIC_ONLY.csv"),
  row.names = FALSE
)

strict_support_counts <- obj@meta.data %>%
  count(
    tech,
    disease_stage,
    strict_cell_identity_support,
    name = "n_cells"
  ) %>%
  group_by(tech, disease_stage) %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(
  strict_support_counts,
  file.path(TAB, "10_strict_identity_support_by_stage_tech_DIAGNOSTIC_ONLY.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 8. Figures - always pure white background
# -------------------------------------------------------------------------
save_white <- function(p, filename, width = 10, height = 7) {
  ggsave(
    file.path(PLOT, filename),
    p,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

p1 <- DimPlot(
  obj,
  reduction = "umap.ec",
  group.by = "EC_recluster",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 0.35
) +
  ggtitle("GSE218166 candidate EC reclustering") +
  theme_classic(base_size = 11)

save_white(p1, "01_EC_recluster_UMAP.png", 10, 8)

p2 <- DimPlot(
  obj,
  reduction = "umap.ec",
  group.by = "EC_cluster_review_flag",
  raster = TRUE,
  pt.size = 0.35
) +
  ggtitle("Cluster-level lineage audit flag") +
  theme_classic(base_size = 11)

save_white(p2, "02_EC_cluster_review_flag_UMAP.png", 10, 8)

feature_cols <- c(
  "ec_core_mean", "immune_mean", "mural_mean", "fibro_mean",
  "author_EC5_mean", "endmt_like_mean"
)

p3 <- FeaturePlot(
  obj,
  reduction = "umap.ec",
  features = feature_cols,
  ncol = 3,
  raster = TRUE
) & theme_classic(base_size = 9)

save_white(p3, "03_lineage_and_state_scores_UMAP.png", 15, 10)

dot_genes <- unique(c(
  "TIE1", "PECAM1", "VWF", "KDR", "CDH5", "CLDN5", "ERG", "ESAM",
  "PTPRC", "LST1", "TYROBP", "C5AR1",
  "CD3D", "CD3E", "NKG7",
  "ACTA2", "TAGLN", "MYH11", "CNN1", "RGS5", "PDGFRB",
  "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM",
  "PPBP", "PF4", "HBB",
  "ACKR1", "PLVAP", "COL15A1",
  "CXCL12", "EFNB2", "GJA4", "SOX17", "HEY1",
  "PROX1", "PDPN", "CCL21",
  "SLCO4A1", "PRDM15", "PROM1",
  "ANGPT2", "PGF", "LOX", "CTHRC1", "THBS1",
  "SULF1", "ITLN1", "MGP", "ELN", "COL8A1", "FN1"
))
dot_genes <- intersect(dot_genes, rownames(obj))

if (length(dot_genes) > 0) {
  p4 <- DotPlot(
    obj,
    features = dot_genes,
    group.by = "EC_recluster",
    assay = "RNA"
  ) +
    RotatedAxis() +
    theme_classic(base_size = 8) +
    ggtitle("EC recluster marker audit")

  save_white(p4, "04_EC_recluster_marker_DotPlot.png", 19, 9)
}

p5 <- ggplot(
  stage_cluster,
  aes(x = disease_stage, y = fraction, fill = EC_recluster)
) +
  geom_col() +
  facet_wrap(~tech) +
  theme_classic(base_size = 11) +
  labs(
    title = "Stage composition after stage-blind EC reclustering",
    subtitle = "Diagnostic only; disease stage was not used to define EC identity",
    x = NULL,
    y = "Fraction of candidate EC",
    fill = "EC recluster"
  )

save_white(p5, "05_stage_composition_DIAGNOSTIC_ONLY.png", 11, 6.5)

# Heatmap-style table for cluster QC
cluster_qc_long <- cluster_qc %>%
  select(
    EC_recluster,
    frac_canonical_EC_positive,
    frac_strict_identity_support,
    frac_hard_nonEC_conflict,
    frac_mural_high,
    frac_fibro_high
  ) %>%
  pivot_longer(
    -EC_recluster,
    names_to = "metric",
    values_to = "value"
  )

p6 <- ggplot(
  cluster_qc_long,
  aes(x = metric, y = factor(EC_recluster), fill = value)
) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(title = "EC recluster lineage QC fractions")

save_white(p6, "06_EC_recluster_QC_heatmap.png", 11, 8)

# -------------------------------------------------------------------------
# 9. Save full review object (NOT final frozen EC yet)
# -------------------------------------------------------------------------
saveRDS(
  obj,
  file.path(OBJ, "GSE218166_step2B_candidate_EC_reclustered_review.rds"),
  compress = FALSE
)

# Cell barcode sets for later manual final freeze
high_clusters <- cluster_qc$EC_recluster[
  cluster_qc$auto_review_flag == "HIGH_CONFIDENCE_EC_REVIEW"
]
possible_clusters <- cluster_qc$EC_recluster[
  cluster_qc$auto_review_flag %in%
    c("HIGH_CONFIDENCE_EC_REVIEW", "POSSIBLE_EC_REVIEW")
]

writeLines(
  colnames(obj)[obj$EC_recluster %in% high_clusters],
  file.path(TAB, "11_auto_high_confidence_cluster_cells.txt")
)

writeLines(
  colnames(obj)[obj$EC_recluster %in% possible_clusters],
  file.path(TAB, "12_auto_high_plus_possible_cluster_cells.txt")
)

# -------------------------------------------------------------------------
# 10. Decision summary
# -------------------------------------------------------------------------
n_input <- ncol(obj)
n_strict_support <- sum(obj$strict_cell_identity_support, na.rm = TRUE)

summary_lines <- c(
  "============================================================",
  "GSE218166 Step 2B EC recluster + lineage audit",
  "============================================================",
  paste0("Run time: ", format(Sys.time())),
  paste0("Input A+B+C candidate EC: ", n_input),
  paste0("Cell-level strict identity support: ", n_strict_support),
  paste0("Number of EC reclusters: ", length(unique(obj$EC_recluster))),
  paste0("Reduction used for clustering: ", reduction_use),
  "",
  "CRITICAL INTERPRETATION RULES",
  "1. Disease stage was NOT used for EC reclustering or identity decisions.",
  "2. A/B/C from Step2A are treated only as an inclusive candidate pool.",
  "3. Mural/fibro expression is NOT an automatic exclusion because remodeling/EndMT-like EC can co-express these programs.",
  "4. Immune/T-NK/epithelial/platelet-RBC conflicts require multiple detected markers before being treated as hard lineage conflicts.",
  "5. The automatic cluster flags are REVIEW FLAGS, not the final biological freeze.",
  "",
  "SAMPLE UNIT AUDIT",
  paste0(
    "donor+position+process+tech units = ",
    unit_def_summary$n_units[unit_def_summary$definition == "donor+position+process+tech"]
  ),
  paste0(
    "donor+position+process+plate+tech units = ",
    unit_def_summary$n_units[unit_def_summary$definition == "donor+position+process+plate+tech"]
  ),
  "Review whether 'plate' resolves the published physical coronary-segment structure before pseudobulk.",
  "",
  "NEXT",
  "Upload the Step2B review bundle.",
  "Final EC clusters will then be frozen in Step2C before any Normal -> Early -> Advanced testing."
)

writeLines(
  summary_lines,
  file.path(OUT, "Step2B_DECISION_SUMMARY.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(LOG, "sessionInfo.txt")
)

cat(paste(summary_lines, collapse = "\n"))
cat("\n\nCompleted:", format(Sys.time()), "\n")
