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
  "GSE218166_step2A_allcells_with_EC_identity_scores.rds"
)

OUT  <- file.path(BASE, "step2B_V3_EC_recluster_lineage_audit")
TAB  <- file.path(OUT, "tables")
PLOT <- file.path(OUT, "plots")
OBJ  <- file.path(OUT, "objects")
LOG  <- file.path(OUT, "logs")

for (d in c(OUT, TAB, PLOT, OBJ, LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

sink(file.path(LOG, "step2B_V3_console.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

cat("============================================================\n")
cat("GSE218166 STEP 2B V3 - STANDALONE HARMONY\n")
cat("Absolute-marker EC rescue + stage-blind EC reclustering\n")
cat("Start:", format(Sys.time()), "\n")
cat("Input:", IN_RDS, "\n")
cat("============================================================\n\n")

if (!file.exists(IN_RDS)) {
  stop("Input RDS not found: ", IN_RDS)
}

obj_all <- readRDS(IN_RDS)
DefaultAssay(obj_all) <- "RNA"

cat("ALL Step2A cells:", ncol(obj_all), "\n")
cat("Genes:", nrow(obj_all), "\n")
cat("Tech:", paste(unique(as.character(obj_all$tech)), collapse = ", "), "\n\n")

required_md <- c(
  "donor_id", "coronary_position", "lesion_process",
  "disease_stage", "tech", "plate", "seurat_clusters",
  "EC_review_tier"
)
missing_md <- setdiff(required_md, colnames(obj_all@meta.data))
if (length(missing_md) > 0) {
  stop("Missing metadata: ", paste(missing_md, collapse = ", "))
}

obj_all$disease_stage <- factor(
  as.character(obj_all$disease_stage),
  levels = c("Normal", "Early", "Advanced"),
  ordered = TRUE
)

# -------------------------------------------------------------------------
# 1. Replicate-unit audit
# -------------------------------------------------------------------------
obj_all$Individual_Level_reconstructed <- paste(
  obj_all$donor_id,
  as.character(obj_all$disease_stage),
  sep = "__"
)

obj_all$segment_process_unit <- paste(
  obj_all$donor_id,
  obj_all$coronary_position,
  obj_all$lesion_process,
  obj_all$tech,
  sep = "__"
)

obj_all$segment_process_plate_unit <- paste(
  obj_all$donor_id,
  obj_all$coronary_position,
  obj_all$lesion_process,
  obj_all$plate,
  obj_all$tech,
  sep = "__"
)

unit_summary <- data.frame(
  definition = c(
    "donor x disease_stage (Individual_Level_reconstructed)",
    "donor x position x lesion_process x tech",
    "donor x position x lesion_process x plate x tech"
  ),
  n_units = c(
    length(unique(obj_all$Individual_Level_reconstructed)),
    length(unique(obj_all$segment_process_unit)),
    length(unique(obj_all$segment_process_plate_unit))
  )
)

write.csv(
  unit_summary,
  file.path(TAB, "01_statistical_unit_audit_summary.csv"),
  row.names = FALSE
)

unit_detail <- obj_all@meta.data %>%
  count(
    Individual_Level_reconstructed,
    donor_id,
    disease_stage,
    tech,
    plate,
    coronary_position,
    lesion_process,
    name = "n_cells"
  ) %>%
  arrange(tech, donor_id, disease_stage, plate, coronary_position)

write.csv(
  unit_detail,
  file.path(TAB, "02_Individual_Level_plate_segment_audit.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 2. Absolute EC identity evidence across ALL Step2A cells
#    disease_stage is NOT used anywhere in identity selection.
# -------------------------------------------------------------------------
obj_all <- NormalizeData(
  obj_all,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

get_layer <- function(x, what = c("data", "counts")) {
  what <- match.arg(what)
  ans <- tryCatch(
    GetAssayData(x, assay = "RNA", layer = what),
    error = function(e1) {
      tryCatch(
        GetAssayData(x, assay = "RNA", slot = what),
        error = function(e2) NULL
      )
    }
  )
  ans
}

rna_data <- get_layer(obj_all, "data")
rna_counts <- get_layer(obj_all, "counts")

if (is.null(rna_data)) stop("RNA normalized data not available.")
if (is.null(rna_counts)) {
  cat("WARNING: counts layer not available; using normalized data > 0.\n")
  rna_counts <- rna_data
}

marker_sets <- list(
  EC_core = c(
    "TIE1","PECAM1","VWF","KDR","CDH5","CLDN5","ESAM",
    "ENG","EMCN","RAMP2","EPAS1","ERG","FLT1"
  ),
  EC_key = c(
    "TIE1","PECAM1","KDR","CDH5","CLDN5","ERG","ESAM","EMCN"
  ),
  Immune = c(
    "PTPRC","LST1","TYROBP","FCER1G","AIF1","CTSS","LYZ","C5AR1"
  ),
  T_NK = c(
    "CD3D","CD3E","CD2","TRBC1","TRBC2","NKG7","GNLY","CD69"
  ),
  Epithelial = c(
    "EPCAM","KRT8","KRT18","KRT19","KRT14"
  ),
  Platelet_RBC = c(
    "PPBP","PF4","NRGN","HBB","HBA1","HBA2"
  ),
  Mural = c(
    "ACTA2","TAGLN","MYH11","CNN1","RGS5","CSPG4","PDGFRB","NOTCH3"
  ),
  Fibroblast = c(
    "COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA"
  ),
  Arterial = c(
    "CXCL12","EFNB2","GJA4","GJA5","SOX17","HEY1","SEMA3G","BMX"
  ),
  Venous_microvascular = c(
    "ACKR1","PLVAP","COL15A1","CA4","RGCC"
  ),
  Lymphatic = c(
    "PROX1","PDPN","CCL21","LYVE1","FLT4"
  ),
  Author_EC5 = c(
    "SLCO4A1","PRDM15","PROM1"
  ),
  Project_stable = c(
    "CXCL12","EFNB2","GJA4","SOX17","HEY1","KLF2","KLF4","THBD","EPAS1"
  ),
  Project_early = c(
    "ANGPT2","CCL14","PGF","LOX","CTHRC1","THBS1"
  ),
  Project_late = c(
    "SULF1","ITLN1","MGP","ELN","COL8A1","FN1","SPARC","COL6A2","IGFBP5"
  ),
  EndMT_like = c(
    "FN1","COL1A1","COL3A1","COL8A1","BGN","MMP2","VIM"
  )
)

presence <- do.call(
  rbind,
  lapply(names(marker_sets), function(nm) {
    data.frame(
      marker_set = nm,
      gene = marker_sets[[nm]],
      present = marker_sets[[nm]] %in% rownames(rna_data),
      stringsAsFactors = FALSE
    )
  })
)
write.csv(
  presence,
  file.path(TAB, "03_marker_panel_presence.csv"),
  row.names = FALSE
)

n_detected <- function(mat, genes) {
  gg <- intersect(genes, rownames(mat))
  if (length(gg) == 0) return(rep(0L, ncol(mat)))
  as.integer(Matrix::colSums(mat[gg, , drop = FALSE] > 0))
}

mean_expr <- function(mat, genes) {
  gg <- intersect(genes, rownames(mat))
  if (length(gg) == 0) return(rep(NA_real_, ncol(mat)))
  Matrix::colMeans(mat[gg, , drop = FALSE])
}

audit <- data.frame(
  cell = colnames(obj_all),
  ec_core_detected = n_detected(rna_counts, marker_sets$EC_core),
  ec_key_detected = n_detected(rna_counts, marker_sets$EC_key),
  immune_detected = n_detected(rna_counts, marker_sets$Immune),
  tnk_detected = n_detected(rna_counts, marker_sets$T_NK),
  epithelial_detected = n_detected(rna_counts, marker_sets$Epithelial),
  platelet_rbc_detected = n_detected(rna_counts, marker_sets$Platelet_RBC),
  mural_detected = n_detected(rna_counts, marker_sets$Mural),
  fibro_detected = n_detected(rna_counts, marker_sets$Fibroblast),
  author_EC5_detected = n_detected(rna_counts, marker_sets$Author_EC5),

  ec_core_mean = mean_expr(rna_data, marker_sets$EC_core),
  arterial_mean = mean_expr(rna_data, marker_sets$Arterial),
  venous_microvascular_mean = mean_expr(rna_data, marker_sets$Venous_microvascular),
  lymphatic_mean = mean_expr(rna_data, marker_sets$Lymphatic),
  author_EC5_mean = mean_expr(rna_data, marker_sets$Author_EC5),
  project_stable_mean = mean_expr(rna_data, marker_sets$Project_stable),
  project_early_mean = mean_expr(rna_data, marker_sets$Project_early),
  project_late_mean = mean_expr(rna_data, marker_sets$Project_late),
  endmt_like_mean = mean_expr(rna_data, marker_sets$EndMT_like),
  stringsAsFactors = FALSE
)

meta_cols <- c(
  "donor_id","coronary_position","lesion_process","disease_stage",
  "tech","plate","seurat_clusters","EC_review_tier",
  "Individual_Level_reconstructed"
)

audit <- cbind(
  audit,
  obj_all@meta.data[audit$cell, meta_cols, drop = FALSE]
)

audit <- audit %>%
  mutate(
    step2A_candidate = as.character(EC_review_tier) != "Not_candidate",

    immune_conflict = immune_detected >= 4,
    tnk_conflict = tnk_detected >= 4,
    epithelial_conflict = epithelial_detected >= 3,
    platelet_rbc_conflict = platelet_rbc_detected >= 3,

    hard_nonEC_conflict =
      immune_conflict |
      tnk_conflict |
      epithelial_conflict |
      platelet_rbc_conflict,

    canonical_EC_support =
      ec_core_detected >= 5 &
      ec_key_detected >= 2,

    inclusive_EC_support =
      ec_core_detected >= 3 &
      ec_key_detected >= 1,

    atypical_EC5_rescue =
      ec_core_detected >= 2 &
      author_EC5_detected >= 2,

    # Inclusive stage-blind rescue.
    # Mural/fibroblast programs are deliberately NOT hard exclusions.
    prefilter_keep =
      (
        step2A_candidate |
        inclusive_EC_support |
        atypical_EC5_rescue
      ) &
      !hard_nonEC_conflict,

    prefilter_origin = case_when(
      prefilter_keep & step2A_candidate & canonical_EC_support
        ~ "Step2A_plus_canonical",
      prefilter_keep & step2A_candidate
        ~ "Step2A_only_or_weak_absolute",
      prefilter_keep & !step2A_candidate & canonical_EC_support
        ~ "RESCUED_canonical_absolute",
      prefilter_keep & !step2A_candidate & atypical_EC5_rescue
        ~ "RESCUED_EC5_like",
      prefilter_keep & !step2A_candidate
        ~ "RESCUED_inclusive_absolute",
      TRUE ~ "Excluded"
    )
  )

write.csv(
  audit,
  gzfile(file.path(TAB, "04_allcell_absolute_EC_lineage_audit.csv.gz")),
  row.names = FALSE
)

origin_summary <- audit %>%
  count(prefilter_origin, name = "n_cells") %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  arrange(desc(n_cells))

write.csv(
  origin_summary,
  file.path(TAB, "05_prefilter_origin_summary.csv"),
  row.names = FALSE
)

rescue_cluster <- audit %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_all = n(),
    n_step2A_candidate = sum(step2A_candidate),
    n_keep = sum(prefilter_keep),
    n_rescued_from_NotCandidate = sum(prefilter_keep & !step2A_candidate),
    frac_keep = mean(prefilter_keep),
    median_ec_core_detected = median(ec_core_detected),
    median_ec_key_detected = median(ec_key_detected),
    frac_hard_nonEC_conflict = mean(hard_nonEC_conflict),
    .groups = "drop"
  ) %>%
  arrange(desc(n_rescued_from_NotCandidate), desc(frac_keep))

write.csv(
  rescue_cluster,
  file.path(TAB, "06_prefilter_rescue_by_broad_cluster.csv"),
  row.names = FALSE
)

# Transfer audit variables into all-cell object.
audit_add <- c(
  "ec_core_detected","ec_key_detected",
  "immune_detected","tnk_detected","epithelial_detected",
  "platelet_rbc_detected","mural_detected","fibro_detected",
  "author_EC5_detected","ec_core_mean","arterial_mean",
  "venous_microvascular_mean","lymphatic_mean","author_EC5_mean",
  "project_stable_mean","project_early_mean","project_late_mean",
  "endmt_like_mean","step2A_candidate","hard_nonEC_conflict",
  "canonical_EC_support","inclusive_EC_support",
  "atypical_EC5_rescue","prefilter_keep","prefilter_origin"
)

for (cc in audit_add) {
  obj_all[[cc]] <- audit[[cc]][match(colnames(obj_all), audit$cell)]
}

keep_cells <- audit$cell[audit$prefilter_keep]

cat("Stage-blind EC candidates retained:", length(keep_cells), "\n")
cat(
  "Rescued from Step2A Not_candidate:",
  sum(audit$prefilter_keep & !audit$step2A_candidate),
  "\n"
)

if (length(keep_cells) < 500) {
  stop("Too few EC candidates after rescue: ", length(keep_cells))
}

obj <- subset(obj_all, cells = keep_cells)
DefaultAssay(obj) <- "RNA"

# -------------------------------------------------------------------------
# 3. EC-only reclustering
# -------------------------------------------------------------------------
obj <- NormalizeData(obj, verbose = FALSE)

obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 2500,
  verbose = FALSE
)

obj <- ScaleData(
  obj,
  features = VariableFeatures(obj),
  verbose = FALSE
)

obj <- RunPCA(
  obj,
  features = VariableFeatures(obj),
  npcs = 30,
  reduction.name = "pca.ec",
  reduction.key = "PCAEC_",
  verbose = FALSE
)

reduction_use <- "pca.ec"

# -------------------------------------------------------------------------
# 4. Harmony: STANDALONE matrix API
#    This intentionally bypasses the Seurat RunHarmony wrapper.
#    Therefore no assay.use / reduction.use legacy argument can be forwarded.
#    Only tech is corrected. disease_stage is NOT corrected.
# -------------------------------------------------------------------------
if (requireNamespace("harmony", quietly = TRUE)) {
  cat("\nRunning Harmony STANDALONE; correction variable = tech ONLY\n")
  cat("Harmony version:", as.character(utils::packageVersion("harmony")), "\n")

  pca_embed <- Seurat::Embeddings(obj, reduction = "pca.ec")
  harmony_meta <- obj@meta.data[rownames(pca_embed), , drop = FALSE]

  if (!identical(rownames(pca_embed), rownames(harmony_meta))) {
    stop("PCA/meta cell order mismatch before Harmony.")
  }

  harmony_embed <- harmony::RunHarmony(
    data_mat = pca_embed,
    meta_data = harmony_meta,
    vars_use = "tech",
    return_object = FALSE,
    verbose = FALSE
  )

  harmony_embed <- as.matrix(harmony_embed)
  rownames(harmony_embed) <- rownames(pca_embed)
  colnames(harmony_embed) <- paste0(
    "HARMONYEC_",
    seq_len(ncol(harmony_embed))
  )

  obj[["harmony.ec"]] <- Seurat::CreateDimReducObject(
    embeddings = harmony_embed,
    key = "HARMONYEC_",
    assay = DefaultAssay(obj)
  )

  reduction_use <- "harmony.ec"

  cat(
    "Harmony standalone completed:",
    nrow(harmony_embed), "cells x",
    ncol(harmony_embed), "dimensions\n"
  )
} else {
  cat("WARNING: harmony package absent; PCA will be used.\n")
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
  reduction.key = "UMAPEC_",
  seed.use = 20260813,
  verbose = FALSE
)

# -------------------------------------------------------------------------
# 5. Cluster-level EC lineage QC
# -------------------------------------------------------------------------
cluster_qc <- obj@meta.data %>%
  group_by(EC_recluster) %>%
  summarise(
    n_cells = n(),
    n_donors = n_distinct(donor_id),
    n_Individual_Level = n_distinct(Individual_Level_reconstructed),
    n_tech = n_distinct(tech),

    median_ec_core_detected = median(ec_core_detected, na.rm = TRUE),
    median_ec_key_detected = median(ec_key_detected, na.rm = TRUE),

    frac_canonical_EC_support = mean(canonical_EC_support, na.rm = TRUE),
    frac_inclusive_EC_support = mean(inclusive_EC_support, na.rm = TRUE),
    frac_hard_nonEC_conflict = mean(hard_nonEC_conflict, na.rm = TRUE),

    frac_mural_4plus = mean(mural_detected >= 4, na.rm = TRUE),
    frac_fibro_3plus = mean(fibro_detected >= 3, na.rm = TRUE),

    median_arterial = median(arterial_mean, na.rm = TRUE),
    median_venous_microvascular = median(venous_microvascular_mean, na.rm = TRUE),
    median_lymphatic = median(lymphatic_mean, na.rm = TRUE),
    median_author_EC5 = median(author_EC5_mean, na.rm = TRUE),

    median_project_stable = median(project_stable_mean, na.rm = TRUE),
    median_project_early = median(project_early_mean, na.rm = TRUE),
    median_project_late = median(project_late_mean, na.rm = TRUE),
    median_EndMT_like = median(endmt_like_mean, na.rm = TRUE),

    frac_rescued_from_Step2A_NotCandidate =
      mean(!step2A_candidate, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  mutate(
    auto_review_flag = case_when(
      frac_canonical_EC_support >= 0.60 &
        frac_hard_nonEC_conflict <= 0.10
        ~ "HIGH_CONFIDENCE_EC_REVIEW",

      frac_inclusive_EC_support >= 0.60 &
        frac_hard_nonEC_conflict <= 0.15
        ~ "POSSIBLE_EC_REVIEW",

      TRUE
        ~ "MIXED_OR_NON_EC_REVIEW"
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
    desc(frac_canonical_EC_support)
  )

write.csv(
  cluster_qc,
  file.path(TAB, "07_EC_recluster_lineage_QC_summary.csv"),
  row.names = FALSE
)

flag_map <- setNames(
  cluster_qc$auto_review_flag,
  cluster_qc$EC_recluster
)
obj$EC_cluster_review_flag <- unname(
  flag_map[as.character(obj$EC_recluster)]
)

# -------------------------------------------------------------------------
# 6. EC cluster markers
# -------------------------------------------------------------------------
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "EC_recluster"

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
    file.path(TAB, "08_EC_recluster_all_positive_markers.csv"),
    row.names = FALSE
  )

  fc_col <- intersect(
    c("avg_log2FC", "avg_logFC"),
    colnames(markers)
  )[1]

  if (length(fc_col) == 1 && !is.na(fc_col)) {
    top_markers <- markers %>%
      group_by(cluster) %>%
      arrange(desc(.data[[fc_col]]), .by_group = TRUE) %>%
      slice_head(n = 20) %>%
      ungroup()

    write.csv(
      top_markers,
      file.path(TAB, "09_EC_recluster_top20_markers.csv"),
      row.names = FALSE
    )
  }
}

# -------------------------------------------------------------------------
# 7. Stage diagnostics AFTER stage-blind identity/reclustering
# -------------------------------------------------------------------------
stage_cluster <- obj@meta.data %>%
  count(tech, disease_stage, EC_recluster, name = "n_cells") %>%
  group_by(tech, disease_stage) %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(
  stage_cluster,
  file.path(
    TAB,
    "10_stage_distribution_by_EC_recluster_DIAGNOSTIC_ONLY.csv"
  ),
  row.names = FALSE
)

replicate_counts <- obj@meta.data %>%
  count(
    Individual_Level_reconstructed,
    donor_id,
    disease_stage,
    tech,
    name = "n_candidate_EC"
  ) %>%
  arrange(tech, disease_stage, donor_id)

write.csv(
  replicate_counts,
  file.path(TAB, "11_candidate_EC_counts_by_Individual_Level.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 8. Figures - PURE WHITE background
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
  pt.size = 0.30
) +
  ggtitle("GSE218166 candidate EC reclustering") +
  theme_classic(base_size = 11)

save_white(p1, "01_EC_recluster_UMAP.png", 10, 8)

p2 <- DimPlot(
  obj,
  reduction = "umap.ec",
  group.by = "prefilter_origin",
  raster = TRUE,
  pt.size = 0.30
) +
  ggtitle("Step2A candidates and absolute-marker rescued cells") +
  theme_classic(base_size = 11)

save_white(p2, "02_prefilter_origin_UMAP.png", 11, 8)

p3 <- DimPlot(
  obj,
  reduction = "umap.ec",
  group.by = "EC_cluster_review_flag",
  raster = TRUE,
  pt.size = 0.30
) +
  ggtitle("Cluster-level EC lineage review flag") +
  theme_classic(base_size = 11)

save_white(p3, "03_EC_cluster_review_flag_UMAP.png", 10, 8)

p4 <- FeaturePlot(
  obj,
  reduction = "umap.ec",
  features = c(
    "ec_core_detected",
    "mural_detected",
    "fibro_detected",
    "author_EC5_mean",
    "project_stable_mean",
    "project_late_mean"
  ),
  ncol = 3,
  raster = TRUE
) & theme_classic(base_size = 9)

save_white(p4, "04_identity_and_project_scores_UMAP.png", 15, 10)

dot_genes <- unique(c(
  "TIE1","PECAM1","VWF","KDR","CDH5","CLDN5","ERG","ESAM",
  "PTPRC","LST1","TYROBP","C5AR1",
  "CD3D","CD3E","NKG7",
  "ACTA2","TAGLN","MYH11","CNN1","RGS5","PDGFRB",
  "COL1A1","COL1A2","COL3A1","DCN","LUM",
  "PPBP","PF4","HBB",
  "ACKR1","PLVAP","COL15A1",
  "CXCL12","EFNB2","GJA4","SOX17","HEY1",
  "PROX1","PDPN","CCL21",
  "SLCO4A1","PRDM15","PROM1",
  "ANGPT2","PGF","LOX","CTHRC1","THBS1",
  "SULF1","ITLN1","MGP","ELN","COL8A1","FN1"
))
dot_genes <- intersect(dot_genes, rownames(obj))

if (length(dot_genes) > 0) {
  p5 <- DotPlot(
    obj,
    features = dot_genes,
    group.by = "EC_recluster",
    assay = "RNA"
  ) +
    RotatedAxis() +
    theme_classic(base_size = 8) +
    ggtitle("EC recluster marker audit")

  save_white(p5, "05_EC_recluster_marker_DotPlot.png", 19, 9)
}

p6 <- ggplot(
  stage_cluster,
  aes(x = disease_stage, y = fraction, fill = EC_recluster)
) +
  geom_col() +
  facet_wrap(~tech) +
  theme_classic(base_size = 11) +
  labs(
    title = "Stage composition after stage-blind EC reclustering",
    subtitle = "Diagnostic only; stage was not used in EC selection or Harmony",
    x = NULL,
    y = "Fraction of candidate EC",
    fill = "EC recluster"
  )

save_white(p6, "06_stage_composition_DIAGNOSTIC_ONLY.png", 11, 6.5)

cluster_qc_long <- cluster_qc %>%
  select(
    EC_recluster,
    frac_canonical_EC_support,
    frac_inclusive_EC_support,
    frac_hard_nonEC_conflict,
    frac_mural_4plus,
    frac_fibro_3plus,
    frac_rescued_from_Step2A_NotCandidate
  ) %>%
  pivot_longer(
    -EC_recluster,
    names_to = "metric",
    values_to = "value"
  )

p7 <- ggplot(
  cluster_qc_long,
  aes(x = metric, y = factor(EC_recluster), fill = value)
) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(title = "EC recluster lineage-QC fractions")

save_white(p7, "07_EC_recluster_QC_heatmap.png", 12, 8)

# -------------------------------------------------------------------------
# 9. Save REVIEW object (not final frozen EC yet)
# -------------------------------------------------------------------------
saveRDS(
  obj,
  file.path(
    OBJ,
    "GSE218166_step2B_V3_candidate_EC_reclustered_review.rds"
  ),
  compress = FALSE
)

high_clusters <- cluster_qc$EC_recluster[
  cluster_qc$auto_review_flag == "HIGH_CONFIDENCE_EC_REVIEW"
]

possible_clusters <- cluster_qc$EC_recluster[
  cluster_qc$auto_review_flag %in%
    c("HIGH_CONFIDENCE_EC_REVIEW", "POSSIBLE_EC_REVIEW")
]

writeLines(
  colnames(obj)[obj$EC_recluster %in% high_clusters],
  file.path(TAB, "12_auto_high_confidence_cluster_cells.txt")
)

writeLines(
  colnames(obj)[obj$EC_recluster %in% possible_clusters],
  file.path(TAB, "13_auto_high_plus_possible_cluster_cells.txt")
)

# -------------------------------------------------------------------------
# 10. Summary
# -------------------------------------------------------------------------
n_rescued <- sum(
  audit$prefilter_keep & !audit$step2A_candidate
)

summary_lines <- c(
  "============================================================",
  "GSE218166 Step 2B V3 - STANDALONE HARMONY",
  "============================================================",
  paste0("Run time: ", format(Sys.time())),
  paste0("ALL Step2A cells examined: ", ncol(obj_all)),
  paste0("Stage-blind EC candidates retained: ", ncol(obj)),
  paste0("Cells rescued from Step2A Not_candidate: ", n_rescued),
  paste0("EC reclusters: ", length(unique(obj$EC_recluster))),
  paste0("Clustering reduction: ", reduction_use),
  "",
  "DESIGN RULES",
  "1. Disease stage was NOT used to select EC cells.",
  "2. Step2A A/B/C tiers were not treated as final truth.",
  "3. Absolute EC marker evidence can rescue Step2A false negatives.",
  "4. Mural/fibroblast programs are not automatic exclusions.",
  "5. Hard immune/T-NK/epithelial/platelet-RBC conflicts are excluded.",
  "6. Harmony corrects tech only, using the standalone embedding API.",
  "",
  "NEXT",
  "Upload GSE218166_step2B_V3_REVIEW_SMALL.tar.gz.",
  "Final EC freeze will be performed in Step2C before stage testing."
)

writeLines(
  summary_lines,
  file.path(OUT, "Step2B_V3_DECISION_SUMMARY.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(LOG, "sessionInfo.txt")
)

cat("\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\n\nCompleted:", format(Sys.time()), "\n")
