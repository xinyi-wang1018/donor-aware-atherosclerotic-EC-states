#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260813)

BASE <- "/root/GSE218166"
IN_RDS <- file.path(
  BASE,
  "step1_build_and_broad_cluster",
  "GSE218166_step1_allcells_clustered.rds"
)
OUT <- file.path(BASE, "step2A_strict_EC_identity_audit")
TAB <- file.path(OUT, "tables")
PLOT <- file.path(OUT, "plots")
OBJ <- file.path(OUT, "objects")
LOG <- file.path(OUT, "logs")

for (d in c(OUT, TAB, PLOT, OBJ, LOG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

sink(file.path(LOG, "step2A_console.log"), split = TRUE)
on.exit({
  try(sink(), silent = TRUE)
}, add = TRUE)

cat("============================================================\n")
cat("GSE218166 Step 2A: strict EC identity audit\n")
cat("Start:", format(Sys.time()), "\n")
cat("Input:", IN_RDS, "\n")
cat("============================================================\n\n")

if (!file.exists(IN_RDS)) {
  stop("Input RDS not found: ", IN_RDS)
}

obj <- readRDS(IN_RDS)
DefaultAssay(obj) <- "RNA"

cat("Cells:", ncol(obj), "\n")
cat("Genes:", nrow(obj), "\n")
cat("Assays:", paste(Assays(obj), collapse = ", "), "\n")
cat("Reductions:", paste(Reductions(obj), collapse = ", "), "\n\n")

# ------------------------------------------------------------------
# 1. Metadata audit
# ------------------------------------------------------------------
md <- obj@meta.data
writeLines(colnames(md), file.path(TAB, "01_object_metadata_columns.txt"))

required_md <- c(
  "donor_id", "coronary_position", "lesion_process",
  "disease_stage", "tech", "seurat_clusters"
)
missing_md <- setdiff(required_md, colnames(md))
if (length(missing_md) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_md, collapse = ", ")
  )
}

# Never use disease stage to define EC identity.
obj$donor_id <- as.character(obj$donor_id)
obj$coronary_position <- as.character(obj$coronary_position)
obj$lesion_process <- as.character(obj$lesion_process)
obj$disease_stage <- factor(
  as.character(obj$disease_stage),
  levels = c("Normal", "Early", "Advanced"),
  ordered = TRUE
)
obj$tech <- as.character(obj$tech)
obj$seurat_clusters <- as.character(obj$seurat_clusters)

# This reconstruction is deliberately conservative and includes pathology process.
# It is used for audit only until a true physical-segment identifier is verified.
obj$sample_unit_reconstructed <- paste(
  obj$donor_id,
  obj$coronary_position,
  obj$lesion_process,
  obj$tech,
  sep = "__"
)

sample_unit_audit <- obj@meta.data %>%
  count(
    sample_unit_reconstructed,
    donor_id,
    coronary_position,
    lesion_process,
    disease_stage,
    tech,
    name = "n_cells"
  ) %>%
  arrange(tech, donor_id, coronary_position, lesion_process)

write.csv(
  sample_unit_audit,
  file.path(TAB, "02_reconstructed_sample_unit_audit.csv"),
  row.names = FALSE
)

unit_stage_check <- sample_unit_audit %>%
  group_by(sample_unit_reconstructed) %>%
  summarise(
    n_stage = n_distinct(disease_stage),
    stages = paste(sort(unique(as.character(disease_stage))), collapse = ";"),
    .groups = "drop"
  )

write.csv(
  unit_stage_check,
  file.path(TAB, "03_reconstructed_unit_stage_consistency.csv"),
  row.names = FALSE
)

if (any(unit_stage_check$n_stage > 1)) {
  stop(
    "A reconstructed sample unit maps to more than one disease stage. ",
    "Do not continue before reviewing tables/03_reconstructed_unit_stage_consistency.csv"
  )
}

stage_tech <- obj@meta.data %>%
  count(tech, disease_stage, name = "n_cells") %>%
  arrange(tech, disease_stage)

write.csv(stage_tech, file.path(TAB, "04_cells_by_tech_stage.csv"), row.names = FALSE)

# Audit original CellInfo files. Do not automatically choose any unknown sample-ID column.
raw_cellinfo <- c(
  FL = file.path(BASE, "GSM6735811_CellInfo.csv.gz"),
  `3DE` = file.path(BASE, "GSM6753958_sepDE3_CellInfo.csv.gz")
)

raw_audit_rows <- list()
candidate_id_rows <- list()

for (nm in names(raw_cellinfo)) {
  fp <- raw_cellinfo[[nm]]
  if (!file.exists(fp)) {
    cat("WARNING: raw CellInfo missing:", fp, "\n")
    next
  }
  x <- read.csv(fp, check.names = FALSE)
  writeLines(
    colnames(x),
    file.path(TAB, paste0("05_", nm, "_raw_CellInfo_columns.txt"))
  )

  raw_audit_rows[[nm]] <- data.frame(
    tech = nm,
    n_rows = nrow(x),
    n_columns = ncol(x),
    stringsAsFactors = FALSE
  )

  id_like <- grep(
    "sample|segment|specimen|section|arter|slice|piece|block|source|id",
    colnames(x),
    ignore.case = TRUE,
    value = TRUE
  )

  if (length(id_like) > 0) {
    tmp <- lapply(id_like, function(cc) {
      vals <- x[[cc]]
      data.frame(
        tech = nm,
        column = cc,
        n_unique = length(unique(vals[!is.na(vals)])),
        example_values = paste(head(unique(as.character(vals[!is.na(vals)])), 8), collapse = " | "),
        stringsAsFactors = FALSE
      )
    })
    candidate_id_rows[[nm]] <- do.call(rbind, tmp)
  }
}

if (length(raw_audit_rows) > 0) {
  write.csv(
    do.call(rbind, raw_audit_rows),
    file.path(TAB, "06_raw_CellInfo_dimensions.csv"),
    row.names = FALSE
  )
}
if (length(candidate_id_rows) > 0) {
  write.csv(
    do.call(rbind, candidate_id_rows),
    file.path(TAB, "07_raw_CellInfo_candidate_ID_columns.csv"),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------
# 2. Retrieve RNA data robustly under Seurat v4/v5
# ------------------------------------------------------------------
get_rna_layer <- function(x, what = c("data", "counts")) {
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

rna_data <- get_rna_layer(obj, "data")
if (is.null(rna_data) || nrow(rna_data) == 0 || ncol(rna_data) == 0) {
  cat("RNA data layer absent; running NormalizeData on the Step1 object copy.\n")
  obj <- NormalizeData(
    obj,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  rna_data <- get_rna_layer(obj, "data")
}
if (is.null(rna_data)) stop("Could not obtain RNA normalized data layer.")

rna_counts <- get_rna_layer(obj, "counts")
if (is.null(rna_counts)) {
  cat("WARNING: counts layer unavailable; marker detection will use normalized data > 0.\n")
  rna_counts <- rna_data
}

# ------------------------------------------------------------------
# 3. Pre-specified identity and diagnostic marker panels
# ------------------------------------------------------------------
marker_sets <- list(
  EC_core = c(
    "PECAM1", "VWF", "KDR", "CDH5", "CLDN5", "ESAM",
    "ENG", "EMCN", "RAMP2", "EPAS1", "ERG", "FLT1"
  ),
  Immune_myeloid = c(
    "PTPRC", "LST1", "TYROBP", "FCER1G", "AIF1", "CTSS", "LYZ"
  ),
  T_NK = c(
    "CD3D", "CD3E", "TRBC1", "TRBC2", "NKG7", "GNLY"
  ),
  Mural = c(
    "ACTA2", "TAGLN", "MYH11", "CNN1", "RGS5", "CSPG4", "PDGFRB", "NOTCH3"
  ),
  Fibroblast = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "COL6A1"
  ),
  Epithelial = c(
    "EPCAM", "KRT8", "KRT18", "KRT19", "KRT14"
  ),
  Platelet_RBC = c(
    "PPBP", "PF4", "NRGN", "HBB", "HBA1", "HBA2"
  ),
  Arterial_diagnostic = c(
    "CXCL12", "EFNB2", "GJA4", "GJA5", "SOX17", "HEY1", "SEMA3G", "BMX"
  ),
  Lymphatic_diagnostic = c(
    "PROX1", "PDPN", "CCL21", "LYVE1", "FLT4"
  ),
  Author_EC5_diagnostic = c(
    "SLCO4A1", "PRDM15", "PROM1"
  ),
  Project_early_activation = c(
    "ANGPT2", "CCL14", "PGF", "LOX", "CTHRC1", "THBS1"
  ),
  Project_late_remodeling = c(
    "SULF1", "ITLN1", "MGP", "ELN", "COL8A1", "FN1",
    "POSTN", "SPARC", "COL6A2", "IGFBP5"
  )
)

present_rows <- do.call(rbind, lapply(names(marker_sets), function(nm) {
  g <- marker_sets[[nm]]
  data.frame(
    marker_set = nm,
    gene = g,
    present = g %in% rownames(rna_data),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  present_rows,
  file.path(TAB, "08_marker_panel_presence.csv"),
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

scores <- data.frame(
  cell = colnames(obj),
  ec_core_mean = mean_expr(rna_data, marker_sets$EC_core),
  ec_core_detected = n_detected(rna_counts, marker_sets$EC_core),
  immune_myeloid_mean = mean_expr(rna_data, marker_sets$Immune_myeloid),
  tnk_mean = mean_expr(rna_data, marker_sets$T_NK),
  mural_mean = mean_expr(rna_data, marker_sets$Mural),
  fibro_mean = mean_expr(rna_data, marker_sets$Fibroblast),
  epithelial_mean = mean_expr(rna_data, marker_sets$Epithelial),
  platelet_rbc_mean = mean_expr(rna_data, marker_sets$Platelet_RBC),
  arterial_mean = mean_expr(rna_data, marker_sets$Arterial_diagnostic),
  lymphatic_mean = mean_expr(rna_data, marker_sets$Lymphatic_diagnostic),
  author_ec5_mean = mean_expr(rna_data, marker_sets$Author_EC5_diagnostic),
  early_activation_mean = mean_expr(rna_data, marker_sets$Project_early_activation),
  late_remodeling_mean = mean_expr(rna_data, marker_sets$Project_late_remodeling),
  stringsAsFactors = FALSE
)

scores <- cbind(
  scores,
  obj@meta.data[scores$cell, c(
    "donor_id", "coronary_position", "lesion_process",
    "disease_stage", "tech", "seurat_clusters",
    "sample_unit_reconstructed"
  ), drop = FALSE]
)

# ------------------------------------------------------------------
# 4. Technology-stratified ranks
#    IMPORTANT: disease_stage is NOT used in any identity threshold.
# ------------------------------------------------------------------
percent_rank <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (sum(ok) <= 1) {
    out[ok] <- 1
    return(out)
  }
  out[ok] <- rank(x[ok], ties.method = "average") / sum(ok)
  out
}

scores <- scores %>%
  group_by(tech) %>%
  mutate(
    ec_pct = percent_rank(ec_core_mean),
    immune_pct = percent_rank(immune_myeloid_mean),
    tnk_pct = percent_rank(tnk_mean),
    mural_pct = percent_rank(mural_mean),
    fibro_pct = percent_rank(fibro_mean),
    epithelial_pct = percent_rank(epithelial_mean),
    platelet_rbc_pct = percent_rank(platelet_rbc_mean)
  ) %>%
  ungroup()

scores$hard_nonEC_pct <- pmax(
  scores$immune_pct,
  scores$tnk_pct,
  scores$epithelial_pct,
  scores$platelet_rbc_pct,
  na.rm = TRUE
)

scores$all_competitor_pct <- pmax(
  scores$immune_pct,
  scores$tnk_pct,
  scores$mural_pct,
  scores$fibro_pct,
  scores$epithelial_pct,
  scores$platelet_rbc_pct,
  na.rm = TRUE
)

scores$ec_vs_hard_nonEC_margin <- scores$ec_pct - scores$hard_nonEC_pct
scores$ec_vs_all_competitor_margin <- scores$ec_pct - scores$all_competitor_pct

# These are REVIEW TIERS, not the final biological freeze.
# A is intentionally conservative.
tier_A <- with(
  scores,
  ec_pct >= 0.85 &
    ec_core_detected >= 5 &
    ec_vs_all_competitor_margin >= 0.05 &
    hard_nonEC_pct < 0.90
)

# B preserves potential EC-to-mesenchymal transition candidates but still excludes
# cells dominated by obvious immune/epithelial/platelet programs.
tier_B <- with(
  scores,
  !tier_A &
    ec_pct >= 0.75 &
    ec_core_detected >= 4 &
    ec_vs_hard_nonEC_margin >= 0.05 &
    hard_nonEC_pct < 0.95
)

# C is deliberately permissive and is only for sensitivity/manual review.
tier_C <- with(
  scores,
  !tier_A & !tier_B &
    ec_pct >= 0.70 &
    ec_core_detected >= 3 &
    ec_vs_hard_nonEC_margin > -0.05 &
    hard_nonEC_pct < 0.98
)

scores$EC_review_tier <- "Not_candidate"
scores$EC_review_tier[tier_C] <- "C_lineage_candidate"
scores$EC_review_tier[tier_B] <- "B_lineage_supported"
scores$EC_review_tier[tier_A] <- "A_ultra_strict"

scores$EC_review_tier <- factor(
  scores$EC_review_tier,
  levels = c(
    "Not_candidate",
    "C_lineage_candidate",
    "B_lineage_supported",
    "A_ultra_strict"
  )
)

# Transfer scores into object metadata.
add_cols <- c(
  "ec_core_mean", "ec_core_detected",
  "immune_myeloid_mean", "tnk_mean", "mural_mean", "fibro_mean",
  "epithelial_mean", "platelet_rbc_mean",
  "arterial_mean", "lymphatic_mean", "author_ec5_mean",
  "early_activation_mean", "late_remodeling_mean",
  "ec_pct", "immune_pct", "tnk_pct", "mural_pct", "fibro_pct",
  "epithelial_pct", "platelet_rbc_pct",
  "ec_vs_hard_nonEC_margin", "ec_vs_all_competitor_margin",
  "EC_review_tier"
)

for (cc in add_cols) {
  obj[[cc]] <- scores[[cc]][match(colnames(obj), scores$cell)]
}

write.csv(
  scores,
  gzfile(file.path(TAB, "09_cell_level_EC_identity_scores.csv.gz")),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 5. Summaries independent of stage selection
# ------------------------------------------------------------------
tier_summary <- scores %>%
  count(tech, EC_review_tier, name = "n_cells") %>%
  group_by(tech) %>%
  mutate(fraction_within_tech = n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(
  tier_summary,
  file.path(TAB, "10_EC_review_tier_counts_by_tech.csv"),
  row.names = FALSE
)

stage_tier_summary <- scores %>%
  count(tech, disease_stage, EC_review_tier, name = "n_cells") %>%
  group_by(tech, disease_stage) %>%
  mutate(fraction_within_tech_stage = n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(
  stage_tier_summary,
  file.path(TAB, "11_EC_review_tier_counts_by_stage_tech.csv"),
  row.names = FALSE
)

cluster_summary <- scores %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells = n(),
    n_tech = n_distinct(tech),
    median_ec_core = median(ec_core_mean, na.rm = TRUE),
    median_ec_pct = median(ec_pct, na.rm = TRUE),
    median_ec_detected = median(ec_core_detected, na.rm = TRUE),
    median_immune_pct = median(immune_pct, na.rm = TRUE),
    median_mural_pct = median(mural_pct, na.rm = TRUE),
    median_fibro_pct = median(fibro_pct, na.rm = TRUE),
    median_ec_vs_all_margin = median(ec_vs_all_competitor_margin, na.rm = TRUE),
    frac_A = mean(EC_review_tier == "A_ultra_strict"),
    frac_B = mean(EC_review_tier == "B_lineage_supported"),
    frac_C = mean(EC_review_tier == "C_lineage_candidate"),
    frac_any_candidate = mean(EC_review_tier != "Not_candidate"),
    .groups = "drop"
  ) %>%
  arrange(desc(frac_A), desc(frac_any_candidate), desc(median_ec_pct))

write.csv(
  cluster_summary,
  file.path(TAB, "12_cluster_EC_identity_summary.csv"),
  row.names = FALSE
)

candidate_counts <- scores %>%
  filter(EC_review_tier != "Not_candidate") %>%
  count(
    donor_id, coronary_position, lesion_process,
    disease_stage, tech, EC_review_tier,
    name = "n_candidate_cells"
  ) %>%
  arrange(tech, donor_id, disease_stage)

write.csv(
  candidate_counts,
  file.path(TAB, "13_candidate_EC_counts_by_donor_segment_stage.csv"),
  row.names = FALSE
)

# Candidate cell lists
writeLines(
  scores$cell[scores$EC_review_tier == "A_ultra_strict"],
  file.path(TAB, "14_A_ultra_strict_cell_barcodes.txt")
)
writeLines(
  scores$cell[scores$EC_review_tier %in% c("A_ultra_strict", "B_lineage_supported")],
  file.path(TAB, "15_AB_candidate_cell_barcodes.txt")
)
writeLines(
  scores$cell[scores$EC_review_tier != "Not_candidate"],
  file.path(TAB, "16_ABC_candidate_cell_barcodes.txt")
)

# ------------------------------------------------------------------
# 6. Diagnostic plots (white background)
# ------------------------------------------------------------------
save_white <- function(p, path, width = 10, height = 7) {
  ggsave(
    path, p,
    width = width, height = height,
    dpi = 300, bg = "white"
  )
}

reduction_use <- if ("umap.harmony" %in% Reductions(obj)) {
  "umap.harmony"
} else if ("umap" %in% Reductions(obj)) {
  "umap"
} else if ("umap.pca" %in% Reductions(obj)) {
  "umap.pca"
} else {
  stop("No UMAP reduction found.")
}

p_tier <- DimPlot(
  obj,
  reduction = reduction_use,
  group.by = "EC_review_tier",
  raster = TRUE,
  pt.size = 0.25
) +
  ggtitle("GSE218166 Step2A provisional EC review tiers") +
  theme_classic(base_size = 11)

save_white(p_tier, file.path(PLOT, "01_UMAP_EC_review_tiers.png"), 11, 8)

p_cluster <- DimPlot(
  obj,
  reduction = reduction_use,
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 0.25
) +
  ggtitle("Broad clusters for EC identity audit") +
  theme_classic(base_size = 11)

save_white(p_cluster, file.path(PLOT, "02_UMAP_broad_clusters_review.png"), 11, 8)

p_ec <- FeaturePlot(
  obj,
  reduction = reduction_use,
  features = c("ec_pct", "ec_vs_all_competitor_margin"),
  ncol = 2,
  raster = TRUE
) & theme_classic(base_size = 10)

save_white(p_ec, file.path(PLOT, "03_UMAP_EC_rank_and_margin.png"), 13, 6)

p_comp <- FeaturePlot(
  obj,
  reduction = reduction_use,
  features = c("immune_pct", "mural_pct", "fibro_pct", "platelet_rbc_pct"),
  ncol = 2,
  raster = TRUE
) & theme_classic(base_size = 10)

save_white(p_comp, file.path(PLOT, "04_UMAP_competing_lineage_ranks.png"), 13, 11)

score_long <- scores %>%
  select(
    tech, ec_pct, immune_pct, mural_pct, fibro_pct,
    epithelial_pct, platelet_rbc_pct
  ) %>%
  tidyr::pivot_longer(
    cols = -tech,
    names_to = "score",
    values_to = "percentile"
  )

p_density <- ggplot(
  score_long,
  aes(x = percentile, linetype = tech)
) +
  geom_density(linewidth = 0.8) +
  facet_wrap(~score, ncol = 3, scales = "free_y") +
  theme_classic(base_size = 10) +
  theme(panel.grid = element_blank()) +
  labs(
    title = "Technology-stratified lineage score percentiles",
    x = "Within-technology percentile",
    y = "Density"
  )

save_white(p_density, file.path(PLOT, "05_lineage_rank_density_by_tech.png"), 12, 9)

dot_genes <- unique(c(
  "PECAM1", "VWF", "KDR", "CDH5", "CLDN5", "ERG",
  "PTPRC", "LST1", "TYROBP",
  "ACTA2", "TAGLN", "MYH11", "RGS5",
  "COL1A1", "COL3A1", "DCN", "LUM",
  "PPBP", "PF4",
  "CXCL12", "EFNB2", "GJA4", "SOX17", "HEY1",
  "ANGPT2", "PGF", "LOX", "CTHRC1", "THBS1",
  "SULF1", "SLCO4A1", "PRDM15"
))
dot_genes <- intersect(dot_genes, rownames(obj))

if (length(dot_genes) > 0) {
  p_dot <- DotPlot(
    obj,
    features = dot_genes,
    group.by = "EC_review_tier",
    assay = "RNA"
  ) +
    RotatedAxis() +
    theme_classic(base_size = 9) +
    ggtitle("Marker audit across provisional EC review tiers")

  save_white(p_dot, file.path(PLOT, "06_marker_DotPlot_by_EC_review_tier.png"), 16, 6.5)
}

# Cluster × marker diagnostics
cluster_dot_genes <- intersect(
  c(
    "PECAM1", "VWF", "KDR", "CDH5", "CLDN5",
    "PTPRC", "LST1", "TYROBP",
    "ACTA2", "TAGLN", "MYH11", "RGS5",
    "COL1A1", "COL3A1", "DCN", "LUM",
    "PPBP", "PF4",
    "CXCL12", "EFNB2", "GJA4", "SOX17", "HEY1",
    "SLCO4A1", "SULF1"
  ),
  rownames(obj)
)

if (length(cluster_dot_genes) > 0) {
  p_cdot <- DotPlot(
    obj,
    features = cluster_dot_genes,
    group.by = "seurat_clusters",
    assay = "RNA"
  ) +
    RotatedAxis() +
    theme_classic(base_size = 8) +
    ggtitle("Broad-cluster EC/non-EC marker audit")

  save_white(p_cdot, file.path(PLOT, "07_marker_DotPlot_by_broad_cluster.png"), 16, 10)
}

# Stage is shown only AFTER identity tiers are assigned.
p_stage_tier <- ggplot(
  stage_tier_summary,
  aes(x = disease_stage, y = fraction_within_tech_stage, fill = EC_review_tier)
) +
  geom_col() +
  facet_wrap(~tech) +
  theme_classic(base_size = 11) +
  labs(
    title = "Provisional EC-tier composition by stage (diagnostic only)",
    x = NULL,
    y = "Fraction of all cells",
    fill = "EC review tier"
  )

save_white(p_stage_tier, file.path(PLOT, "08_EC_tier_fraction_by_stage_tech.png"), 10, 6)

# ------------------------------------------------------------------
# 7. Save provisional objects; not yet the frozen final EC object
# ------------------------------------------------------------------
saveRDS(
  obj,
  file.path(OBJ, "GSE218166_step2A_allcells_with_EC_identity_scores.rds"),
  compress = FALSE
)

ab_cells <- scores$cell[
  scores$EC_review_tier %in% c("A_ultra_strict", "B_lineage_supported")
]
abc_cells <- scores$cell[scores$EC_review_tier != "Not_candidate"]

if (length(ab_cells) > 0) {
  saveRDS(
    subset(obj, cells = ab_cells),
    file.path(OBJ, "GSE218166_step2A_AB_candidate_EC.rds"),
    compress = FALSE
  )
}
if (length(abc_cells) > 0) {
  saveRDS(
    subset(obj, cells = abc_cells),
    file.path(OBJ, "GSE218166_step2A_ABC_candidate_EC.rds"),
    compress = FALSE
  )
}

# ------------------------------------------------------------------
# 8. Human-readable decision report
# ------------------------------------------------------------------
n_units <- nrow(sample_unit_audit)
nA <- sum(scores$EC_review_tier == "A_ultra_strict")
nB <- sum(scores$EC_review_tier == "B_lineage_supported")
nC <- sum(scores$EC_review_tier == "C_lineage_candidate")

summary_lines <- c(
  "============================================================",
  "GSE218166 Step 2A strict EC identity audit summary",
  "============================================================",
  paste0("Run time: ", format(Sys.time())),
  paste0("Cells: ", ncol(obj)),
  paste0("Genes: ", nrow(obj)),
  paste0("UMAP used: ", reduction_use),
  "",
  "IMPORTANT METHOD RULES",
  "1. Disease stage was NOT used to define EC identity.",
  "2. Broad cluster identity alone was NOT used to freeze EC.",
  "3. A/B/C are provisional review tiers, not final frozen EC labels.",
  "4. Tier A is intentionally conservative.",
  "5. Tier B/C retain possible EC-lineage transition cells for sensitivity review.",
  "",
  "PROVISIONAL TIER COUNTS",
  paste0("A_ultra_strict: ", nA),
  paste0("B_lineage_supported: ", nB),
  paste0("C_lineage_candidate: ", nC),
  paste0("Any A/B/C candidate: ", nA + nB + nC),
  "",
  "SAMPLE-UNIT AUDIT",
  paste0("Reconstructed donor-position-process-tech units: ", n_units),
  "These reconstructed units are NOT automatically accepted as the final pseudobulk units.",
  "Review raw CellInfo candidate ID columns before Step 3.",
  "",
  "NEXT DECISION",
  "Review 01-08 plots and tables 08-13.",
  "Then freeze the main high-confidence EC set in Step 2B.",
  "Only after that should Normal -> Early -> Advanced testing begin.",
  ""
)

writeLines(
  summary_lines,
  file.path(OUT, "Step2A_DECISION_SUMMARY.txt")
)

capture.output(sessionInfo(), file = file.path(LOG, "sessionInfo.txt"))

cat("\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\nCompleted:", format(Sys.time()), "\n")
