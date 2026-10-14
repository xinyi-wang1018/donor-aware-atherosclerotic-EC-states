#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(20260813)

BASE <- "/root/GSE218166"

# Re-open ALL cells. Step2D is a fresh, stage-blind lineage audit.
IN_RDS <- file.path(
  BASE,
  "step2A_strict_EC_identity_audit",
  "objects",
  "GSE218166_step2A_allcells_with_EC_identity_scores.rds"
)

OUT  <- file.path(BASE, "step2D_AUTHOR_ALIGNED_EC_AUDIT")
TAB  <- file.path(OUT, "tables")
PLOT <- file.path(OUT, "plots")
OBJ  <- file.path(OUT, "objects")
LOG  <- file.path(OUT, "logs")

for (d in c(OUT, TAB, PLOT, OBJ, LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

sink(file.path(LOG, "step2D_console.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

cat("============================================================\n")
cat("GSE218166 STEP 2D - AUTHOR-ALIGNED EC LINEAGE AUDIT\n")
cat("Technology-stratified relative-expression scoring\n")
cat("Disease stage is NOT used for cell identity.\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

if (!file.exists(IN_RDS)) stop("Input not found: ", IN_RDS)

obj <- readRDS(IN_RDS)
DefaultAssay(obj) <- "RNA"

required_md <- c(
  "tech", "seurat_clusters", "donor_id", "disease_stage",
  "coronary_position", "lesion_process", "plate"
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

obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ----------------------------------------------------------------------
# 1. Marker panels
#
# Author-aligned major-type anchors:
# EC: TIE1 / PECAM1 (Extended Data Fig.1)
# VSMC: CNN1 / MYH11
# Myeloid: CD14 / C5AR1
# Lymphoid: CD69 / CD2
#
# Additional canonical markers are used to stabilize the lineage score.
# ----------------------------------------------------------------------
modules <- list(
  EC_author_anchor = c("TIE1", "PECAM1"),
  EC_canonical = c(
    "TIE1","PECAM1","VWF","KDR","CDH5","CLDN5",
    "ESAM","EMCN","FLT1","TEK","PTPRB","ADGRL4","ERG"
  ),

  VSMC_author_anchor = c("CNN1", "MYH11"),
  VSMC_contractile = c(
    "CNN1","MYH11","ACTA2","TAGLN","RGS5","CSPG4",
    "PDGFRB","NOTCH3","MYOCD"
  ),

  Fibroblast = c(
    "COL1A1","COL1A2","COL3A1","DCN","LUM",
    "PDGFRA","DPT","PI16","MFAP5"
  ),

  Myeloid_author_anchor = c("CD14","C5AR1"),
  Myeloid = c(
    "CD14","C5AR1","PTPRC","LST1","TYROBP",
    "FCER1G","AIF1","CTSS","LYZ"
  ),

  Lymphoid_author_anchor = c("CD69","CD2"),
  Lymphoid = c(
    "CD69","CD2","CD3D","CD3E","TRBC1","TRBC2",
    "NKG7","GNLY"
  ),

  Platelet_RBC = c(
    "PPBP","PF4","HBB","HBA1","HBA2"
  ),

  EC1_dflow_candidate = c(
    "KLF2","KLF4","THBD","CXCL12","EFNB2","GJA4","GJA5",
    "SOX17","HEY1","SEMA3G"
  ),
  EC2_ACKR1_candidate = c(
    "ACKR1","PLVAP","VWF","EMCN","ADGRL4"
  ),
  EC4_EndMT_candidate = c(
    "FN1","COL8A1","COL3A1","BGN","MMP2","VIM",
    "ACTA2","TAGLN"
  ),
  EC5_SLCO4A1_candidate = c(
    "SLCO4A1","PRDM15","PROM1"
  ),

  Project_early = c(
    "ANGPT2","CCL14","PGF","LOX","CTHRC1","THBS1"
  ),
  Project_late = c(
    "SULF1","ITLN1","MGP","ELN","COL8A1",
    "FN1","SPARC","COL6A2","IGFBP5"
  )
)

presence <- bind_rows(lapply(names(modules), function(nm) {
  data.frame(
    module = nm,
    gene = modules[[nm]],
    present = modules[[nm]] %in% rownames(obj),
    stringsAsFactors = FALSE
  )
}))

write.csv(
  presence,
  file.path(TAB, "01_marker_presence.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 2. Technology-stratified gene z-scores
#
# Why:
# FL and 3DE have very different detection depth.
# "gene detected / not detected" is therefore not a fair lineage criterion.
#
# We standardize each marker gene WITHIN technology, then average marker
# z-scores into lineage scores.
# ----------------------------------------------------------------------
all_marker_genes <- unique(unlist(modules))
all_marker_genes <- intersect(all_marker_genes, rownames(obj))

if (length(all_marker_genes) < 10) {
  stop("Too few marker genes found.")
}

techs <- unique(as.character(obj$tech))
score_parts <- list()

for (tt in techs) {
  cat("Scoring technology:", tt, "\n")

  cells_tt <- colnames(obj)[as.character(obj$tech) == tt]
  x <- subset(obj, cells = cells_tt)
  DefaultAssay(x) <- "RNA"

  # Scale only the curated marker genes, not the whole transcriptome.
  x <- ScaleData(
    x,
    features = all_marker_genes,
    assay = "RNA",
    do.center = TRUE,
    do.scale = TRUE,
    verbose = FALSE
  )

  scaled <- GetAssayData(x, assay = "RNA", layer = "scale.data")

  module_score <- function(gs) {
    gg <- intersect(gs, rownames(scaled))
    if (length(gg) == 0) return(rep(NA_real_, ncol(scaled)))
    Matrix::colMeans(scaled[gg, , drop = FALSE])
  }

  z <- data.frame(
    cell = colnames(x),
    tech = tt,
    stringsAsFactors = FALSE
  )

  for (nm in names(modules)) {
    z[[paste0(nm, "_z")]] <- module_score(modules[[nm]])
  }

  score_parts[[tt]] <- z
}

scores <- bind_rows(score_parts)

# Metadata
meta_cols <- c(
  "seurat_clusters","donor_id","disease_stage",
  "coronary_position","lesion_process","plate"
)

scores <- cbind(
  scores,
  obj@meta.data[scores$cell, meta_cols, drop = FALSE]
)

# ----------------------------------------------------------------------
# 3. Relative lineage dominance
# ----------------------------------------------------------------------
scores <- scores %>%
  mutate(
    EC_score = EC_canonical_z,

    nonEC_max = pmax(
      VSMC_contractile_z,
      Fibroblast_z,
      Myeloid_z,
      Lymphoid_z,
      Platelet_RBC_z,
      na.rm = TRUE
    ),

    EC_margin = EC_score - nonEC_max,

    # Winner-takes-highest-score is for DIAGNOSTIC annotation only.
    relative_major_type = case_when(
      EC_score >= VSMC_contractile_z &
        EC_score >= Fibroblast_z &
        EC_score >= Myeloid_z &
        EC_score >= Lymphoid_z &
        EC_score >= Platelet_RBC_z ~ "EC_relative_winner",

      VSMC_contractile_z >= Fibroblast_z &
        VSMC_contractile_z >= Myeloid_z &
        VSMC_contractile_z >= Lymphoid_z &
        VSMC_contractile_z >= Platelet_RBC_z ~ "VSMC_relative_winner",

      Fibroblast_z >= Myeloid_z &
        Fibroblast_z >= Lymphoid_z &
        Fibroblast_z >= Platelet_RBC_z ~ "Fibroblast_relative_winner",

      Myeloid_z >= Lymphoid_z &
        Myeloid_z >= Platelet_RBC_z ~ "Myeloid_relative_winner",

      Lymphoid_z >= Platelet_RBC_z ~ "Lymphoid_relative_winner",
      TRUE ~ "Platelet_RBC_relative_winner"
    ),

    # Conservative audit tiers based on relative expression, not raw detection.
    author_aligned_EC_tier = case_when(
      EC_author_anchor_z > 0 &
        EC_canonical_z > 0.25 &
        EC_margin > 0.25
        ~ "A_relative_strict_EC",

      EC_canonical_z > 0 &
        EC_margin > 0
        ~ "B_relative_EC",

      EC_canonical_z > 0 &
        (
          EC4_EndMT_candidate_z > 0.5 |
          EC5_SLCO4A1_candidate_z > 0.5
        )
        ~ "C_remodeling_EC_candidate",

      TRUE ~ "Not_EC_by_relative_score"
    )
  )

write.csv(
  scores,
  gzfile(file.path(TAB, "02_cell_level_author_aligned_scores.csv.gz")),
  row.names = FALSE
)

# Add score fields into object
add_cols <- setdiff(
  colnames(scores),
  c("cell","tech","seurat_clusters","donor_id","disease_stage",
    "coronary_position","lesion_process","plate")
)

for (cc in add_cols) {
  obj[[cc]] <- scores[[cc]][match(colnames(obj), scores$cell)]
}

# ----------------------------------------------------------------------
# 4. Broad cluster audit, technology-specific and pooled
# ----------------------------------------------------------------------
cluster_tech <- scores %>%
  group_by(tech, seurat_clusters) %>%
  summarise(
    n_cells = n(),

    median_EC_score = median(EC_score, na.rm = TRUE),
    median_EC_anchor = median(EC_author_anchor_z, na.rm = TRUE),
    median_VSMC = median(VSMC_contractile_z, na.rm = TRUE),
    median_FB = median(Fibroblast_z, na.rm = TRUE),
    median_myeloid = median(Myeloid_z, na.rm = TRUE),
    median_lymphoid = median(Lymphoid_z, na.rm = TRUE),
    median_EC_margin = median(EC_margin, na.rm = TRUE),

    frac_EC_relative_winner =
      mean(relative_major_type == "EC_relative_winner"),
    frac_A =
      mean(author_aligned_EC_tier == "A_relative_strict_EC"),
    frac_AB =
      mean(author_aligned_EC_tier %in%
             c("A_relative_strict_EC","B_relative_EC")),
    frac_remodel_candidate =
      mean(author_aligned_EC_tier == "C_remodeling_EC_candidate"),

    median_EC1 = median(EC1_dflow_candidate_z, na.rm = TRUE),
    median_EC2_ACKR1 = median(EC2_ACKR1_candidate_z, na.rm = TRUE),
    median_EC4_EndMT = median(EC4_EndMT_candidate_z, na.rm = TRUE),
    median_EC5 = median(EC5_SLCO4A1_candidate_z, na.rm = TRUE),

    median_project_early = median(Project_early_z, na.rm = TRUE),
    median_project_late = median(Project_late_z, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  arrange(
    tech,
    desc(median_EC_margin),
    desc(frac_EC_relative_winner)
  )

write.csv(
  cluster_tech,
  file.path(TAB, "03_broad_cluster_relative_lineage_audit_by_tech.csv"),
  row.names = FALSE
)

cluster_pooled <- scores %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells = n(),
    n_tech = n_distinct(tech),

    median_EC_score = median(EC_score, na.rm = TRUE),
    median_EC_anchor = median(EC_author_anchor_z, na.rm = TRUE),
    median_VSMC = median(VSMC_contractile_z, na.rm = TRUE),
    median_FB = median(Fibroblast_z, na.rm = TRUE),
    median_myeloid = median(Myeloid_z, na.rm = TRUE),
    median_lymphoid = median(Lymphoid_z, na.rm = TRUE),
    median_EC_margin = median(EC_margin, na.rm = TRUE),

    frac_EC_relative_winner =
      mean(relative_major_type == "EC_relative_winner"),
    frac_A =
      mean(author_aligned_EC_tier == "A_relative_strict_EC"),
    frac_AB =
      mean(author_aligned_EC_tier %in%
             c("A_relative_strict_EC","B_relative_EC")),

    median_EC1 = median(EC1_dflow_candidate_z, na.rm = TRUE),
    median_EC2_ACKR1 = median(EC2_ACKR1_candidate_z, na.rm = TRUE),
    median_EC4_EndMT = median(EC4_EndMT_candidate_z, na.rm = TRUE),
    median_EC5 = median(EC5_SLCO4A1_candidate_z, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  arrange(
    desc(median_EC_margin),
    desc(frac_EC_relative_winner)
  )

write.csv(
  cluster_pooled,
  file.path(TAB, "04_broad_cluster_relative_lineage_audit_pooled.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 5. Special audit of clusters previously implicated by EC markers
# ----------------------------------------------------------------------
ec_focus_clusters <- unique(as.character(
  head(cluster_pooled$seurat_clusters, 10)
))

focus <- cluster_tech %>%
  filter(as.character(seurat_clusters) %in% ec_focus_clusters)

write.csv(
  focus,
  file.path(TAB, "05_top10_EC_like_clusters_detailed_audit.csv"),
  row.names = FALSE
)

# Stage is viewed ONLY after the identity scores were computed.
stage_by_cluster <- scores %>%
  filter(as.character(seurat_clusters) %in% ec_focus_clusters) %>%
  count(
    tech,
    seurat_clusters,
    disease_stage,
    name = "n_cells"
  ) %>%
  group_by(tech, seurat_clusters) %>%
  mutate(frac_within_cluster = n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(
  stage_by_cluster,
  file.path(TAB, "06_top_EC_like_cluster_stage_composition_DIAGNOSTIC_ONLY.csv"),
  row.names = FALSE
)

tier_stage <- scores %>%
  count(
    tech,
    disease_stage,
    author_aligned_EC_tier,
    name = "n_cells"
  )

write.csv(
  tier_stage,
  file.path(TAB, "07_relative_EC_tier_counts_by_stage_tech_DIAGNOSTIC_ONLY.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 6. Co-expression conflict audit
#
# This asks whether canonical EC and contractile-VSMC signals occur in
# the same cells, or whether cluster averages merely mix two populations.
# ----------------------------------------------------------------------
data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")

det <- function(g) {
  if (!g %in% rownames(data_mat)) {
    return(rep(FALSE, ncol(data_mat)))
  }
  as.numeric(data_mat[g, ]) > 0
}

coexpr <- data.frame(
  cell = colnames(obj),
  PECAM1 = det("PECAM1"),
  TIE1 = det("TIE1"),
  CDH5 = det("CDH5"),
  CLDN5 = det("CLDN5"),
  KDR = det("KDR"),
  MYH11 = det("MYH11"),
  CNN1 = det("CNN1"),
  ACTA2 = det("ACTA2"),
  RGS5 = det("RGS5"),
  stringsAsFactors = FALSE
)

coexpr <- cbind(
  coexpr,
  obj@meta.data[
    coexpr$cell,
    c("tech","seurat_clusters","disease_stage"),
    drop = FALSE
  ]
)

coexpr <- coexpr %>%
  mutate(
    EC_anchor_any = PECAM1 | TIE1,
    EC_key_3plus =
      (PECAM1 + TIE1 + CDH5 + CLDN5 + KDR) >= 3,
    VSMC_anchor_any = MYH11 | CNN1,
    VSMC_contractile_3plus =
      (MYH11 + CNN1 + ACTA2 + RGS5) >= 3,
    EC_VSMC_double_program =
      EC_key_3plus & VSMC_contractile_3plus
  )

coexpr_summary <- coexpr %>%
  group_by(tech, seurat_clusters) %>%
  summarise(
    n_cells = n(),
    pct_EC_anchor_any = mean(EC_anchor_any) * 100,
    pct_EC_key_3plus = mean(EC_key_3plus) * 100,
    pct_VSMC_anchor_any = mean(VSMC_anchor_any) * 100,
    pct_VSMC_contractile_3plus =
      mean(VSMC_contractile_3plus) * 100,
    pct_EC_VSMC_double_program =
      mean(EC_VSMC_double_program) * 100,
    .groups = "drop"
  ) %>%
  arrange(
    tech,
    desc(pct_EC_key_3plus)
  )

write.csv(
  coexpr_summary,
  file.path(TAB, "08_EC_VSMC_single_cell_coexpression_audit.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 7. Figures - pure white background
# ----------------------------------------------------------------------
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

reduction_use <- if ("umap.harmony" %in% Reductions(obj)) {
  "umap.harmony"
} else if ("umap" %in% Reductions(obj)) {
  "umap"
} else if ("umap.pca" %in% Reductions(obj)) {
  "umap.pca"
} else {
  NULL
}

if (!is.null(reduction_use)) {
  p1 <- FeaturePlot(
    obj,
    reduction = reduction_use,
    features = c(
      "EC_canonical_z",
      "VSMC_contractile_z",
      "Fibroblast_z",
      "EC_margin"
    ),
    ncol = 2,
    raster = TRUE
  ) & theme_classic(base_size = 9)

  save_white(
    p1,
    "01_author_aligned_relative_lineage_scores_UMAP.png",
    13,
    10
  )

  p2 <- DimPlot(
    obj,
    reduction = reduction_use,
    group.by = "author_aligned_EC_tier",
    raster = TRUE,
    pt.size = 0.25
  ) +
    ggtitle("Author-aligned relative EC audit") +
    theme_classic(base_size = 10)

  save_white(
    p2,
    "02_author_aligned_EC_tiers_UMAP.png",
    11,
    8
  )
}

# Cluster relative-score heatmap
heat_cols <- c(
  "median_EC_anchor","median_EC_score",
  "median_VSMC","median_FB",
  "median_myeloid","median_lymphoid",
  "median_EC_margin",
  "median_EC1","median_EC2_ACKR1",
  "median_EC4_EndMT","median_EC5"
)

heat_df <- cluster_tech %>%
  select(tech, seurat_clusters, all_of(heat_cols)) %>%
  pivot_longer(
    cols = all_of(heat_cols),
    names_to = "metric",
    values_to = "value"
  )

p3 <- ggplot(
  heat_df,
  aes(
    x = metric,
    y = factor(seurat_clusters),
    fill = value
  )
) +
  geom_tile() +
  facet_wrap(~tech, scales = "free_y") +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(
    title = "Technology-stratified relative lineage scores by broad cluster"
  )

save_white(
  p3,
  "03_cluster_relative_lineage_score_heatmap.png",
  14,
  11
)

# EC vs VSMC relative score scatter
plot_sample <- scores
if (nrow(plot_sample) > 12000) {
  set.seed(20260813)
  plot_sample <- plot_sample[
    sample(seq_len(nrow(plot_sample)), 12000),
    ,
    drop = FALSE
  ]
}

p4 <- ggplot(
  plot_sample,
  aes(
    x = EC_canonical_z,
    y = VSMC_contractile_z,
    shape = tech
  )
) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_classic(base_size = 10) +
  labs(
    title = "EC versus VSMC relative program score",
    x = "EC canonical score (within-tech z)",
    y = "VSMC contractile score (within-tech z)"
  )

save_white(
  p4,
  "04_EC_vs_VSMC_relative_score_scatter.png",
  9,
  7
)

# Top EC-like broad clusters
top_plot <- cluster_pooled %>%
  slice_head(n = 15)

p5 <- ggplot(
  top_plot,
  aes(
    x = reorder(as.factor(seurat_clusters), median_EC_margin),
    y = median_EC_margin
  )
) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 10) +
  labs(
    title = "Broad clusters ranked by EC lineage margin",
    x = "Broad cluster",
    y = "Median EC minus strongest non-EC program"
  )

save_white(
  p5,
  "05_broad_clusters_ranked_by_EC_margin.png",
  9,
  7
)

# ----------------------------------------------------------------------
# 8. Save scored object for the next final freeze
# ----------------------------------------------------------------------
saveRDS(
  obj,
  file.path(
    OBJ,
    "GSE218166_step2D_allcells_author_aligned_lineage_scored.rds"
  ),
  compress = FALSE
)

# ----------------------------------------------------------------------
# 9. Summary
# ----------------------------------------------------------------------
top5 <- head(cluster_pooled, 5)

summary_lines <- c(
  "============================================================",
  "GSE218166 STEP2D AUTHOR-ALIGNED EC AUDIT",
  "============================================================",
  paste0("All cells audited: ", ncol(obj)),
  paste0("Technologies: ", paste(techs, collapse = ", ")),
  "",
  "WHY STEP2D WAS NEEDED",
  "Raw marker detection rates are not comparable between FL and 3DE.",
  "The prior Main_strict_EC set retained strong VSMC/fibroblast signals.",
  "The source paper itself distinguishes EC using TIE1/PECAM1 and VSMC using CNN1/MYH11.",
  "Therefore lineage identity is re-audited here using within-technology relative expression.",
  "",
  "IMPORTANT",
  "Disease stage was NOT used to calculate identity scores or choose EC candidates.",
  "This step is an audit, not the final disease-stage test.",
  "",
  "Top five broad clusters by pooled EC margin:",
  paste(
    paste0(
      top5$seurat_clusters,
      " [margin=", sprintf("%.3f", top5$median_EC_margin),
      ", ECwinner=", sprintf("%.2f", top5$frac_EC_relative_winner),
      "]"
    ),
    collapse = "; "
  ),
  "",
  "NEXT",
  "Upload GSE218166_step2D_REVIEW_SMALL.tar.gz.",
  "Then perform the FINAL author-aligned EC freeze and only afterward start Normal-Early-Advanced validation."
)

writeLines(
  summary_lines,
  file.path(OUT, "Step2D_DECISION_SUMMARY.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(LOG, "sessionInfo.txt")
)

cat(paste(summary_lines, collapse = "\n"))
cat("\nCompleted:", format(Sys.time()), "\n")
