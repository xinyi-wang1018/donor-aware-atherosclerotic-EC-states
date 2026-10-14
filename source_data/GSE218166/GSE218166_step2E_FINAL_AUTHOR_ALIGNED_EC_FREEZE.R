#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(20260813)

BASE <- "/root/GSE218166"

IN_RDS <- file.path(
  BASE,
  "step2D_AUTHOR_ALIGNED_EC_AUDIT",
  "objects",
  "GSE218166_step2D_allcells_author_aligned_lineage_scored.rds"
)

OUT  <- file.path(BASE, "step2E_FINAL_AUTHOR_ALIGNED_EC_FREEZE")
TAB  <- file.path(OUT, "tables")
PLOT <- file.path(OUT, "plots")
OBJ  <- file.path(OUT, "objects")
LOG  <- file.path(OUT, "logs")

for (d in c(OUT, TAB, PLOT, OBJ, LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

sink(file.path(LOG, "step2E_console.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

cat("============================================================\n")
cat("GSE218166 STEP 2E - FINAL AUTHOR-ALIGNED EC FREEZE\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

if (!file.exists(IN_RDS)) stop("Input not found: ", IN_RDS)

obj <- readRDS(IN_RDS)
DefaultAssay(obj) <- "RNA"

required_md <- c(
  "author_aligned_EC_tier",
  "EC_canonical_z",
  "EC_author_anchor_z",
  "EC_margin",
  "VSMC_contractile_z",
  "Fibroblast_z",
  "Myeloid_z",
  "Lymphoid_z",
  "tech",
  "seurat_clusters",
  "donor_id",
  "disease_stage",
  "coronary_position",
  "lesion_process",
  "plate"
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

# ----------------------------------------------------------------------
# FINAL PRE-SPECIFIED FREEZE
#
# MAIN:
# A_relative_strict_EC
#   EC_author_anchor_z > 0
#   EC_canonical_z > 0.25
#   EC_margin > 0.25 over ALL major non-EC programs
#
# EXPANDED:
# A + B_relative_EC
#
# REMODELING C:
# exploratory only; not used as primary EC identity set.
#
# Disease stage is NEVER used here.
# ----------------------------------------------------------------------

main_cells <- colnames(obj)[
  obj$author_aligned_EC_tier == "A_relative_strict_EC"
]

expanded_cells <- colnames(obj)[
  obj$author_aligned_EC_tier %in%
    c("A_relative_strict_EC", "B_relative_EC")
]

remodel_cells <- colnames(obj)[
  obj$author_aligned_EC_tier == "C_remodeling_EC_candidate"
]

if (length(main_cells) == 0) stop("No main strict EC cells found.")

obj$final_EC_freeze_set <- "Excluded_from_final_EC"
obj$final_EC_freeze_set[colnames(obj) %in% remodel_cells] <-
  "Exploratory_remodeling_candidate"
obj$final_EC_freeze_set[colnames(obj) %in% expanded_cells] <-
  "Expanded_AB_EC"
obj$final_EC_freeze_set[colnames(obj) %in% main_cells] <-
  "Main_A_strict_EC"

obj$final_EC_freeze_set <- factor(
  obj$final_EC_freeze_set,
  levels = c(
    "Excluded_from_final_EC",
    "Exploratory_remodeling_candidate",
    "Expanded_AB_EC",
    "Main_A_strict_EC"
  )
)

main_obj <- subset(obj, cells = main_cells)
expanded_obj <- subset(obj, cells = expanded_cells)
remodel_obj <- subset(obj, cells = remodel_cells)

main_obj$final_analysis_role <- "PRIMARY_EC"
expanded_obj$final_analysis_role <- "EXPANDED_SENSITIVITY_EC"
remodel_obj$final_analysis_role <- "EXPLORATORY_REMODELING_ONLY"

# donor x stage unit
for (xname in c("main_obj", "expanded_obj", "remodel_obj")) {
  x <- get(xname)
  x$donor_stage_unit <- paste(
    x$donor_id,
    as.character(x$disease_stage),
    sep = "__"
  )
  assign(xname, x)
}

# ----------------------------------------------------------------------
# 1. Save final objects
# ----------------------------------------------------------------------
saveRDS(
  main_obj,
  file.path(OBJ, "GSE218166_FINAL_MAIN_A_STRICT_EC.rds"),
  compress = FALSE
)

saveRDS(
  expanded_obj,
  file.path(OBJ, "GSE218166_FINAL_EXPANDED_AB_EC.rds"),
  compress = FALSE
)

saveRDS(
  remodel_obj,
  file.path(OBJ, "GSE218166_EXPLORATORY_REMODELING_EC_CANDIDATES.rds"),
  compress = FALSE
)

# ----------------------------------------------------------------------
# 2. Freeze summary
# ----------------------------------------------------------------------
freeze_summary <- data.frame(
  analysis_set = c(
    "MAIN_A_STRICT_EC",
    "EXPANDED_AB_EC",
    "EXPLORATORY_REMODELING_CANDIDATE"
  ),
  n_cells = c(
    ncol(main_obj),
    ncol(expanded_obj),
    ncol(remodel_obj)
  ),
  definition = c(
    "A_relative_strict_EC: author EC anchor >0; canonical EC >0.25; EC margin >0.25",
    "A_relative_strict_EC + B_relative_EC",
    "C_remodeling_EC_candidate only; NOT primary identity set"
  ),
  downstream_role = c(
    "PRIMARY Normal-Early-Advanced validation",
    "Sensitivity analysis",
    "Exploratory state mapping only"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  freeze_summary,
  file.path(TAB, "01_final_EC_freeze_summary.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 3. Stage and donor counts
# ----------------------------------------------------------------------
summarise_set <- function(sobj, set_name) {
  sobj@meta.data %>%
    count(
      tech,
      disease_stage,
      donor_id,
      donor_stage_unit,
      name = "n_cells"
    ) %>%
    mutate(analysis_set = set_name)
}

donor_stage_counts <- bind_rows(
  summarise_set(main_obj, "MAIN_A_STRICT_EC"),
  summarise_set(expanded_obj, "EXPANDED_AB_EC"),
  summarise_set(remodel_obj, "EXPLORATORY_REMODELING_CANDIDATE")
)

write.csv(
  donor_stage_counts,
  file.path(TAB, "02_EC_counts_by_set_tech_donor_stage.csv"),
  row.names = FALSE
)

stage_summary <- donor_stage_counts %>%
  group_by(analysis_set, tech, disease_stage) %>%
  summarise(
    total_cells = sum(n_cells),
    n_donors = n_distinct(donor_id),
    n_donor_stage_units = n(),
    median_cells_per_donor_stage = median(n_cells),
    min_cells_per_donor_stage = min(n_cells),
    max_cells_per_donor_stage = max(n_cells),
    n_units_ge10 = sum(n_cells >= 10),
    n_units_ge20 = sum(n_cells >= 20),
    n_units_ge30 = sum(n_cells >= 30),
    .groups = "drop"
  )

write.csv(
  stage_summary,
  file.path(TAB, "03_EC_stage_replicate_adequacy.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 4. Broad-cluster provenance
# ----------------------------------------------------------------------
cluster_origin <- bind_rows(
  main_obj@meta.data %>%
    count(tech, seurat_clusters, name = "n_cells") %>%
    mutate(analysis_set = "MAIN_A_STRICT_EC"),
  expanded_obj@meta.data %>%
    count(tech, seurat_clusters, name = "n_cells") %>%
    mutate(analysis_set = "EXPANDED_AB_EC"),
  remodel_obj@meta.data %>%
    count(tech, seurat_clusters, name = "n_cells") %>%
    mutate(analysis_set = "EXPLORATORY_REMODELING_CANDIDATE")
) %>%
  group_by(analysis_set, tech) %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(
  cluster_origin,
  file.path(TAB, "04_final_EC_broad_cluster_provenance.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 5. Relative lineage purity summary
# ----------------------------------------------------------------------
score_cols <- c(
  "EC_author_anchor_z",
  "EC_canonical_z",
  "EC_margin",
  "VSMC_contractile_z",
  "Fibroblast_z",
  "Myeloid_z",
  "Lymphoid_z",
  "EC1_dflow_candidate_z",
  "EC2_ACKR1_candidate_z",
  "EC4_EndMT_candidate_z",
  "EC5_SLCO4A1_candidate_z",
  "Project_early_z",
  "Project_late_z"
)

score_summary_one <- function(sobj, set_name) {
  md <- sobj@meta.data
  bind_rows(lapply(score_cols, function(cc) {
    data.frame(
      analysis_set = set_name,
      score = cc,
      median = median(md[[cc]], na.rm = TRUE),
      mean = mean(md[[cc]], na.rm = TRUE),
      q25 = quantile(md[[cc]], 0.25, na.rm = TRUE),
      q75 = quantile(md[[cc]], 0.75, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

score_summary <- bind_rows(
  score_summary_one(main_obj, "MAIN_A_STRICT_EC"),
  score_summary_one(expanded_obj, "EXPANDED_AB_EC"),
  score_summary_one(remodel_obj, "EXPLORATORY_REMODELING_CANDIDATE")
)

write.csv(
  score_summary,
  file.path(TAB, "05_final_EC_relative_lineage_score_summary.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 6. Main-set stage diagnostics ONLY
# No inferential stage testing yet.
# ----------------------------------------------------------------------
main_stage_scores <- main_obj@meta.data %>%
  group_by(tech, disease_stage) %>%
  summarise(
    n_cells = n(),
    n_donors = n_distinct(donor_id),
    median_EC1 = median(EC1_dflow_candidate_z, na.rm = TRUE),
    median_EC2_ACKR1 = median(EC2_ACKR1_candidate_z, na.rm = TRUE),
    median_EC4_EndMT = median(EC4_EndMT_candidate_z, na.rm = TRUE),
    median_EC5 = median(EC5_SLCO4A1_candidate_z, na.rm = TRUE),
    median_project_early = median(Project_early_z, na.rm = TRUE),
    median_project_late = median(Project_late_z, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  main_stage_scores,
  file.path(
    TAB,
    "06_MAIN_A_stage_score_descriptives_DIAGNOSTIC_ONLY.csv"
  ),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 7. Step3 recommendation
# ----------------------------------------------------------------------
main_adequacy <- stage_summary %>%
  filter(analysis_set == "MAIN_A_STRICT_EC")

fl_rows <- main_adequacy %>% filter(tech == "FL")
de_rows <- main_adequacy %>% filter(tech == "3DE")

fl_primary_ok <- (
  nrow(fl_rows) == 3 &&
  all(fl_rows$n_donors >= 3) &&
  all(fl_rows$total_cells >= 50)
)

de_three_stage_ok <- (
  nrow(de_rows) == 3 &&
  all(de_rows$n_donors >= 3) &&
  all(de_rows$total_cells >= 50)
)

recommendation <- data.frame(
  component = c(
    "Primary stage validation",
    "Expanded sensitivity",
    "3DE technology validation",
    "Exploratory remodeling candidates"
  ),
  dataset_or_set = c(
    "FL + MAIN_A_STRICT_EC",
    "FL + EXPANDED_AB_EC",
    "3DE + MAIN_A_STRICT_EC",
    "C_remodeling candidates"
  ),
  recommendation = c(
    ifelse(
      fl_primary_ok,
      "PROCEED: donor x disease-stage pseudobulk",
      "REVIEW BEFORE INFERENCE"
    ),
    "Run same endpoints as sensitivity analysis",
    ifelse(
      de_three_stage_ok,
      "Three-stage sensitivity feasible",
      "Do not treat as equal-strength three-stage validation; Early is underpowered"
    ),
    "Use for state interpretation only, not primary EC identity"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  recommendation,
  file.path(TAB, "07_STEP3_RECOMMENDATION.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 8. Figures - pure white background
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
  p1 <- DimPlot(
    obj,
    reduction = reduction_use,
    group.by = "final_EC_freeze_set",
    raster = TRUE,
    pt.size = 0.25
  ) +
    ggtitle("Final author-aligned EC freeze") +
    theme_classic(base_size = 10)

  save_white(
    p1,
    "01_FINAL_author_aligned_EC_freeze_UMAP.png",
    11,
    8
  )
}

p2 <- ggplot(
  stage_summary %>%
    filter(
      analysis_set %in%
        c("MAIN_A_STRICT_EC", "EXPANDED_AB_EC")
    ),
  aes(
    x = disease_stage,
    y = total_cells,
    fill = analysis_set
  )
) +
  geom_col(position = "dodge") +
  facet_wrap(~tech) +
  theme_classic(base_size = 10) +
  labs(
    title = "Final EC cell counts by stage and technology",
    x = NULL,
    y = "Number of EC cells"
  )

save_white(
  p2,
  "02_FINAL_EC_counts_by_stage_tech.png",
  11,
  6.5
)

p3 <- ggplot(
  donor_stage_counts %>%
    filter(analysis_set == "MAIN_A_STRICT_EC"),
  aes(
    x = disease_stage,
    y = n_cells,
    shape = tech
  )
) +
  geom_jitter(width = 0.15, height = 0, size = 2.5) +
  theme_classic(base_size = 10) +
  labs(
    title = "Main strict EC cells per donor-stage unit",
    x = NULL,
    y = "EC cells"
  )

save_white(
  p3,
  "03_MAIN_A_EC_cells_per_donor_stage.png",
  9,
  6
)

# ----------------------------------------------------------------------
# 9. Final summary
# ----------------------------------------------------------------------
summary_lines <- c(
  "============================================================",
  "GSE218166 STEP2E FINAL AUTHOR-ALIGNED EC FREEZE",
  "============================================================",
  paste0("Main A strict EC cells: ", ncol(main_obj)),
  paste0("Expanded A+B EC cells: ", ncol(expanded_obj)),
  paste0("Exploratory remodeling-C cells: ", ncol(remodel_obj)),
  "",
  "FINAL PRIMARY EC DEFINITION",
  "A_relative_strict_EC only.",
  "This requires positive author EC anchor, strong canonical EC score,",
  "and EC program dominance by >0.25 over the strongest competing non-EC program.",
  "",
  "IMPORTANT",
  "Disease stage was NOT used to define the final EC population.",
  "Broad cluster membership alone was NOT used as the final freeze rule.",
  "Cluster 18 is therefore only partly retained: EC-dominant cells are kept,",
  "while cells without sufficient EC dominance are excluded.",
  "",
  paste0(
    "FL primary three-stage adequacy: ",
    ifelse(fl_primary_ok, "PASS", "REVIEW")
  ),
  paste0(
    "3DE equal-strength three-stage adequacy: ",
    ifelse(de_three_stage_ok, "PASS", "FAIL/WEAK")
  ),
  "",
  "NEXT",
  "Proceed to Step3 Normal-Early-Advanced validation using FL Main A strict EC.",
  "Repeat key endpoints in Expanded A+B EC as sensitivity.",
  "Use 3DE only as limited technology sensitivity if its Early group remains sparse."
)

writeLines(
  summary_lines,
  file.path(OUT, "Step2E_FINAL_DECISION_SUMMARY.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(LOG, "sessionInfo.txt")
)

cat(paste(summary_lines, collapse = "\n"))
cat("\nCompleted:", format(Sys.time()), "\n")
