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

IN_RDS <- file.path(
  BASE,
  "step2B_V3_EC_recluster_lineage_audit",
  "objects",
  "GSE218166_step2B_V3_candidate_EC_reclustered_review.rds"
)

OUT  <- file.path(BASE, "step2C_FINAL_EC_FREEZE")
TAB  <- file.path(OUT, "tables")
PLOT <- file.path(OUT, "plots")
OBJ  <- file.path(OUT, "objects")
LOG  <- file.path(OUT, "logs")

for (d in c(OUT, TAB, PLOT, OBJ, LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

sink(file.path(LOG, "step2C_console.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

cat("============================================================\n")
cat("GSE218166 STEP 2C - FINAL STRICT EC FREEZE\n")
cat("Start:", format(Sys.time()), "\n")
cat("Input:", IN_RDS, "\n")
cat("============================================================\n\n")

if (!file.exists(IN_RDS)) {
  stop("Input Step2B V3 RDS not found: ", IN_RDS)
}

obj <- readRDS(IN_RDS)
DefaultAssay(obj) <- "RNA"

required_md <- c(
  "EC_recluster",
  "inclusive_EC_support",
  "canonical_EC_support",
  "hard_nonEC_conflict",
  "donor_id",
  "disease_stage",
  "tech",
  "plate",
  "coronary_position",
  "lesion_process",
  "seurat_clusters",
  "prefilter_origin",
  "Individual_Level_reconstructed"
)

missing_md <- setdiff(required_md, colnames(obj@meta.data))
if (length(missing_md) > 0) {
  stop("Missing metadata: ", paste(missing_md, collapse = ", "))
}

obj$EC_recluster <- as.character(obj$EC_recluster)
obj$disease_stage <- factor(
  as.character(obj$disease_stage),
  levels = c("Normal", "Early", "Advanced"),
  ordered = TRUE
)

# -------------------------------------------------------------------------
# 1. Manual cluster decision from Step2B marker review
# -------------------------------------------------------------------------
manual_review <- data.frame(
  EC_recluster = as.character(0:8),
  manual_identity = c(
    "Mural_SMC_pericyte",
    "Fibroblast_stromal",
    "Endothelial_high_confidence",
    "Borderline_mixed_EC_like",
    "Activated_stromal_immune_like",
    "T_NK_immune",
    "Schwann_perineurial_like",
    "Mural_pericyte_epithelial_like",
    "Mast_cell"
  ),
  main_freeze_decision = c(
    "Exclude",
    "Exclude",
    "Include",
    "Exclude_from_main_keep_for_sensitivity",
    "Exclude",
    "Exclude",
    "Exclude",
    "Exclude",
    "Exclude"
  ),
  key_reason = c(
    "RGS5/MYH11/ACTA2/MYOCD-dominant mural program",
    "DCN/LUM/DPT/PI16/MFAP5 fibroblast program",
    "ACKR1/ADGRL4/EMCN/F8/TEK/FLT1/PECAM1/CDH5/PTPRB endothelial program",
    "EC markers present but atypical/mixed technical and lineage features; retain only as expanded sensitivity",
    "CCL19/CXCL1/IL6/PDPN/CD248 plus immunoglobulin/stromal signals",
    "CD3D/CD3G/TRBC2/PTPRC/NKG7-related immune program",
    "PLP1/PMP2/S100B/CHL1 neural-Schwann/perineurial program",
    "ACTA2/TAGLN/RGS5/MYH11 plus epithelial-like signals",
    "TPSAB1/TPSB2/MS4A2/HDC/CPA3/KIT mast-cell program"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  manual_review,
  file.path(TAB, "01_manual_cluster_identity_and_freeze_decision.csv"),
  row.names = FALSE
)

identity_map <- setNames(
  manual_review$manual_identity,
  manual_review$EC_recluster
)

obj$manual_identity <- unname(
  identity_map[as.character(obj$EC_recluster)]
)

# -------------------------------------------------------------------------
# 2. FINAL FREEZE SETS
#
# MAIN STRICT:
#   cluster 2 + inclusive EC marker evidence + no hard non-EC conflict
#
# ULTRA STRICT SENSITIVITY:
#   cluster 2 + canonical EC support
#
# EXPANDED SENSITIVITY:
#   clusters 2 or 3 + inclusive EC support
#
# Disease stage is NOT used in any freeze rule.
# -------------------------------------------------------------------------
main_cells <- colnames(obj)[
  obj$EC_recluster == "2" &
    obj$inclusive_EC_support &
    !obj$hard_nonEC_conflict
]

ultra_cells <- colnames(obj)[
  obj$EC_recluster == "2" &
    obj$canonical_EC_support &
    !obj$hard_nonEC_conflict
]

expanded_cells <- colnames(obj)[
  obj$EC_recluster %in% c("2", "3") &
    obj$inclusive_EC_support &
    !obj$hard_nonEC_conflict
]

obj$final_freeze_label <- "Excluded_from_main"
obj$final_freeze_label[colnames(obj) %in% expanded_cells] <-
  "Expanded_sensitivity_EC"
obj$final_freeze_label[colnames(obj) %in% main_cells] <-
  "Main_strict_EC"
obj$final_freeze_label[colnames(obj) %in% ultra_cells] <-
  "Ultra_strict_EC"

obj$final_freeze_label <- factor(
  obj$final_freeze_label,
  levels = c(
    "Excluded_from_main",
    "Expanded_sensitivity_EC",
    "Main_strict_EC",
    "Ultra_strict_EC"
  )
)

freeze_summary <- data.frame(
  analysis_set = c(
    "Main_strict_EC",
    "Ultra_strict_EC",
    "Expanded_sensitivity_EC"
  ),
  definition = c(
    "EC_recluster==2 AND inclusive_EC_support AND no hard nonEC conflict",
    "EC_recluster==2 AND canonical_EC_support AND no hard nonEC conflict",
    "EC_recluster in {2,3} AND inclusive_EC_support AND no hard nonEC conflict"
  ),
  n_cells = c(
    length(main_cells),
    length(ultra_cells),
    length(expanded_cells)
  ),
  role = c(
    "PRIMARY downstream Normal-Early-Advanced validation",
    "Sensitivity analysis for maximal identity purity",
    "Sensitivity analysis testing whether borderline EC-like cluster 3 changes conclusions"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  freeze_summary,
  file.path(TAB, "02_final_freeze_set_summary.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 3. Create frozen objects
# -------------------------------------------------------------------------
main_obj <- subset(obj, cells = main_cells)
ultra_obj <- subset(obj, cells = ultra_cells)
expanded_obj <- subset(obj, cells = expanded_cells)

main_obj$freeze_set <- "Main_strict_EC"
ultra_obj$freeze_set <- "Ultra_strict_EC"
expanded_obj$freeze_set <- "Expanded_sensitivity_EC"

saveRDS(
  main_obj,
  file.path(OBJ, "GSE218166_FINAL_main_strict_EC.rds"),
  compress = FALSE
)

saveRDS(
  ultra_obj,
  file.path(OBJ, "GSE218166_FINAL_ultra_strict_EC_sensitivity.rds"),
  compress = FALSE
)

saveRDS(
  expanded_obj,
  file.path(OBJ, "GSE218166_FINAL_expanded_EC_sensitivity.rds"),
  compress = FALSE
)

# -------------------------------------------------------------------------
# 4. Main strict EC sample / stage audit
# -------------------------------------------------------------------------
main_meta <- main_obj@meta.data

stage_tech_counts <- main_meta %>%
  count(tech, disease_stage, name = "n_cells") %>%
  arrange(tech, disease_stage)

write.csv(
  stage_tech_counts,
  file.path(TAB, "03_main_strict_EC_cells_by_tech_stage.csv"),
  row.names = FALSE
)

donor_stage_counts <- main_meta %>%
  count(
    tech,
    donor_id,
    disease_stage,
    Individual_Level_reconstructed,
    name = "n_cells"
  ) %>%
  arrange(tech, disease_stage, donor_id)

write.csv(
  donor_stage_counts,
  file.path(TAB, "04_main_strict_EC_cells_by_donor_stage.csv"),
  row.names = FALSE
)

replicate_summary <- donor_stage_counts %>%
  group_by(tech, disease_stage) %>%
  summarise(
    n_donor_stage_units = n(),
    n_unique_donors = n_distinct(donor_id),
    total_cells = sum(n_cells),
    median_cells_per_unit = median(n_cells),
    min_cells_per_unit = min(n_cells),
    max_cells_per_unit = max(n_cells),
    n_units_ge_10_cells = sum(n_cells >= 10),
    n_units_ge_20_cells = sum(n_cells >= 20),
    n_units_ge_30_cells = sum(n_cells >= 30),
    .groups = "drop"
  ) %>%
  mutate(
    three_stage_primary_eligibility =
      ifelse(
        n_unique_donors >= 3 &
          total_cells >= 50 &
          n_units_ge_10_cells >= 3,
        "PASS",
        "WEAK"
      )
  )

write.csv(
  replicate_summary,
  file.path(TAB, "05_main_strict_EC_replicate_adequacy_by_tech_stage.csv"),
  row.names = FALSE
)

# repeated donors across stages
donor_stage_matrix <- donor_stage_counts %>%
  mutate(present = 1L) %>%
  select(tech, donor_id, disease_stage, present) %>%
  distinct() %>%
  pivot_wider(
    names_from = disease_stage,
    values_from = present,
    values_fill = 0
  )

write.csv(
  donor_stage_matrix,
  file.path(TAB, "06_main_strict_EC_donor_stage_presence_matrix.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 5. Provenance audit: where frozen ECs came from
# -------------------------------------------------------------------------
origin_summary <- main_meta %>%
  count(prefilter_origin, name = "n_cells") %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  arrange(desc(n_cells))

write.csv(
  origin_summary,
  file.path(TAB, "07_main_strict_EC_prefilter_origin.csv"),
  row.names = FALSE
)

broad_cluster_summary <- main_meta %>%
  count(seurat_clusters, name = "n_cells") %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  arrange(desc(n_cells))

write.csv(
  broad_cluster_summary,
  file.path(TAB, "08_main_strict_EC_original_broad_cluster_origin.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 6. Marker purity audit in final sets
# -------------------------------------------------------------------------
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

main_obj <- NormalizeData(main_obj, verbose = FALSE)
ultra_obj <- NormalizeData(ultra_obj, verbose = FALSE)
expanded_obj <- NormalizeData(expanded_obj, verbose = FALSE)

marker_groups <- list(
  EC_core = c(
    "TIE1","PECAM1","VWF","KDR","CDH5","CLDN5","ERG","ESAM",
    "EMCN","FLT1","TEK","PTPRB","ADGRL4"
  ),
  Immune = c(
    "PTPRC","LST1","TYROBP","FCER1G","LYZ","C5AR1"
  ),
  T_NK = c(
    "CD3D","CD3E","TRBC2","NKG7","GNLY"
  ),
  Mural = c(
    "ACTA2","TAGLN","MYH11","RGS5","PDGFRB"
  ),
  Fibroblast = c(
    "COL1A1","COL1A2","COL3A1","DCN","LUM"
  ),
  Platelet_RBC = c(
    "PPBP","PF4","HBB","HBA1","HBA2"
  ),
  Arterial = c(
    "CXCL12","EFNB2","GJA4","GJA5","SOX17","HEY1","SEMA3G","BMX"
  ),
  Venous_microvascular = c(
    "ACKR1","PLVAP","COL15A1","CA4","RGCC"
  ),
  Author_EC5 = c(
    "SLCO4A1","PRDM15","PROM1"
  ),
  Project_early = c(
    "ANGPT2","CCL14","PGF","LOX","CTHRC1","THBS1"
  ),
  Project_late = c(
    "SULF1","ITLN1","MGP","ELN","COL8A1","FN1","SPARC","COL6A2","IGFBP5"
  )
)

marker_summary_one <- function(sobj, set_name) {
  mat <- get_layer(sobj, "data")
  if (is.null(mat)) stop("No RNA data layer for ", set_name)

  out <- lapply(names(marker_groups), function(grp) {
    gg <- intersect(marker_groups[[grp]], rownames(mat))
    if (length(gg) == 0) return(NULL)

    data.frame(
      analysis_set = set_name,
      marker_group = grp,
      gene = gg,
      mean_expression = as.numeric(Matrix::rowMeans(mat[gg, , drop = FALSE])),
      percent_expressing = as.numeric(
        Matrix::rowMeans(mat[gg, , drop = FALSE] > 0) * 100
      ),
      stringsAsFactors = FALSE
    )
  })

  bind_rows(out)
}

marker_summary <- bind_rows(
  marker_summary_one(main_obj, "Main_strict_EC"),
  marker_summary_one(ultra_obj, "Ultra_strict_EC"),
  marker_summary_one(expanded_obj, "Expanded_sensitivity_EC")
)

write.csv(
  marker_summary,
  file.path(TAB, "09_final_sets_marker_expression_summary.csv"),
  row.names = FALSE
)

marker_group_summary <- marker_summary %>%
  group_by(analysis_set, marker_group) %>%
  summarise(
    mean_of_mean_expression = mean(mean_expression),
    mean_percent_expressing = mean(percent_expressing),
    .groups = "drop"
  )

write.csv(
  marker_group_summary,
  file.path(TAB, "10_final_sets_marker_group_summary.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 7. FL vs 3DE downstream recommendation
# -------------------------------------------------------------------------
eligibility_wide <- replicate_summary %>%
  select(
    tech, disease_stage,
    n_unique_donors, total_cells,
    n_units_ge_10_cells,
    three_stage_primary_eligibility
  )

fl_pass <- eligibility_wide %>%
  filter(tech == "FL") %>%
  summarise(
    all_three_stages_pass =
      all(three_stage_primary_eligibility == "PASS")
  ) %>%
  pull(all_three_stages_pass)

de3_pass <- eligibility_wide %>%
  filter(tech == "3DE") %>%
  summarise(
    all_three_stages_pass =
      all(three_stage_primary_eligibility == "PASS")
  ) %>%
  pull(all_three_stages_pass)

recommendation <- data.frame(
  analysis = c(
    "Primary three-stage validation",
    "Independent technology sensitivity",
    "Cell-level descriptive visualization"
  ),
  recommended_dataset = c(
    ifelse(fl_pass, "FL", "REVIEW"),
    "3DE",
    "FL + 3DE shown separately"
  ),
  rule = c(
    "Use donor x disease-stage pseudobulk; account for repeated donors across stages",
    ifelse(
      de3_pass,
      "Three-stage sensitivity feasible",
      "Do NOT use as equal-strength three-stage test if Early replication/cell count is sparse; use limited sensitivity only"
    ),
    "Never treat cells as biological replicates"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  recommendation,
  file.path(TAB, "11_step3_analysis_recommendation.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 8. Figures - pure white background
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

# Original Step2B UMAP showing manual identities
p1 <- DimPlot(
  obj,
  reduction = "umap.ec",
  group.by = "manual_identity",
  raster = TRUE,
  pt.size = 0.30,
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Manual lineage review of Step2B candidate clusters") +
  theme_classic(base_size = 10)

save_white(p1, "01_manual_cluster_identity_UMAP.png", 12, 8)

# Main strict location
obj$main_strict_flag <- ifelse(
  colnames(obj) %in% main_cells,
  "Main_strict_EC",
  "Not_main"
)

p2 <- DimPlot(
  obj,
  reduction = "umap.ec",
  group.by = "main_strict_flag",
  raster = TRUE,
  pt.size = 0.30
) +
  ggtitle("Final main strict EC freeze") +
  theme_classic(base_size = 11)

save_white(p2, "02_FINAL_main_strict_EC_UMAP.png", 10, 8)

# Stage counts in main strict set
p3 <- ggplot(
  stage_tech_counts,
  aes(x = disease_stage, y = n_cells, fill = tech)
) +
  geom_col(position = "dodge") +
  theme_classic(base_size = 11) +
  labs(
    title = "Main strict EC counts by technology and disease stage",
    x = NULL,
    y = "Number of frozen ECs"
  )

save_white(p3, "03_main_strict_EC_counts_by_stage_tech.png", 9, 6)

# donor-stage cell counts
p4 <- ggplot(
  donor_stage_counts,
  aes(
    x = disease_stage,
    y = n_cells,
    shape = tech
  )
) +
  geom_jitter(width = 0.15, height = 0, size = 2.5) +
  theme_classic(base_size = 11) +
  labs(
    title = "Frozen EC counts per donor-stage unit",
    x = NULL,
    y = "EC cells per donor-stage unit"
  )

save_white(p4, "04_main_strict_EC_donor_stage_counts.png", 9, 6)

# marker dot plot across manual cluster identities
dot_genes <- unique(c(
  "TIE1","PECAM1","VWF","KDR","CDH5","CLDN5","ERG","ESAM","EMCN",
  "ACKR1","PLVAP","TEK","PTPRB",
  "PTPRC","LST1","TYROBP",
  "CD3D","NKG7",
  "ACTA2","TAGLN","MYH11","RGS5","PDGFRB",
  "COL1A1","COL3A1","DCN","LUM",
  "PPBP","PF4","HBB",
  "CXCL12","EFNB2","GJA4","SOX17","HEY1",
  "SLCO4A1","PRDM15","PROM1",
  "ANGPT2","PGF","LOX","CTHRC1","THBS1",
  "SULF1","ITLN1","MGP","ELN","COL8A1","FN1"
))
dot_genes <- intersect(dot_genes, rownames(obj))

if (length(dot_genes) > 0) {
  p5 <- DotPlot(
    obj,
    features = dot_genes,
    group.by = "manual_identity",
    assay = "RNA"
  ) +
    RotatedAxis() +
    theme_classic(base_size = 7.5) +
    ggtitle("Manual lineage decision marker audit")

  save_white(p5, "05_manual_identity_marker_DotPlot.png", 20, 8)
}

# -------------------------------------------------------------------------
# 9. Human-readable decision summary
# -------------------------------------------------------------------------
fl_stage <- replicate_summary %>% filter(tech == "FL")
de_stage <- replicate_summary %>% filter(tech == "3DE")

summary_lines <- c(
  "============================================================",
  "GSE218166 STEP 2C - FINAL STRICT EC FREEZE",
  "============================================================",
  paste0("Run time: ", format(Sys.time())),
  "",
  "FINAL DECISION",
  "Main EC identity cluster = Step2B EC_recluster 2.",
  "Borderline EC-like cluster 3 is excluded from the primary freeze and retained only for sensitivity analysis.",
  "Clusters 0,1,4,5,6,7,8 are excluded as non-EC/mixed non-EC populations based on marker programs.",
  "",
  paste0("Main strict EC cells: ", length(main_cells)),
  paste0("Ultra-strict EC sensitivity cells: ", length(ultra_cells)),
  paste0("Expanded cluster2+3 EC sensitivity cells: ", length(expanded_cells)),
  "",
  "IMPORTANT",
  "Disease stage was NOT used in the freeze rule.",
  "The automatic Step2B HIGH_CONFIDENCE label was NOT used as final truth.",
  "Manual lineage marker review overrides the automatic review flag.",
  "",
  "STEP 3 DESIGN",
  paste0(
    "FL three-stage eligibility: ",
    ifelse(fl_pass, "PASS", "REVIEW/FAIL")
  ),
  paste0(
    "3DE three-stage eligibility: ",
    ifelse(de3_pass, "PASS", "WEAK/FAIL")
  ),
  "Primary stage testing should use donor x disease-stage pseudobulk rather than cell-level P values.",
  "Repeated donors appearing in more than one stage must be accounted for in the statistical design.",
  "3DE should not be given equal weight to FL if its Early group remains sparse after strict EC freeze.",
  ""
)

writeLines(
  summary_lines,
  file.path(OUT, "Step2C_FINAL_DECISION_SUMMARY.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(LOG, "sessionInfo.txt")
)

cat(paste(summary_lines, collapse = "\n"))
cat("\nCompleted:", format(Sys.time()), "\n")
