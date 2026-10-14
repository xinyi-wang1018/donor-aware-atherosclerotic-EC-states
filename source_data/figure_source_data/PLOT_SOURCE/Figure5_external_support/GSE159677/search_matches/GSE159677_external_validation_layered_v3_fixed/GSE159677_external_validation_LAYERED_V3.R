options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 40 * 1024^3)
set.seed(20260731)

BASE <- "/root/GSE159677_external_validation_layered_v3_fixed"
OUT_DIR <- file.path(BASE, "results")
OBJ_DIR <- file.path(BASE, "objects")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)

path_lines <- readLines(file.path(BASE, "input_paths.txt"))
path_map <- setNames(sub("^[^=]+=" , "", path_lines), sub("=.*$", "", path_lines))
FULL_RDS <- unname(path_map[["FULL_RDS"]])
REF_RDS <- unname(path_map[["REF_RDS"]])

required <- c(
  "Seurat", "SeuratObject", "Matrix", "ggplot2", "dplyr",
  "tidyr", "patchwork", "pheatmap", "harmony"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required R packages: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(pheatmap)
})

white_theme <- theme_classic(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )
theme_set(white_theme)

save_png <- function(plot, filename, width = 10, height = 7, dpi = 300) {
  ggsave(
    filename = file.path(OUT_DIR, filename), plot = plot,
    width = width, height = height, dpi = dpi,
    bg = "white", limitsize = FALSE
  )
}

join_layers_safe <- function(obj, assay = "RNA") {
  tryCatch(
    SeuratObject::JoinLayers(obj, assay = assay),
    error = function(e1) tryCatch(
      Seurat::JoinLayers(obj, assay = assay),
      error = function(e2) obj
    )
  )
}

get_assay_matrix <- function(obj, assay = "RNA", layer = "data") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
}

run_harmony_strict <- function(obj, group_var, reduction_name = "pca", save_name = "harmony") {
  message("Running Harmony by ", group_var, " ...")
  out <- tryCatch(
    harmony::RunHarmony(
      object = obj,
      group.by.vars = group_var,
      reduction.use = reduction_name,
      reduction.save = save_name,
      verbose = TRUE
    ),
    error = function(e1) {
      message("First Harmony call failed: ", conditionMessage(e1))
      tryCatch(
        harmony::RunHarmony(
          object = obj,
          group.by.vars = group_var,
          reduction.use = reduction_name,
          verbose = TRUE
        ),
        error = function(e2) {
          stop(
            "Harmony failed.\nFirst error: ", conditionMessage(e1),
            "\nSecond error: ", conditionMessage(e2)
          )
        }
      )
    }
  )
  if (!save_name %in% Reductions(out) && !"harmony" %in% Reductions(out)) {
    stop("Harmony returned without a harmony reduction.")
  }
  out
}

marker_count <- function(count_mat, genes) {
  genes <- intersect(genes, rownames(count_mat))
  if (length(genes) == 0) return(rep(0L, ncol(count_mat)))
  as.integer(Matrix::colSums(count_mat[genes, , drop = FALSE] > 0))
}

gene_z_score <- function(data_mat, genes, cap = 5) {
  genes <- intersect(genes, rownames(data_mat))
  if (length(genes) == 0) return(rep(NA_real_, ncol(data_mat)))
  x <- as.matrix(data_mat[genes, , drop = FALSE])
  mu <- rowMeans(x)
  s <- apply(x, 1, stats::sd)
  s[!is.finite(s) | s < 1e-8] <- 1
  z <- sweep(x, 1, mu, "-")
  z <- sweep(z, 1, s, "/")
  z[z > cap] <- cap
  z[z < -cap] <- -cap
  colMeans(z)
}

mean_log_score <- function(data_mat, genes) {
  genes <- intersect(genes, rownames(data_mat))
  if (length(genes) == 0) return(rep(NA_real_, ncol(data_mat)))
  as.numeric(Matrix::colMeans(data_mat[genes, , drop = FALSE]))
}

safe_num <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  x
}

paired_stats <- function(df, value_col, feature_col, expected_map = NULL, eligible_patients = NULL) {
  features <- unique(df[[feature_col]])
  out <- lapply(features, function(feat) {
    dd <- df[df[[feature_col]] == feat, c("patient_id", "condition", value_col), drop = FALSE]
    if (!is.null(eligible_patients)) dd <- dd[dd$patient_id %in% eligible_patients, , drop = FALSE]
    wide <- tidyr::pivot_wider(dd, names_from = condition, values_from = all_of(value_col))
    if (!all(c("PA", "AC") %in% colnames(wide))) return(NULL)
    wide <- wide[complete.cases(wide[, c("PA", "AC")]), , drop = FALSE]
    if (nrow(wide) == 0) return(NULL)
    diff <- safe_num(wide$AC - wide$PA)
    diff <- diff[!is.na(diff)]
    if (length(diff) == 0) return(NULL)
    p <- if (length(diff) >= 2) {
      tryCatch(
        wilcox.test(wide$AC, wide$PA, paired = TRUE, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    } else NA_real_
    dz <- if (length(diff) >= 2 && is.finite(sd(diff)) && sd(diff) > 0) {
      mean(diff) / sd(diff)
    } else NA_real_
    expected <- if (!is.null(expected_map) && feat %in% names(expected_map)) expected_map[[feat]] else NA_character_
    n_expected <- if (is.na(expected)) {
      NA_integer_
    } else if (expected == "AC>PA") {
      sum(diff > 0)
    } else {
      sum(diff < 0)
    }
    get_diff <- function(pid) {
      idx <- match(pid, wide$patient_id)
      if (is.na(idx)) NA_real_ else wide$AC[idx] - wide$PA[idx]
    }
    data.frame(
      feature = feat,
      n_pairs = length(diff),
      mean_AC_minus_PA = mean(diff),
      median_AC_minus_PA = median(diff),
      paired_Cohen_dz = dz,
      wilcoxon_p = p,
      expected_direction = expected,
      n_pairs_in_expected_direction = n_expected,
      fraction_in_expected_direction = ifelse(is.na(n_expected), NA_real_, n_expected / length(diff)),
      P1_AC_minus_PA = get_diff("P1"),
      P2_AC_minus_PA = get_diff("P2"),
      P3_AC_minus_PA = get_diff("P3")
    )
  })
  bind_rows(out)
}

normalize_cluster_label <- function(x) {
  x <- as.character(x)
  x <- gsub("^cluster[_ ]*", "", x, ignore.case = TRUE)
  x <- gsub("^C", "", x, ignore.case = TRUE)
  x
}

find_reference_cluster_labels <- function(ref) {
  wanted <- c("1", "2", "3", "4", "5", "7")
  valid_cluster_vector <- function(v) {
    vv <- normalize_cluster_label(v)
    uu <- unique(vv[!is.na(vv) & vv != ""])
    length(uu) >= length(wanted) && length(uu) <= 40 &&
      all(grepl("^[0-9]+$", uu)) && all(wanted %in% uu)
  }
  idv <- normalize_cluster_label(Idents(ref))
  if (valid_cluster_vector(idv)) return(idv)
  preferred <- c(
    "seurat_clusters", "cluster", "clusters", "cluster_id",
    "EC_cluster", "ec_cluster", "sling_cluster", "slingshot_cluster"
  )
  candidates <- unique(c(preferred[preferred %in% colnames(ref@meta.data)], colnames(ref@meta.data)))
  for (cn in candidates) {
    vv <- ref@meta.data[[cn]]
    if (valid_cluster_vector(vv)) {
      message("Using reference cluster metadata column: ", cn)
      return(normalize_cluster_label(vv))
    }
  }
  stop("Could not safely identify reference labels 1,2,3,4,5,7")
}

# =============================================================================
# A. 读取 V2 已完成 Harmony 的全细胞对象
# =============================================================================
if (!file.exists(FULL_RDS)) stop("Full V2 object not found: ", FULL_RDS)
full <- readRDS(FULL_RDS)
if (!"RNA" %in% Assays(full)) stop("Full object has no RNA assay.")
full <- join_layers_safe(full, "RNA")
DefaultAssay(full) <- "RNA"

required_meta <- c("sample_id", "patient_id", "condition")
if (!all(required_meta %in% colnames(full@meta.data))) {
  stop("Full object lacks required metadata: ", paste(setdiff(required_meta, colnames(full@meta.data)), collapse = ", "))
}
if (!"full_cluster" %in% colnames(full@meta.data)) {
  if (!"seurat_clusters" %in% colnames(full@meta.data)) stop("No cluster metadata found.")
  full$full_cluster <- as.character(full$seurat_clusters)
}
full$condition <- factor(as.character(full$condition), levels = c("PA", "AC"))

# 重新标准化 RNA data，确保后续评分一致；原始 counts 不改变。
full <- NormalizeData(full, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
full_counts <- get_assay_matrix(full, "RNA", "counts")
full_data <- get_assay_matrix(full, "RNA", "data")

# =============================================================================
# B. 分层筛选：筛选基因与待验证 C2/C4 程序分离，避免循环论证
# =============================================================================
marker_sets <- list(
  EC_selection_anchor = c(
    "PECAM1", "VWF", "CDH5", "KDR", "ESAM",
    "EMCN", "ENG", "PLVAP", "EGFL7"
  ),
  EC_support_only = c("CLDN5", "RAMP2", "TEK", "ERG", "EPAS1"),
  Immune = c(
    "PTPRC", "LST1", "TYROBP", "FCER1G", "LYZ",
    "CD3D", "CD3E", "MS4A1", "CD79A", "NKG7"
  ),
  SMC_Pericyte_generic = c(
    "ACTA2", "TAGLN", "MYH11", "CNN1", "RGS5",
    "CSPG4", "MCAM", "PDGFRB", "DES"
  ),
  Fibroblast_generic = c(
    "COL1A1", "COL1A2", "DCN", "LUM", "COL3A1",
    "PDGFRA", "PI16", "CCDC80", "COL6A1"
  )
)

full$n_EC_anchor <- marker_count(full_counts, marker_sets$EC_selection_anchor)
full$n_EC_support <- marker_count(full_counts, marker_sets$EC_support_only)
full$n_immune <- marker_count(full_counts, marker_sets$Immune)
full$n_SMC_generic <- marker_count(full_counts, marker_sets$SMC_Pericyte_generic)
full$n_fibro_generic <- marker_count(full_counts, marker_sets$Fibroblast_generic)
full$PTPRC_detected <- if ("PTPRC" %in% rownames(full_counts)) as.numeric(full_counts["PTPRC", ] > 0) else 0

full$EC_anchor_z <- gene_z_score(full_data, marker_sets$EC_selection_anchor)
full$Immune_z_v3 <- gene_z_score(full_data, marker_sets$Immune)
full$SMC_z_v3 <- gene_z_score(full_data, marker_sets$SMC_Pericyte_generic)
full$Fibro_z_v3 <- gene_z_score(full_data, marker_sets$Fibroblast_generic)

cluster_audit <- full@meta.data %>%
  group_by(full_cluster) %>%
  summarise(
    n_cells = n(),
    frac_anchor_2plus = mean(n_EC_anchor >= 2),
    frac_anchor_3plus = mean(n_EC_anchor >= 3),
    frac_PTPRC_positive = mean(PTPRC_detected > 0),
    frac_immune_2plus = mean(n_immune >= 2),
    median_EC_anchor_z = median(EC_anchor_z, na.rm = TRUE),
    median_Immune_z = median(Immune_z_v3, na.rm = TRUE),
    median_SMC_z = median(SMC_z_v3, na.rm = TRUE),
    median_Fibro_z = median(Fibro_z_v3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    candidate_EC_cluster =
      frac_anchor_2plus >= 0.15 &
      frac_PTPRC_positive <= 0.20 &
      frac_immune_2plus <= 0.25 &
      median_EC_anchor_z >= median_Immune_z + 0.05
  )
write.csv(cluster_audit, file.path(OUT_DIR, "01_cluster_audit_LAYERED_V3.csv"), row.names = FALSE)

candidate_clusters <- cluster_audit$full_cluster[cluster_audit$candidate_EC_cluster]
if (length(candidate_clusters) == 0) stop("No candidate EC clusters passed V3 audit.")
message("V3 candidate EC clusters: ", paste(candidate_clusters, collapse = ", "))

max_nonEC_z <- pmax(full$Immune_z_v3, full$SMC_z_v3, full$Fibro_z_v3, na.rm = TRUE)
max_mesenchymal_z <- pmax(full$SMC_z_v3, full$Fibro_z_v3, na.rm = TRUE)

# Tier 1：高特异性经典 EC。
classical <-
  full$full_cluster %in% candidate_clusters &
  full$n_EC_anchor >= 3 &
  full$PTPRC_detected == 0 &
  full$n_immune <= 1 &
  full$n_SMC_generic <= 1 &
  full$n_fibro_generic <= 1 &
  full$EC_anchor_z >= max_nonEC_z + 0.10

# Tier 2 总入口：仍保留至少 2 个内皮谱系锚点，严格排除免疫，
# 但允许其间质/平滑肌分数略高于经典 EC。
lineage_supported <-
  full$full_cluster %in% candidate_clusters &
  full$n_EC_anchor >= 2 &
  full$PTPRC_detected == 0 &
  full$n_immune <= 1 &
  full$EC_anchor_z >= full$Immune_z_v3 + 0.10 &
  full$EC_anchor_z >= max_mesenchymal_z - 0.75

transition_like <-
  lineage_supported & !classical &
  (
    full$n_SMC_generic >= 1 |
    full$n_fibro_generic >= 1 |
    max_mesenchymal_z > -0.10
  )

low_marker_ec <- lineage_supported & !classical & !transition_like

full$EC_layer_v3 <- "Excluded_nonEC"
full$EC_layer_v3[classical] <- "Tier1_classical_EC"
full$EC_layer_v3[transition_like] <- "Tier2_transition_candidate"
full$EC_layer_v3[low_marker_ec] <- "Tier2_low_marker_EC"
full$EC_layer_v3 <- factor(
  full$EC_layer_v3,
  levels = c(
    "Tier1_classical_EC",
    "Tier2_transition_candidate",
    "Tier2_low_marker_EC",
    "Excluded_nonEC"
  )
)
full$expanded_EC_v3 <- full$EC_layer_v3 != "Excluded_nonEC"

layer_counts <- full@meta.data %>%
  count(sample_id, patient_id, condition, EC_layer_v3, name = "n_cells") %>%
  tidyr::complete(
    sample_id,
    EC_layer_v3 = levels(full$EC_layer_v3),
    fill = list(n_cells = 0)
  ) %>%
  left_join(
    full@meta.data %>% distinct(sample_id, patient_id, condition),
    by = "sample_id",
    suffix = c("", ".meta")
  ) %>%
  mutate(
    patient_id = coalesce(patient_id, patient_id.meta),
    condition = coalesce(as.character(condition), as.character(condition.meta))
  ) %>%
  select(sample_id, patient_id, condition, EC_layer_v3, n_cells)
write.csv(layer_counts, file.path(OUT_DIR, "02_EC_layer_counts_by_sample_V3.csv"), row.names = FALSE)

expanded_counts <- full@meta.data[full$expanded_EC_v3, , drop = FALSE] %>%
  count(sample_id, patient_id, condition, name = "n_expanded_EC") %>%
  mutate(sample_eligible_30 = n_expanded_EC >= 30)
write.csv(expanded_counts, file.path(OUT_DIR, "03_expanded_EC_counts_and_eligibility_V3.csv"), row.names = FALSE)

all_samples <- sort(unique(full$sample_id))
missing_expanded <- setdiff(all_samples, expanded_counts$sample_id)
if (length(missing_expanded) > 0) stop("No expanded EC retained in: ", paste(missing_expanded, collapse = ", "))
if (sum(full$expanded_EC_v3) < 500) stop("Expanded EC retained fewer than 500 cells.")

pair_eligibility <- expanded_counts %>%
  select(patient_id, condition, n_expanded_EC, sample_eligible_30) %>%
  pivot_wider(
    names_from = condition,
    values_from = c(n_expanded_EC, sample_eligible_30),
    values_fill = list(n_expanded_EC = 0, sample_eligible_30 = FALSE)
  ) %>%
  mutate(pair_eligible_30 = sample_eligible_30_PA & sample_eligible_30_AC)
write.csv(pair_eligibility, file.path(OUT_DIR, "04_patient_pair_eligibility_V3.csv"), row.names = FALSE)
eligible_patients <- pair_eligibility$patient_id[pair_eligibility$pair_eligible_30]

p_layer <- DimPlot(full, reduction = "umap", group.by = "EC_layer_v3", pt.size = 0.12) +
  ggtitle("V3 layered EC selection on full GSE159677 atlas")
save_png(p_layer, "01_full_UMAP_EC_layers_V3.png", 11, 8)

p_layer_split <- DimPlot(
  full, reduction = "umap", group.by = "EC_layer_v3", split.by = "condition", pt.size = 0.10
) + ggtitle("V3 layered EC selection by region")
save_png(p_layer_split, "02_full_UMAP_EC_layers_by_condition_V3.png", 15, 7)

selection_panel <- intersect(unique(c(
  marker_sets$EC_selection_anchor,
  marker_sets$EC_support_only,
  marker_sets$Immune,
  marker_sets$SMC_Pericyte_generic,
  marker_sets$Fibroblast_generic
)), rownames(full))
Idents(full) <- full$EC_layer_v3
p_layer_dot <- DotPlot(full, features = selection_panel, group.by = "EC_layer_v3") +
  RotatedAxis() + ggtitle("Lineage evidence and non-EC programs across V3 layers")
save_png(p_layer_dot, "03_EC_layer_marker_dotplot_V3.png", 18, 7)

# =============================================================================
# C. Expanded EC 重新降维、Harmony、聚类
# =============================================================================
expanded_cells <- rownames(full@meta.data)[full$expanded_EC_v3]
ec <- subset(full, cells = expanded_cells)
ec <- join_layers_safe(ec, "RNA")
DefaultAssay(ec) <- "RNA"
ec <- NormalizeData(ec, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
ec <- FindVariableFeatures(ec, selection.method = "vst", nfeatures = min(2500, nrow(ec)), verbose = FALSE)
ec <- ScaleData(ec, features = VariableFeatures(ec), verbose = FALSE)
ec <- RunPCA(ec, features = VariableFeatures(ec), npcs = 40, verbose = FALSE)
ec <- run_harmony_strict(ec, group_var = "patient_id", reduction_name = "pca", save_name = "harmony")

dims_ec <- seq_len(min(30, ncol(Embeddings(ec, "harmony"))))
ec <- RunUMAP(ec, reduction = "harmony", dims = dims_ec, seed.use = 20260731, verbose = FALSE)
ec <- FindNeighbors(ec, reduction = "harmony", dims = dims_ec, verbose = FALSE)
ec <- FindClusters(ec, resolution = 0.4, random.seed = 20260731, verbose = FALSE)
ec$expanded_EC_cluster <- as.character(ec$seurat_clusters)

p_ec_condition <- DimPlot(ec, group.by = "condition", pt.size = 0.25) +
  ggtitle("Expanded EC: AC versus PA")
p_ec_sample <- DimPlot(ec, group.by = "sample_id", pt.size = 0.22) +
  ggtitle("Expanded EC: Harmony corrected by patient")
p_ec_layer <- DimPlot(ec, group.by = "EC_layer_v3", pt.size = 0.25) +
  ggtitle("Expanded EC: Tier 1 and Tier 2 distribution")
p_ec_cluster <- DimPlot(ec, group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 0.22) +
  ggtitle("Expanded EC reclustering")
save_png(p_ec_condition, "04_expanded_EC_UMAP_condition_V3.png", 9, 7)
save_png(p_ec_sample, "05_expanded_EC_UMAP_sample_Harmony_V3.png", 10, 7)
save_png(p_ec_layer, "06_expanded_EC_UMAP_layer_V3.png", 10, 7)
save_png(p_ec_cluster, "07_expanded_EC_UMAP_clusters_V3.png", 9, 7)

# =============================================================================
# D. 固定程序评分：不用于细胞筛选
# =============================================================================
signatures <- list(
  C7_arterial_identity = c("EFNB2", "GJA4", "GJA5", "SOX17", "HEY1", "SEMA3G", "CXCL12"),
  C7_homeostatic_protection = c("KLF2", "KLF4", "NOS3", "THBD", "TEK", "ERG", "CLDN5", "RAMP2", "EPAS1"),
  C3_ANGPT2_PGF_activation = c("ANGPT2", "PGF"),
  C1_C5_inflammatory_adhesion = c("VCAM1", "SELE", "ICAM1", "CXCL2", "CXCL3", "NFKBIA", "IFITM3", "IFIT1", "ISG15", "STAT1"),
  C2_vascular_remodeling = c("ITLN1", "SULF1", "MGP", "CLU", "ELN", "COL8A1", "COL8A2", "EFEMP1", "MFAP4"),
  C4_ECM_remodeling_primary = c("POSTN", "IGFBP5", "SPARC", "COL6A2", "TIMP1", "MAFB", "NR2F2", "FN1", "VCAN"),
  C4_broad_ECM_sensitivity = c("POSTN", "IGFBP5", "SPARC", "COL6A2", "TIMP1", "MAFB", "NR2F2", "FN1", "VCAN", "COL1A1", "COL1A2", "COL3A1"),
  Generic_fibro_QC = marker_sets$Fibroblast_generic,
  Generic_SMC_QC = marker_sets$SMC_Pericyte_generic
)

expected_program <- c(
  C7_arterial_identity = NA_character_,
  C7_homeostatic_protection = "PA>AC",
  C3_ANGPT2_PGF_activation = "AC>PA",
  C1_C5_inflammatory_adhesion = "AC>PA",
  C2_vascular_remodeling = "AC>PA",
  C4_ECM_remodeling_primary = "AC>PA",
  C4_broad_ECM_sensitivity = "AC>PA",
  Generic_fibro_QC = NA_character_,
  Generic_SMC_QC = NA_character_
)

signatures_available <- lapply(signatures, function(x) unique(intersect(x, rownames(ec))))
sig_table <- bind_rows(lapply(names(signatures_available), function(nm) {
  data.frame(
    program = nm,
    expected_direction = expected_program[[nm]],
    gene = signatures_available[[nm]]
  )
}))
write.csv(sig_table, file.path(OUT_DIR, "05_fixed_validation_signatures_V3.csv"), row.names = FALSE)

ec_data <- get_assay_matrix(ec, "RNA", "data")
for (nm in names(signatures_available)) {
  ec[[paste0(nm, "_meanlog")]] <- mean_log_score(ec_data, signatures_available[[nm]])
  ec[[paste0(nm, "_z")]] <- gene_z_score(ec_data, signatures_available[[nm]])
}

# =============================================================================
# E. 参考标签转移
# =============================================================================
if (!file.exists(REF_RDS)) stop("Reference object not found: ", REF_RDS)
ref <- readRDS(REF_RDS)
if (!"RNA" %in% Assays(ref)) stop("Reference object has no RNA assay.")
ref <- join_layers_safe(ref, "RNA")
DefaultAssay(ref) <- "RNA"
ref$reference_cluster_num <- find_reference_cluster_labels(ref)
ref <- subset(ref, subset = reference_cluster_num %in% c("1", "2", "3", "4", "5", "7"))
ref$reference_state <- paste0("C", ref$reference_cluster_num)
Idents(ref) <- ref$reference_state
ref <- NormalizeData(ref, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = min(3000, nrow(ref)), verbose = FALSE)
ref <- ScaleData(ref, features = VariableFeatures(ref), verbose = FALSE)
ref <- RunPCA(ref, features = VariableFeatures(ref), npcs = 30, verbose = FALSE)

common_features <- intersect(VariableFeatures(ref), rownames(ec))
if (length(common_features) < 500) stop("Too few common variable genes for mapping: ", length(common_features))
ndim_transfer <- min(30, ncol(Embeddings(ref, "pca")), ncol(Embeddings(ec, "pca")))
anchors <- FindTransferAnchors(
  reference = ref,
  query = ec,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  features = common_features,
  dims = seq_len(ndim_transfer),
  verbose = FALSE
)
pred <- TransferData(
  anchorset = anchors,
  refdata = ref$reference_state,
  dims = seq_len(ndim_transfer),
  verbose = FALSE
)
score_cols <- grep("^prediction\\.score\\.C", colnames(pred), value = TRUE)
if (length(score_cols) >= 2) {
  sm <- as.matrix(pred[, score_cols, drop = FALSE])
  pred$prediction.margin <- apply(sm, 1, function(x) {
    sx <- sort(x, decreasing = TRUE)
    sx[1] - sx[2]
  })
} else {
  pred$prediction.margin <- NA_real_
}
ec <- AddMetaData(ec, pred)
ec$mapped_state <- ifelse(
  ec$prediction.score.max >= 0.50 &
    (is.na(ec$prediction.margin) | ec$prediction.margin >= 0.05),
  as.character(ec$predicted.id),
  "Uncertain"
)

p_state <- DimPlot(ec, group.by = "mapped_state", pt.size = 0.24) +
  ggtitle("Expanded EC mapped to discovery states")
p_conf <- FeaturePlot(ec, features = "prediction.score.max", pt.size = 0.24) +
  ggtitle("Reference mapping confidence")
save_png(p_state, "08_expanded_EC_UMAP_mapped_state_V3.png", 10, 7)
save_png(p_conf, "09_expanded_EC_mapping_confidence_V3.png", 9, 7)

state_panel <- intersect(unique(c(
  marker_sets$EC_selection_anchor,
  signatures$C7_arterial_identity,
  signatures$C7_homeostatic_protection,
  signatures$C3_ANGPT2_PGF_activation,
  "CCL14",
  signatures$C1_C5_inflammatory_adhesion,
  signatures$C2_vascular_remodeling,
  signatures$C4_ECM_remodeling_primary,
  marker_sets$SMC_Pericyte_generic,
  marker_sets$Fibroblast_generic
)), rownames(ec))
p_state_dot <- DotPlot(ec, features = state_panel, group.by = "mapped_state") +
  RotatedAxis() + ggtitle("Fixed discovery programs in mapped expanded EC states")
save_png(p_state_dot, "10_expanded_EC_state_marker_dotplot_V3.png", 20, 9)

layer_state_counts <- ec@meta.data %>%
  count(sample_id, patient_id, condition, EC_layer_v3, mapped_state, name = "n_cells") %>%
  group_by(sample_id, EC_layer_v3) %>%
  mutate(proportion_within_layer = n_cells / sum(n_cells)) %>%
  ungroup()
write.csv(layer_state_counts, file.path(OUT_DIR, "06_mapped_state_counts_by_layer_V3.csv"), row.names = FALSE)

# =============================================================================
# F. 样本级程序评分：Expanded 主分析 + Tier 1/Tier 2 敏感性
# =============================================================================
program_names <- names(signatures_available)
program_cols <- paste0(program_names, "_z")

sample_program_layer <- ec@meta.data %>%
  group_by(sample_id, patient_id, condition, EC_layer_v3) %>%
  summarise(
    n_cells = n(),
    across(all_of(program_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )
colnames(sample_program_layer)[match(program_cols, colnames(sample_program_layer))] <- program_names
write.csv(sample_program_layer, file.path(OUT_DIR, "07_sample_program_scores_by_layer_V3.csv"), row.names = FALSE)

sample_program_expanded <- ec@meta.data %>%
  group_by(sample_id, patient_id, condition) %>%
  summarise(
    n_expanded_EC = n(),
    across(all_of(program_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )
colnames(sample_program_expanded)[match(program_cols, colnames(sample_program_expanded))] <- program_names
write.csv(sample_program_expanded, file.path(OUT_DIR, "08_sample_program_scores_expanded_V3.csv"), row.names = FALSE)

expanded_long <- sample_program_expanded %>%
  pivot_longer(cols = all_of(program_names), names_to = "program", values_to = "mean_z_score")
layer_long <- sample_program_layer %>%
  pivot_longer(cols = all_of(program_names), names_to = "program", values_to = "mean_z_score")

p_expanded_program <- ggplot(
  expanded_long,
  aes(condition, mean_z_score, group = patient_id, colour = patient_id)
) +
  geom_line(linewidth = 0.75, alpha = 0.85) +
  geom_point(size = 2.6) +
  facet_wrap(~ program, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "Mean standardized score", title = "Expanded EC paired validation: PA versus AC")
save_png(p_expanded_program, "11_paired_program_scores_expanded_V3.png", 13, 14)

p_layer_program <- ggplot(
  layer_long %>% filter(program %in% c(
    "C7_homeostatic_protection", "C3_ANGPT2_PGF_activation",
    "C1_C5_inflammatory_adhesion", "C2_vascular_remodeling",
    "C4_ECM_remodeling_primary"
  )),
  aes(condition, mean_z_score, group = patient_id, colour = patient_id)
) +
  geom_line(linewidth = 0.65, alpha = 0.8) +
  geom_point(size = 2.2) +
  facet_grid(EC_layer_v3 ~ program, scales = "free_y") +
  labs(x = NULL, y = "Mean standardized score", title = "Program direction across EC layers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_png(p_layer_program, "12_program_scores_by_layer_V3.png", 19, 10)

stats_all <- paired_stats(expanded_long, "mean_z_score", "program", expected_program)
stats_eligible <- paired_stats(
  expanded_long, "mean_z_score", "program", expected_program,
  eligible_patients = eligible_patients
)
write.csv(stats_all, file.path(OUT_DIR, "09_paired_program_statistics_all_pairs_V3.csv"), row.names = FALSE)
write.csv(stats_eligible, file.path(OUT_DIR, "10_paired_program_statistics_eligible_pairs_V3.csv"), row.names = FALSE)

# =============================================================================
# G. 关键基因样本级表达
# =============================================================================
key_genes <- intersect(unique(c(
  signatures$C7_arterial_identity,
  signatures$C7_homeostatic_protection,
  signatures$C3_ANGPT2_PGF_activation,
  "CCL14",
  signatures$C1_C5_inflammatory_adhesion,
  signatures$C2_vascular_remodeling,
  signatures$C4_ECM_remodeling_primary
)), rownames(ec))

key_rows <- lapply(unique(ec$sample_id), function(sid) {
  cells <- colnames(ec)[ec$sample_id == sid]
  md <- ec@meta.data[cells[1], c("sample_id", "patient_id", "condition"), drop = FALSE]
  vals <- Matrix::rowMeans(ec_data[key_genes, cells, drop = FALSE])
  data.frame(
    sample_id = as.character(md$sample_id),
    patient_id = as.character(md$patient_id),
    condition = as.character(md$condition),
    gene = names(vals),
    mean_log_normalized_expression = as.numeric(vals)
  )
})
key_gene_sample <- bind_rows(key_rows)
write.csv(key_gene_sample, file.path(OUT_DIR, "11_key_gene_sample_expression_expanded_V3.csv"), row.names = FALSE)

key_expected <- setNames(rep(NA_character_, length(key_genes)), key_genes)
key_expected[intersect(names(key_expected), signatures$C7_homeostatic_protection)] <- "PA>AC"
key_expected[intersect(names(key_expected), unique(c(
  signatures$C3_ANGPT2_PGF_activation,
  signatures$C1_C5_inflammatory_adhesion,
  signatures$C2_vascular_remodeling,
  signatures$C4_ECM_remodeling_primary
)))] <- "AC>PA"
key_stats_all <- paired_stats(
  key_gene_sample,
  "mean_log_normalized_expression",
  "gene",
  key_expected
)
key_stats_eligible <- paired_stats(
  key_gene_sample,
  "mean_log_normalized_expression",
  "gene",
  key_expected,
  eligible_patients = eligible_patients
)
write.csv(key_stats_all, file.path(OUT_DIR, "12_key_gene_paired_statistics_all_pairs_V3.csv"), row.names = FALSE)
write.csv(key_stats_eligible, file.path(OUT_DIR, "13_key_gene_paired_statistics_eligible_pairs_V3.csv"), row.names = FALSE)

p_key <- ggplot(
  key_gene_sample,
  aes(condition, mean_log_normalized_expression, group = patient_id, colour = patient_id)
) +
  geom_line(linewidth = 0.55, alpha = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  labs(x = NULL, y = "Mean log-normalized expression", title = "Expanded EC key genes: paired PA versus AC")
save_png(p_key, "13_key_gene_paired_expression_expanded_V3.png", 17, 15)

# =============================================================================
# H. Tier 2 过渡候选的组成与 C4/C2 稳健性审计
# =============================================================================
tier2_summary <- ec@meta.data %>%
  mutate(
    is_tier2_transition = EC_layer_v3 == "Tier2_transition_candidate",
    is_tier1 = EC_layer_v3 == "Tier1_classical_EC"
  ) %>%
  group_by(sample_id, patient_id, condition) %>%
  summarise(
    n_expanded_EC = n(),
    n_tier1 = sum(is_tier1),
    n_tier2_transition = sum(is_tier2_transition),
    transition_fraction = mean(is_tier2_transition),
    .groups = "drop"
  )
write.csv(tier2_summary, file.path(OUT_DIR, "14_transition_candidate_fraction_by_sample_V3.csv"), row.names = FALSE)

p_transition_fraction <- ggplot(
  tier2_summary,
  aes(condition, transition_fraction, group = patient_id, colour = patient_id)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  labs(x = NULL, y = "Tier 2 transition-candidate fraction", title = "Transition-candidate fraction among expanded EC")
save_png(p_transition_fraction, "14_transition_candidate_fraction_V3.png", 8, 6)

c4_cor <- data.frame(
  comparison = c(
    "C4_primary_vs_generic_fibro",
    "C4_primary_vs_generic_SMC",
    "C2_vs_generic_fibro",
    "C2_vs_generic_SMC"
  ),
  spearman_rho = c(
    suppressWarnings(cor(ec$C4_ECM_remodeling_primary_z, ec$Generic_fibro_QC_z, method = "spearman", use = "complete.obs")),
    suppressWarnings(cor(ec$C4_ECM_remodeling_primary_z, ec$Generic_SMC_QC_z, method = "spearman", use = "complete.obs")),
    suppressWarnings(cor(ec$C2_vascular_remodeling_z, ec$Generic_fibro_QC_z, method = "spearman", use = "complete.obs")),
    suppressWarnings(cor(ec$C2_vascular_remodeling_z, ec$Generic_SMC_QC_z, method = "spearman", use = "complete.obs"))
  )
)
write.csv(c4_cor, file.path(OUT_DIR, "15_remodeling_vs_nonEC_score_correlations_V3.csv"), row.names = FALSE)

# =============================================================================
# I. 映射状态比例
# =============================================================================
all_states <- c("C1", "C2", "C3", "C4", "C5", "C7", "Uncertain")
# 每个 sample_id 只保留一行样本信息。不要先在 state_counts 中保留
# patient_id/condition 后再重复连接，否则 dplyr 会生成 condition.x/condition.y，
# ggplot2 中 aes(condition, ...) 就会报 “object 'condition' not found”。
sample_meta_state <- ec@meta.data %>%
  transmute(
    sample_id = as.character(sample_id),
    patient_id = as.character(patient_id),
    condition = factor(as.character(condition), levels = c("PA", "AC"))
  ) %>%
  distinct(sample_id, .keep_all = TRUE)

if (anyDuplicated(sample_meta_state$sample_id)) {
  stop("sample_meta_state contains duplicated sample_id values.")
}
if (any(is.na(sample_meta_state$condition))) {
  stop("condition contains values other than PA/AC in sample metadata.")
}

state_counts <- ec@meta.data %>%
  transmute(
    sample_id = as.character(sample_id),
    mapped_state = as.character(mapped_state)
  ) %>%
  count(sample_id, mapped_state, name = "n_cells")

sample_totals <- ec@meta.data %>%
  transmute(sample_id = as.character(sample_id)) %>%
  count(sample_id, name = "n_total_expanded_EC")

state_diag <- tidyr::expand_grid(
  sample_id = sort(unique(as.character(ec$sample_id))),
  mapped_state = all_states
) %>%
  left_join(sample_meta_state, by = "sample_id") %>%
  left_join(state_counts, by = c("sample_id", "mapped_state")) %>%
  mutate(n_cells = tidyr::replace_na(n_cells, 0L)) %>%
  left_join(sample_totals, by = "sample_id") %>%
  mutate(
    condition = factor(as.character(condition), levels = c("PA", "AC")),
    proportion = ifelse(n_total_expanded_EC > 0, n_cells / n_total_expanded_EC, NA_real_)
  )

required_state_cols <- c(
  "sample_id", "patient_id", "condition", "mapped_state",
  "n_cells", "n_total_expanded_EC", "proportion"
)
if (!all(required_state_cols %in% colnames(state_diag))) {
  stop(
    "state_diag lacks required columns: ",
    paste(setdiff(required_state_cols, colnames(state_diag)), collapse = ", ")
  )
}
write.csv(state_diag, file.path(OUT_DIR, "16_mapped_state_counts_proportions_expanded_V3.csv"), row.names = FALSE)

p_state_pair <- ggplot(
  state_diag,
  aes(condition, proportion, group = patient_id, colour = patient_id)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  facet_wrap(~ mapped_state, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Proportion among expanded EC", title = "Mapped-state proportions in expanded EC")
save_png(p_state_pair, "15_paired_mapped_state_proportions_V3.png", 13, 9)

# =============================================================================
# J. Expanded EC pseudobulk；只有 3 对均达到 30 细胞才做正式 edgeR
# =============================================================================
ec_counts <- get_assay_matrix(ec, "RNA", "counts")
samples <- sort(unique(ec$sample_id))
pb_counts <- sapply(samples, function(sid) {
  cells <- colnames(ec)[ec$sample_id == sid]
  Matrix::rowSums(ec_counts[, cells, drop = FALSE])
})
if (!is.matrix(pb_counts)) pb_counts <- as.matrix(pb_counts)
colnames(pb_counts) <- samples
write.csv(
  data.frame(gene = rownames(pb_counts), pb_counts, check.names = FALSE),
  file.path(OBJ_DIR, "expanded_EC_pseudobulk_raw_counts_V3.csv"),
  row.names = FALSE
)

formal_edgeR <- FALSE
if (length(eligible_patients) == 3) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    writeLines(
      "All 3 pairs were cell-count eligible, but edgeR is not installed.",
      file.path(OUT_DIR, "17_edgeR_NOT_RUN_V3.txt")
    )
  } else {
    formal_edgeR <- tryCatch({
      suppressPackageStartupMessages(library(edgeR))
      pb_meta <- ec@meta.data %>%
        distinct(sample_id, patient_id, condition) %>%
        filter(sample_id %in% colnames(pb_counts)) %>%
        arrange(match(sample_id, colnames(pb_counts)))
      pb_meta$condition <- factor(as.character(pb_meta$condition), levels = c("PA", "AC"))
      pb_meta$patient_id <- factor(pb_meta$patient_id)
      design <- model.matrix(~ patient_id + condition, data = pb_meta)
      if (qr(design)$rank < ncol(design)) stop("Paired design matrix is not full rank.")
      y <- edgeR::DGEList(counts = pb_counts[, pb_meta$sample_id, drop = FALSE])
      keep <- edgeR::filterByExpr(y, design = design, min.count = 5)
      y <- y[keep, , keep.lib.sizes = FALSE]
      y <- edgeR::calcNormFactors(y)
      y <- edgeR::estimateDisp(y, design, robust = TRUE)
      fit <- edgeR::glmQLFit(y, design, robust = TRUE)
      coef_idx <- which(colnames(design) == "conditionAC")
      if (length(coef_idx) != 1) stop("conditionAC coefficient not found.")
      qlf <- edgeR::glmQLFTest(fit, coef = coef_idx)
      tab <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table
      tab$gene <- rownames(tab)
      tab <- tab[, c("gene", setdiff(colnames(tab), "gene"))]
      write.csv(tab, file.path(OBJ_DIR, "expanded_EC_edgeR_all_genes_V3.csv"), row.names = FALSE)
      write.csv(
        tab[tab$gene %in% key_genes, , drop = FALSE],
        file.path(OUT_DIR, "17_expanded_EC_edgeR_key_genes_V3.csv"),
        row.names = FALSE
      )
      TRUE
    }, error = function(e) {
      writeLines(
        paste("edgeR failed:", conditionMessage(e)),
        file.path(OUT_DIR, "17_edgeR_NOT_RUN_V3.txt")
      )
      FALSE
    })
  }
} else {
  writeLines(
    c(
      "Formal paired edgeR was intentionally skipped.",
      paste0("Eligible patient pairs with >=30 expanded EC in both PA and AC: ", length(eligible_patients), "/3."),
      "Pseudobulk raw counts were saved for transparency, but no inferential p-values are reported."
    ),
    file.path(OUT_DIR, "17_edgeR_NOT_RUN_V3.txt")
  )
}

# =============================================================================
# K. 自动证据等级：严格区分描述性与正式验证
# =============================================================================
primary_programs <- c(
  "C7_homeostatic_protection",
  "C3_ANGPT2_PGF_activation",
  "C1_C5_inflammatory_adhesion",
  "C2_vascular_remodeling",
  "C4_ECM_remodeling_primary"
)
summary_program <- stats_all %>%
  filter(feature %in% primary_programs) %>%
  mutate(
    n_eligible_pairs = length(eligible_patients),
    evidence_grade = case_when(
      length(eligible_patients) < 3 ~ "descriptive_only_incomplete_cell_support",
      n_pairs_in_expected_direction == 3 ~ "directionally_consistent_all_3_eligible_pairs",
      n_pairs_in_expected_direction == 2 ~ "partially_consistent_2_of_3_eligible_pairs",
      TRUE ~ "not_consistent"
    )
  ) %>%
  select(
    validation_item = feature,
    expected_direction,
    n_eligible_pairs,
    mean_AC_minus_PA,
    paired_Cohen_dz,
    n_pairs_in_expected_direction,
    fraction_in_expected_direction,
    evidence_grade
  )
write.csv(summary_program, file.path(OUT_DIR, "18_primary_validation_summary_LAYERED_V3.csv"), row.names = FALSE)

# 保存对象；不放入轻量包。
saveRDS(full, file.path(OBJ_DIR, "GSE159677_full_with_EC_layers_V3.rds"), compress = "gzip")
saveRDS(ec, file.path(OBJ_DIR, "GSE159677_expanded_EC_layered_validation_V3.rds"), compress = "gzip")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_LAYERED_V3.txt"))

summary_text <- c(
  "GSE159677 V3 layered EC validation completed (FIXED).",
  paste0("Full cells: ", ncol(full)),
  paste0("Tier 1 classical EC: ", sum(full$EC_layer_v3 == "Tier1_classical_EC")),
  paste0("Tier 2 transition candidates: ", sum(full$EC_layer_v3 == "Tier2_transition_candidate")),
  paste0("Tier 2 low-marker EC: ", sum(full$EC_layer_v3 == "Tier2_low_marker_EC")),
  paste0("Expanded EC total: ", ncol(ec)),
  paste0("Eligible patient pairs (both PA and AC >=30 expanded EC): ", length(eligible_patients), "/3"),
  paste0("Eligible patients: ", ifelse(length(eligible_patients) > 0, paste(eligible_patients, collapse = ", "), "none")),
  paste0("Harmony on expanded EC: ", "harmony" %in% Reductions(ec)),
  paste0("Formal paired edgeR completed: ", formal_edgeR),
  paste0("Reference mapping uncertain fraction: ", round(mean(ec$mapped_state == "Uncertain"), 4)),
  "Tier definitions were fixed without using C2 or C4 validation genes.",
  "Tier 1 is high-specificity; Tier 2 is a sensitivity layer and does not prove EndMT.",
  "Patient, not cell, is the biological replicate.",
  "PA is adjacent diseased artery, not a truly healthy control."
)
writeLines(summary_text, file.path(OUT_DIR, "README_LAYERED_V3_results.txt"))
cat(paste(summary_text, collapse = "\n"), "\n")
